
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

!> List of Parameters see Parameter list in mumpstools.hpp

module mumpstools_common
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only : &
            stdout=>output_unit, &
            stderr=>error_unit
    implicit none

    include 'mpif.h'
    include 'dmumps_struc.h'

#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

    type ( DMUMPS_STRUC ), dimension( : ), allocatable, target, save :: gSolvers
    integer, dimension( : ), allocatable, save :: gOccupied
    ! FIXED at 8. The setter that used to change it was removed 2026-08-29:
    ! it had no caller anywhere, it did not reallocate the pool, and a call
    ! after the first create would have left gMaxNumSolvers and the actual
    ! array size disagreeing -- a scan past the end, or a deallocate with
    ! high slots still occupied. Raising the pool size is a reallocation
    ! feature; add it as one if it is ever wanted, not as a bare setter
    integer( int_t ), save :: gMaxNumSolvers = 8
    integer( int_t ), save :: gNumSolvers = 0

end module mumpstools_common

!-------------------------------------------------------------------------------
!> Creates a new MUMPS solver instance
!> Initializes the solver structure and reserves a slot in the global solver array
!> @param[out] aSolverID Unique ID of the created solver (1-based), -1 on failure
!> @param[out] aInfo     Info status (0 = success, -1000 = max solvers reached)
!> @param[in]  aHostIsWorking Passed through to the MUMPS control parameter PAR
!> @param[in]  aSymmetryMode  Passed through to the MUMPS control parameter SYM
!-------------------------------------------------------------------------------

subroutine mumpstools_create_solver( aSolverID, aInfo, aHostIsWorking, aSymmetryMode ) bind( c )
    use mumpstools_common
    implicit none
    integer( int_t ), intent( out ) :: aSolverID
    integer( int_t ), intent( out ) :: aInfo
    integer( int_t ), intent( in )  :: aHostIsWorking
    integer( int_t ), intent( in )  :: aSymmetryMode

    type ( DMUMPS_STRUC ), pointer :: tMUMPS

    integer :: k

    ! MPI_BARRIER's status argument is a DEFAULT integer. aInfo is
    ! integer( int_t ), which is 64-bit under BELFEM_INT64, so passing it
    ! there was a kind mismatch -- and its value was overwritten by the
    ! INFO(1) copy two lines later anyway, so the barrier status was never
    ! read. mumpstools_free_solver already does this correctly
    integer :: ierr

    if ( .not. allocated( gSolvers ) .and. gMaxNumSolvers .gt. 0 ) then
        allocate( gSolvers( gMaxNumSolvers ) )
        allocate( gOccupied( gMaxNumSolvers ) )
        forall( k=1:gMaxNumSolvers  ) gOccupied( k ) = 0
        gNumSolvers = 0
    end if

    ! first-free scan: the occupancy table is the authority, not a
    ! high-water counter. A freed slot is reused, so a consumer that
    ! creates and frees repeatedly while another instance stays alive
    ! cannot exhaust the pool.
    aSolverID = -1
    do k = 1, gMaxNumSolvers
        if ( gOccupied( k ) .eq. 0 ) then
            aSolverID = k
            exit
        end if
    end do

    if ( aSolverID .lt. 0 ) then
        aInfo = -1000

        ! the exhaustion arm is collective too. Both arms of an agreeing pool
        ! now carry a barrier, so a create no longer has one outcome that
        ! synchronizes and one that does not.
        !
        ! This is hygiene, NOT a fix for ranks that DISAGREE about which slot
        ! is free. On disagreement one rank returns here while another calls
        ! DMUMPS, which is itself collective on MPI_COMM_WORLD -- the hang
        ! simply moves inside MUMPS. Making that a loud error needs the chosen
        ! slot allreduced before DMUMPS; it is deliberately not done here,
        ! because a barrier placement is not that fix and must not be sold as
        ! one. No in-tree caller can drive the pool to mixed occupancy today
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    else
        gNumSolvers = gNumSolvers + 1
        tMUMPS => gSolvers( aSolverID )
        gOccupied( aSolverID ) = 1

        tMUMPS%COMM = MPI_COMM_WORLD
        tMUMPS%PAR = aHostIsWorking
        tMUMPS%SYM = aSymmetryMode
        tMUMPS%JOB = -1

        ! kept immediately before DMUMPS, where it belongs: it synchronizes
        ! the ranks for the collective JOB = -1, not for the slot decision
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        call DMUMPS( tMUMPS )

        aInfo = tMUMPS%INFO( 1 )

        ! roll the reservation back if the instance was never built. The
        ! occupancy and the counter are written BEFORE the call above, so a
        ! failed JOB = -1 used to leave the slot marked occupied forever --
        ! one lost slot per failure, out of eight. Returning -1 also makes
        ! this subroutine's own header comment true, and lets the existing
        ! C++ mSolverID <= 0 gate catch the case with no second code path.
        !
        ! The pointers are nullified as the free path does: the slot may be
        ! handed to the next tenant without the pool being reallocated, and
        ! its JOB = -1 must not see this attempt's user pointers.
        !
        ! NOT ATTEMPTED: a JOB = -2 on the failed instance. Whether MUMPS
        ! wants one after a failed JOB = -1 is not settled here -- the 5.5.1
        ! guide was not available to this session. Occupancy rollback is the
        ! right POOL accounting either way; a MUMPS-internal leak, if there
        ! is one, is invisible from this side
        if ( tMUMPS%INFO( 1 ) .lt. 0 ) then
            gOccupied( aSolverID ) = 0
            gNumSolvers = gNumSolvers - 1

            nullify( tMUMPS%irn )
            nullify( tMUMPS%jcn )
            nullify( tMUMPS%A )
            nullify( tMUMPS%rhs )

            aSolverID = -1
        end if
    end if

end subroutine mumpstools_create_solver

!-------------------------------------------------------------------------------
!> Frees a MUMPS solver instance and releases its resources
!> Calls MUMPS with JOB=-2 to destroy the solver instance
!> If all solvers are freed, deallocates the global solver arrays
!> @param[in]  aSolverID ID of the solver to free
!> @param[out] aInfo     Info status from MUMPS
!-------------------------------------------------------------------------------

subroutine mumpstools_free_solver( aSolverID, aInfo ) bind( c )
    use mumpstools_common
    implicit none
    integer( int_t ), intent( in ) :: aSolverID
    integer( int_t ), intent( out ), dimension( 80 ) :: aInfo

    integer :: k, count
    integer :: ierr

    type ( DMUMPS_STRUC ), pointer :: tMUMPS

    ! default to "no error" so the full INFO array is well-defined even when
    ! the instance was already freed
    aInfo = 0

    ! the full drain below deallocates the pool arrays, so a free on a stale
    ! ID after a drain must not index them -- guards the whole body, because
    ! the recount further down indexes gOccupied again
    if( .not. allocated( gOccupied ) ) then
        return
    end if

    ! ... and only then may the id be range-checked, because size() on an
    ! unallocated allocatable is undefined. Order matters here and the guard
    ! above must stay first.
    !
    ! size( gOccupied ) is the bound, NOT gMaxNumSolvers. They agree today,
    ! now that the pool size cannot be changed after allocation, but the
    ! array is what is actually indexed and the check should stay right if a
    ! reallocation path is ever added
    if( aSolverID .lt. 1 .or. aSolverID .gt. size( gOccupied ) ) then
        return
    end if

    if( gOccupied( aSolverID ) .ne. 0 ) then
        tMUMPS => gSolvers( aSolverID )
        tMUMPS%JOB = -2

        call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        call DMUMPS( tMUMPS )
        gOccupied( aSolverID ) = 0

        ! return the full INFO(1:80) so the C++ side can decode INFO(2) details
        forall( k = 1:80 ) aInfo( k ) = tMUMPS%INFO( k )

        ! the slot may be recycled without the full-drain reallocation below,
        ! so the next tenant re-initializes a used struct. JOB = -1 must not
        ! see the previous tenant's user pointers
        nullify( tMUMPS%irn )
        nullify( tMUMPS%jcn )
        nullify( tMUMPS%A )
        nullify( tMUMPS%rhs )
    end if

    count = 0
    do k = 1, gMaxNumSolvers
        if ( gOccupied( k ) .eq. 1 ) then
            count = count + 1
        end if
    end do

    ! occupancy count, kept in step with the table on every free -- not a
    ! high-water mark that only resets on full drain
    gNumSolvers = count

    if( count .eq. 0 ) then
        deallocate( gSolvers )
        deallocate( gOccupied )
    end if

end subroutine mumpstools_free_solver

!-------------------------------------------------------------------------------
!> Number of live solver instances ( pool occupancy ) on this process
!> @param[out] aCount occupancy; identical on every rank because create and
!>                    free are collective
!-------------------------------------------------------------------------------

subroutine mumpstools_num_solvers( aCount ) bind( c )
    use mumpstools_common
    implicit none
    integer( int_t ), intent( out ) :: aCount

    aCount = gNumSolvers

end subroutine mumpstools_num_solvers

!-------------------------------------------------------------------------------
!> Solves a linear system using MUMPS
!> Supports analysis, factorization, and solution phases based on JOB parameter
!> JOB=1: analysis, JOB=2: factorization, JOB=3: solve, JOB=5: factorization+solve with same structure
!> JOB=6: analysis+factorization+solve (all-in-one)
!-------------------------------------------------------------------------------

subroutine mumpstools_solve( &
        aIParameters, & !> list of integer parameters, see Parameter enum in mumpstools.hpp
        aRParameters, & !> list of real parameters (compression tolerance, etc.)
        aN,          & !> matrix dimension
        aNNZ,        & !> number of nonzeros
        aNRHS,       & !> number of columns on right hand side
        aRowIndices, & !> row indices of matrix (1-based, COO format)
        aColIndices, & !> column indices of matrix (1-based, COO format)
        aValues,     & !> matrix values in COO format
        aX,          & !> solution vector (output)
        aY,          & !> right hand side vector (input)
        aInfo,       & !> MUMPS INFO array (80 elements) for error/diagnostic information
        aInfoG,      & !> MUMPS INFOG array (80 elements) - rank-uniform error/diagnostic information
        aRInfoG      & !> MUMPS RINFOG array (20 elements) for real diagnostic information
        ) bind( c )
    
    use mumpstools_common
    implicit none
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! ARGUMENTS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    integer( int_t ), intent( in    ), dimension( 14 )                 :: aIParameters
    real( c_double ), intent( in ),    dimension( 1 )                  :: aRParameters
    integer( int_t ), intent( in    )                                :: aN
    integer( int_t ), intent( in    )                                :: aNNZ
    integer( int_t ), intent( in    )                                :: aNRHS
    integer( int_t ), intent( in    ), dimension( aNNZ ), target     :: aRowIndices
    integer( int_t ), intent( in    ), dimension( aNNZ ), target     :: aColIndices
    real( c_double ), intent( in    ), dimension( aNNZ ), target       :: aValues
    real( c_double ), intent( inout ), dimension( aN * aNRHS ), target :: aX
    real( c_double ), intent( in    ), dimension( aN * aNRHS ), target :: aY
    integer( int_t ), intent( out   ), dimension( 80  )                :: aInfo
    integer( int_t ), intent( out   ), dimension( 80  )                :: aInfoG
    real( c_double ), intent(out ), dimension( 20  )                   :: aRInfoG

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Shortcuts for better code readability of user parameters
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ) :: tSolverID
    integer( int_t ) :: tJob
    integer( int_t ) :: tMasterRank
    integer( int_t ) :: tHostIsWorking
    integer( int_t ) :: tSymmetryMode
    integer( int_t ) :: tInfoLevel
    integer( int_t ) :: tSerialPermutationOrdering
    integer( int_t ) :: tParallelPermutationOrdering
    integer( int_t ) :: tNumRefinementSteps
    integer( int_t ) :: tComputeDeterminant
    integer( int_t ) :: tMemoryRelaxation
    integer( int_t ) :: tMemoryBudget
    integer( int_t ) :: tCompressionMode
    integer( int_t ) :: tErrorAnalysisMode

    double precision :: tEpsilon

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! MUMPS stuff
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    !> the parameters object
    type ( DMUMPS_STRUC ), pointer :: tMUMPS

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Other parameters
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! iteration index
    integer( int_t ) :: k

    ! length of rhs memory
    integer( int_t ) :: tCapacity

    ! MPI status
    integer( int_t ) :: tStatus

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! POPULATE DEFINED PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    tSolverID                    = aIParameters(  1 )
    tJob                         = aIParameters(  2 )

    tMasterRank                  = aIParameters(  3 )
    tInfoLevel                   = aIParameters(  4 )
    tErrorAnalysisMode           = aIParameters(  5 )
    
    tHostIsWorking               = aIParameters(  6 )
    tSymmetryMode                = aIParameters(  7 )
   
    tSerialPermutationOrdering   = aIParameters(  8 )
    tParallelPermutationOrdering = aIParameters(  9 )
    tCompressionMode             = aIParameters( 10 )
    tNumRefinementSteps          = aIParameters( 11 )
    tMemoryRelaxation            = aIParameters( 12 )
    tComputeDeterminant          = aIParameters( 13 )
    tMemoryBudget                = aIParameters( 14 )

    tEpsilon                     = aRParameters( 1 )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! INITIALIZE MUMPS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    tMUMPS => gSolvers( tSolverID )

    ! define the communicator
    tMUMPS%COMM = MPI_COMM_WORLD

    ! type of parallelism (PAR=1 host working, PAR=0 host not working)
    tMUMPS%PAR = tHostIsWorking

    ! symmetry setting (SYM=0 Unsymmetric, SYM=1 Sym. Positive Definite, SYM=2 General Symmetric)
    tMUMPS%SYM = tSymmetryMode


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Set PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Message streams: all THREE are suppressed at normal verbosity -- the
    ! C++ wrapper decodes INFO/INFOG into its own messages ( hard fail: the
    ! error box; soft fail: a box-styled line ), so MUMPS' own lines only
    ! break the log layout. -v 4 and up restores them.
    !
    ! Silencing the streams is the only lever that works here: ICNTL(4) = 0
    ! does NOT quiet everything the manual implies. In dmumps_driver.F the
    ! ERRORG section that emits
    !
    !     On return from DMUMPS, INFOG(1)=              -9
    !     On return from DMUMPS, INFOG(2)=         2327741
    !
    ! is guarded by MPG .gt. 0 alone -- i.e. by ICNTL(3), with no ICNTL(4)
    ! test -- and five further warning sites in that file share the pattern.
    ! Setting ICNTL(1) alone therefore left every soft -9 printing two raw
    ! lines through the middle of the timestep box.
    if( tInfoLevel .ge. 4 ) then
        ! set stream for errors
        tMUMPS%ICNTL( 1 ) = stderr

        ! set stream for diagnostics
        tMUMPS%ICNTL( 2 ) = stdout

        ! set stream for global information
        tMUMPS%ICNTL( 3 ) = stdout
    else
        tMUMPS%ICNTL( 1 ) = -1
        tMUMPS%ICNTL( 2 ) = -1
        tMUMPS%ICNTL( 3 ) = -1
    end if

    ! statistics
    select case( tInfoLevel )
    case( 5 )
        ! set info level
        tMUMPS%ICNTL( 4 ) = 2
    case( 4 )
        ! set info level
        tMUMPS%ICNTL( 4 ) = 2
    case default
        ! set info level
        tMUMPS%ICNTL( 4 ) = 0
    end select

    ! Matrix is already assembled
    tMUMPS%ICNTL( 5 ) = 0

    ! automatic permutation
    tMUMPS%ICNTL( 6 ) = 7

    tMUMPS%ICNTL( 7 ) = tSerialPermutationOrdering

    ! scaling strategy
    tMUMPS%ICNTL( 8 ) = 77 ! automatic choice

    ! refinement steps
    tMUMPS%ICNTL( 10 ) = tNumRefinementSteps

    ! statistics mode
    tMUMPS%ICNTL( 11 ) = tErrorAnalysisMode

    if( tMemoryRelaxation .gt. 0 ) then
        tMUMPS%ICNTL( 14 ) = tMemoryRelaxation
    end if

    ! ICNTL(23): per-process working-memory cap in MB. Written on EVERY
    ! rank with the same value ( the C++ side reduces it to a rank-uniform
    ! number first ), so MUMPS interprets it locally and identically. 0
    ! leaves MUMPS sizing from its analysis estimate
    if( tMemoryBudget .gt. 0 ) then
        tMUMPS%ICNTL( 23 ) = tMemoryBudget
    end if

    ! Matrix is centralized on the host
    tMUMPS%ICNTL( 18 ) = 0

    ! Right and side is always dense
    tMUMPS%ICNTL( 20 ) = 0

    ! Parallel ordering tool ( PT-SCOTCH / ParMETIS )
    tMUMPS%ICNTL( 29 ) = tParallelPermutationOrdering

    ! Flag to be set if determinant shall be computed
    tMUMPS%ICNTL( 33 ) = tComputeDeterminant

    ! for block low-ranking
    tMUMPS%ICNTL( 35 ) = tCompressionMode
    tMUMPS%CNTL( 7 )   = tEpsilon

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! LINK TO MATRIX DATA
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( tMUMPS%MYID .eq. tMasterRank ) then
        tMUMPS%n    = aN
        tMUMPS%nz   = aNNZ
        tMUMPS%nrhs = aNRHS
        tMUMPS%lrhs = aN
        tMUMPS%irn  => aRowIndices
        tMUMPS%jcn  => aColIndices
        tMUMPS%A    => aValues
        tMUMPS%rhs  => aX
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! SOLVE THE PROBLEM
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! compute memory capacity
    tCapacity = aN * aNRHS

    ! copy solution into X-Vector, MUMPS wants that
    forall ( k=1:tCapacity ) aX( k ) = aY( k )

    ! set job to calculate
    tMUMPS%JOB = tJob

    ! wait for other procs before we begin
    call MPI_BARRIER( MPI_COMM_WORLD, tStatus )

    ! call mumps
    call DMUMPS( tMUMPS )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! FINALIZE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! copy local info into output. INFO is rank-LOCAL ( a failing rank holds
    ! the true code, the others the propagated -1 ); INFOG carries the SAME
    ! code and supplementary value on every rank ( MUMPS user guide: "all
    ! processors would return with INFOG(1) = -8 and INFOG(2)=1000" ), which
    ! is what lets the C++ side key a collective retry on it
    ! all 80 entries: the BLR memory estimate INFOG(36) sits past the old
    ! 40-entry copy
    forall( k=1:80 ) aInfo( k )   = tMUMPS%INFO( k )
    forall( k=1:80 ) aInfoG( k )  = tMUMPS%INFOG( k )
    forall( k=1:20 ) aRInfoG( k ) = tMUMPS%RINFOG( k )

end subroutine mumpstools_solve

!------------------------------------------------------------------------------
