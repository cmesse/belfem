!> List of Parameters
!> 1: Master Rank         : 0 ( default )
!> 2: Parallel mode       : 0 - host not working
!>                        : 1 - host working
!> 3: Symmetry mode       : 0 - Unsymmetric
!>                        : 1 - Positive Definite
!>                        : 2 - General Symmetric
!> 4: Info Level          : 0 - silent
!>                        : 4 - some information
!>                        : 5 - all information
!> 5: Serial Reordering   : 0 - AMD
!>                        : 1 - user pivot
!>                        : 2 - AMF
!>                        : 3 - SCOTCH
!>                        : 4 - PORD
!>                        : 5 - METIS
!>                        : 6 - QAMD
!>                        : 7 - AUTOMATIC
!> 6: Parallel Reordering : 0 - AUTOMATIC
!>                        : 1 - PT-SCOTCH
!>                        : 2 - PARMETIS
!> 7: Refinement          :<0 - fixed number of steps
!>                        : 0 - none
!>                        :>0 - maxumum number of refinement steps
!> 8: Compute Determinant : 0 - no
!>                        : 1 - yes
!> 11: Error Analysis    : 0 - none
!>                       : 1 - full
!>                       : 2 - main

subroutine mumpstools_solve( &
        aIParameters, & !> list of integer parameters, see above
        aRParameters, & !> list of real parameters
        aN,          & !> size of matrix
        aNNZ,        & !> number of nonzeros
        aNRHS,       & !> number of cols on right hand side
        aRowIndices, & !> row indices of matrix
        aColIndices, & !> column indices of matrix
        aValues,     & !> data in matrix
        aX,          & !> left hand side
        aY,          & !> right hand side
        aInfo,       & !> information for debugging
        aRInfoG      & !> information for debugging
        ) bind( c )
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only : &
            stdout=>output_unit, &
            stderr=>error_unit
    implicit none
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! ARGUMENTS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    integer( c_int ), intent( in    ), dimension( 11 )                 :: aIParameters
    real( c_double ), intent( in ),    dimension( 1 )                  :: aRParameters
    integer( c_int ), intent( in    )                                  :: aN
    integer( c_int ), intent( in    )                                  :: aNNZ
    integer( c_int ), intent( in    )                                  :: aNRHS
    integer( c_int ), intent( in    ), dimension( aNNZ ), target       :: aRowIndices
    integer( c_int ), intent( in    ), dimension( aNNZ ), target       :: aColIndices
    real( c_double ), intent( in    ), dimension( aNNZ ), target       :: aValues
    real( c_double ), intent( inout ), dimension( aN * aNRHS ), target :: aX
    real( c_double ), intent( in    ), dimension( aN * aNRHS ), target :: aY
    integer( c_int ), intent( out   ), dimension( 40  )                :: aInfo
    real( c_double ), intent(out ), dimension( 20  )                   :: aRInfoG

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Shortcuts for better code readability of user parameters
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    integer :: tMasterRank
    integer :: tParallelMode
    integer :: tSymmetryMode
    integer :: tInfoLevel
    integer :: tSerialPermutationOrdering
    integer :: tParallelPermutationOrdering
    integer :: tNumRefinementSteps
    integer :: tComputeDeterminant
    integer :: tMemoryRelaxation
    integer :: tBlrMode
    integer :: tErrorAnalysisMode
    double precision :: tEpsilon

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! MUMPS stuff
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    include 'mpif.h'
    include 'dmumps_struc.h'

    !> the parameters object
    type ( DMUMPS_STRUC ):: tMUMPS

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Other parameters
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! iteration index
    integer :: k

    ! length of rhs memory
    integer :: tCapacity

    ! MPI status
    integer :: tStatus

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! POPULATE DEFINED PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    tMasterRank                  = aIParameters(  1 )
    tParallelMode                = aIParameters(  2 )
    tSymmetryMode                = aIParameters(  3 )
    tInfoLevel                   = aIParameters(  4 )
    tSerialPermutationOrdering   = aIParameters(  5 )
    tParallelPermutationOrdering = aIParameters(  6 )
    tNumRefinementSteps          = aIParameters(  7 )
    tComputeDeterminant          = aIParameters(  8 )
    tMemoryRelaxation            = aIParameters(  9 )
    tBlrMode                     = aIParameters( 10 )
    tErrorAnalysisMode           = aIParameters( 11 )

    tEpsilon                     = aRParameters( 1 )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! INITIALIZE MUMPS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! define the communicator
    tMUMPS%COMM = MPI_COMM_WORLD

    ! type of parallelism (PAR=1 host working, PAR=0 host not working)
    tMUMPS%PAR = tParallelMode

    ! symmetry setting (SYM=0 Unsymmetric, SYM=1 Sym. Positive Definite, SYM=2 General Symmetric)
    tMUMPS%SYM = tSymmetryMode

    ! procedure to call
    tMUMPS%JOB = -1

    ! initialize MUMPS
    call MPI_BARRIER( MPI_COMM_WORLD, tStatus )
    call DMUMPS( tMUMPS )


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Set PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! set stream for errors
    tMUMPS%ICNTL( 1 ) = stderr

    ! set stream for diagnostics
    tMUMPS%ICNTL( 2 ) = stdout

    ! set stream for global information
    tMUMPS%ICNTL( 3 ) = stdout

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

    ! no zero-free permutation needed
    tMUMPS%ICNTL( 6 ) = 0

    ! permutation ordering
    tMUMPS%ICNTL( 7 ) = tSerialPermutationOrdering

    ! scaling strategy
    tMUMPS%ICNTL( 8 ) = 0

    ! refinement steps
    tMUMPS%ICNTL( 10 ) = tNumRefinementSteps

    ! statistics mode
    tMUMPS%ICNTL( 11 ) = tErrorAnalysisMode

    if( tMemoryRelaxation .gt. 0 ) then
        tMUMPS%ICNTL( 14 ) = tMemoryRelaxation
    end if

    ! Matrix is centralized on the host
    tMUMPS%ICNTL( 18 ) = 0

    ! Right and side is always dense
    tMUMPS%ICNTL( 20 ) = 0

    ! Parallel refinement
    tMUMPS%ICNTL( 29 ) = tParallelPermutationOrdering

    ! Flag to be set if determinant shall be computed
    tMUMPS%ICNTL( 33 ) = tComputeDeterminant

    ! for block low-ranking
    tMUMPS%ICNTL( 35 ) = tBlrMode
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
    tMUMPS%JOB = 6

    ! wait for other procs before we begin
    call MPI_BARRIER( MPI_COMM_WORLD, tStatus )

    ! call mumps
    call DMUMPS( tMUMPS )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! FINALIZE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! copy local info into output
    forall( k=1:40 ) aInfo( k )   = tMUMPS%INFO( k )
    forall( k=1:20 ) aRInfoG( k ) = tMUMPS%RINFOG( k )

    ! tidy up
    tMUMPS%JOB = -2
    call DMUMPS( tMUMPS )

end subroutine mumpstools_solve

!------------------------------------------------------------------------------