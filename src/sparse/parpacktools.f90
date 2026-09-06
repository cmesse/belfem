
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

! Distributed counterpart of arpacktools.f90. The serial driver there owns the
! whole matrix on one rank ; here every rank owns a contiguous block of ROWS
! and PARPACK distributes the Arnoldi basis with it, so the n x ncv basis --
! the largest object in the algorithm -- is split as well.
!
! Row distribution contract ( the caller owns it ):
!   - pointers  are LOCAL and one-based : pointers( 1 ) == 1,
!               pointers( nloc + 1 ) == nnz + 1
!   - indices   are GLOBAL and one-based. They are NOT rebased : a local row
!               references columns anywhere in the global vector, which is why
!               x is assembled below
!   - the row blocks are contiguous and ordered by rank, so rank p owns the
!     rows that follow rank p-1. MPI_ALLGATHERV relies on that ordering
!
! CALLER CONTRACT -- cross rank uniformity ( NOT enforced here, by decision ):
!   PARPACK executes one algorithm whose control flow runs redundantly on every
!   rank ; the ranks stay in lockstep only because they take the same branches
!   and therefore reach the same reductions in the same order. Every rank MUST
!   pass identical nglobal, nev, job, tol and maxit. Only nloc, nnz, values,
!   indices and pointers are legitimately per rank. A divergent nev or job
!   changes the loop trip count inside pdnaupd and hangs the job rather than
!   failing. The one exception is nglobal, which falls out of the row map check
!   below for free.

subroutine parpack_standard_eigen( nloc, nglobal, nnz, values, indices, pointers, &
        job, nev, ncvmin, tol, maxit, sigma, lambdareal, lambdaimag, info ) bind( c )
#ifdef BELFEM_OMP
    use omp_lib
#endif
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

!-----------------------------------------------------------------------
!   input/output parameters
!-----------------------------------------------------------------------

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix data ( the LOCAL row block )
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: nloc     ! rows owned by this rank
    integer( int_t ), intent( in ) :: nglobal  ! rows of the whole matrix
    integer( int_t ), intent( in ) :: nnz      ! nonzeros in the local block
    real( c_double ), dimension( nnz ),      intent( in ) :: values    ! local values
    integer( int_t ), dimension( nnz ),      intent( in ) :: indices   ! GLOBAL column indices
    integer( int_t ), dimension( nloc + 1 ), intent( in ) :: pointers  ! LOCAL row pointers

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! settings
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: job      ! 0: smallest values, 1: largest values
    integer( int_t ), intent( in ) :: nev      ! number of values to compute
    integer( int_t ), intent( in ) :: ncvmin   ! floor for the Krylov subspace
    real( c_double ), intent( in ) :: tol      ! epsilon environment
    integer( int_t ), intent( in ) :: maxit    ! maximum number of restart iterations

    ! spectral fold, identical in meaning to the serial driver. sigma = 0 runs
    ! OP = A ; any other value runs OP = sigma*I - A, whose eigenvalues are
    ! sigma - lambda. Folding about a sigma just above max|lambda| turns the
    ! SMALL end of A into the LARGE end of OP, which is the end Arnoldi
    ! converges. Every rank must be given the SAME sigma -- it changes the
    ! operator, so a divergent value would put the ranks on different problems
    real( c_double ), intent( in ) :: sigma

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! output values ( identical on every rank -- pdneupd is collective )
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdareal
    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdaimag

    ! status array, laid out exactly as in arpacktools.f90 so one decoder
    ! serves both drivers
    !   info( 1 ) : pdnaupd error flag, or a driver code :
    !                 100 unexpected reverse communication request
    !                 101 the row blocks do not sum to nglobal, or nglobal
    !                     is not the same on every rank
    !                 102 at least one rank owns no rows
    !                 103 nnz does not match pointers( nloc + 1 ) - 1
    !                 104 job is neither 0 nor 1
    !               For 101-104 nothing was computed and info( 3:7 ) stay zero;
    !               for 100 the loop ran, so info( 3:7 ) hold the counters
    !   info( 2 ) : pdneupd error flag
    !   info( 3 ) : NCONV,  number of converged Ritz values ( iparam( 5 ) )
    !   info( 4 ) : MXITER, restart iterations taken        ( iparam( 3 ) )
    !   info( 5 ) : NUMOP,  number of OP*x operations       ( iparam( 9 ) )
    !   info( 6 ) : NUMREO, number of re-orthogonalizations ( iparam( 11 ) )
    !   info( 7 ) : NCV,    size of the Krylov subspace actually used
    integer( int_t ), intent( inout ), dimension( 7 ) :: info

!-----------------------------------------------------------------------
!   variables required for the pdnaupd loop
!-----------------------------------------------------------------------
    character, parameter :: bmat = 'I'     ! standard problem
    character( len=2 ) :: which
    integer( int_t ), dimension( 11 ) :: iparam
    integer( int_t ), dimension( 14 ) :: ipntr   ! 14 for the NON-symmetric driver

    integer( int_t ) :: ido

    ! PARPACK OVERWRITES tol when it is <= 0 : pdnaupd.f:564 assigns
    ! tol = pdlamch10( comm, 'EpsMach' ). The dummy is intent( in ), and the C++
    ! side hands over a pointer to a const member, so the library is given this
    ! writable copy instead -- the prescribed tolerance cannot drift
    real( c_double ) :: tTol

    real( c_double ), dimension( : ),    allocatable         :: residuals
    real( c_double ), dimension( :, : ), allocatable         :: vectors

    real( c_double ), dimension( : ),    allocatable, target :: workd
    real( c_double ), dimension( : ),    allocatable         :: workl

    ! contiguous: both point into workd with unit stride. Saying so lets the
    ! matvec address them directly and keeps the MPI_ALLGATHERV send buffer
    ! free of a copy in / copy out through the implicit mpif.h interface
    real( c_double ), dimension( : ), pointer, contiguous :: x   ! local input slice
    real( c_double ), dimension( : ), pointer, contiguous :: y   ! local output slice

    ! the assembled input vector. A local row may reference any global column,
    ! so the full vector is needed on every rank
    real( c_double ), dimension( : ), allocatable :: xglobal

    integer( int_t ) :: ncv
    integer( int_t ) :: lworkl
    integer( int_t ) :: i, j
    real( c_double ) :: acc          ! row accumulator, see the matvec below

    ! +1 leaves the matrix term as it is ( OP = A ), -1 flips it for the fold
    ! ( OP = sigma*I - A ). Branched once, outside the loop
    real( c_double ) :: tSign

!-----------------------------------------------------------------------
!   variables required for the pdneupd
!-----------------------------------------------------------------------
    ! eigenVALUES only. rvec = .false. skips the Ritz vector computation, and
    ! pdneupd then does not reference Z at all ( "If RVEC = .FALSE. or
    ! HOWMNY = 'P', then Z is not referenced", and "In any case, LDZ >= 1" ),
    ! so z collapses to a placeholder instead of a second n x ncv block.
    ! To get eigenvectors back: rvec = .true., z( nloc, ncv ), ldz = nloc.
    ! The type is logical because pdneupd declares it so ( pdneupd.f:335 ) --
    ! this constant is internal and never crosses the C++ boundary, so the
    ! integer-across-the-boundary rule does not apply to it
    logical, parameter :: rvec = .false.
    integer( int_t ), parameter :: ldz = 1

    character, parameter :: howmny = 'A'
    real( c_double ), dimension( :, : ), allocatable :: z
    real( c_double ), dimension( : ),    allocatable :: workev

    ! dneupd documents dr/di as NEV+1, but its own info = 1 remedy is "increase
    ! DR and DI to at least NCV". Sizing them at NCV from the start costs a few
    ! kB and removes that failure mode ; the caller still only gets NEV+1
    real( c_double ), dimension( : ), allocatable :: dr
    real( c_double ), dimension( : ), allocatable :: di

    real( c_double ), parameter :: sigmar = 0.0d0
    real( c_double ), parameter :: sigmai = 0.0d0

    logical, dimension( : ), allocatable :: slct

!-----------------------------------------------------------------------
!   MPI bookkeeping
!-----------------------------------------------------------------------
    ! MPI counts are default INTEGER, which is NOT int_t when BELFEM is built
    ! with 64 bit indices -- keep the two apart and convert explicitly
    integer :: tComm, tNumProcs, tError
    integer :: tMyRows, tSumRows
    integer :: tLocalError, tGlobalError
    integer, dimension( : ), allocatable :: tCounts   ! rows per rank
    integer, dimension( : ), allocatable :: tDispls   ! offset of each rank
    integer :: p

    tComm = MPI_COMM_WORLD
    call MPI_COMM_SIZE( tComm, tNumProcs, tError )

    info = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! row map. Gathered once, then reused for every matvec
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    allocate( tCounts( tNumProcs ) )
    allocate( tDispls( tNumProcs ) )

    tMyRows = int( nloc )
    call MPI_ALLGATHER( tMyRows, 1, MPI_INTEGER, &
                        tCounts, 1, MPI_INTEGER, tComm, tError )

    tDispls( 1 ) = 0
    tSumRows     = tCounts( 1 )
    do p = 2, tNumProcs
        tDispls( p ) = tDispls( p - 1 ) + tCounts( p - 1 )
        tSumRows     = tSumRows + tCounts( p )
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preconditions. Each is tested LOCALLY, then reduced with MPI_MAX so the
    ! decision to bail is identical on every rank. A local branch here would
    ! be the very hazard it is meant to catch: one rank returning while the
    ! others walk into the collective pdnaupd hangs the job
    !
    ! nglobal is checked against the reduced row map, which also catches a
    ! caller that passes a different nglobal on different ranks -- that rank
    ! raises 101 and MPI_MAX hands it to everyone
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    tLocalError = 0

    if ( nloc .le. 0 ) then
        ! pdnaupd rejects n <= 0 ( pdnaupd.f:522 ) and returns with ido = 99,
        ! so an empty rank would leave the reverse communication loop while
        ! the others are still inside collective calls
        tLocalError = 102
    else if ( pointers( nloc + 1 ) - 1 .ne. nnz ) then
        tLocalError = 103
    else if ( job .ne. 0 .and. job .ne. 1 ) then
        tLocalError = 104
    else if ( tSumRows .ne. int( nglobal ) ) then
        tLocalError = 101
    end if

    if ( tError .ne. MPI_SUCCESS ) tLocalError = 105

    call MPI_ALLREDUCE( tLocalError, tGlobalError, 1, MPI_INTEGER, MPI_MAX, &
                        tComm, tError )

    if ( tGlobalError .ne. 0 ) then
        info( 1 ) = tGlobalError
        deallocate( tDispls )
        deallocate( tCounts )
        return
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! size of the Krylov subspace. PARPACK only rejects ncv <= nev + 1, it
    ! does NOT bound ncv by the n it is handed -- and the n handed to pdnaupd
    ! is the LOCAL row count. ( The serial dnaupd DOES bound it, dnaupd.f:504
    ! versus pdnaupd.f:526, which is why the two drivers clamp differently. )
    ! The Krylov space lives in the GLOBAL space, so the clamp is against
    ! nglobal ; clamping against nloc would shrink the subspace as ranks are
    ! added, and nothing in the library would complain.
    ! ncvmin is the caller's floor -- raising it at fixed nev usually reduces
    ! the total OP*x count ( remark 4 of dnaupd ), at the price of a larger
    ! basis: vectors( nloc, ncv ) grows with it on every rank
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    ncv = 2*nev + 1
    if ( ncv < ncvmin ) ncv = ncvmin
    if ( ncv > nglobal ) ncv = nglobal

    lworkl = ncv * ( 3 * ncv + 8 )

    ! every array is LOCAL apart from xglobal
    allocate( residuals( nloc ) )
    allocate( vectors( nloc, ncv ) )
    allocate( workd( 3 * nloc ) )
    allocate( workl( lworkl ) )
    allocate( z( ldz, 1 ) )         ! placeholder, not referenced for rvec = .false.
    allocate( workev( 3 * ncv ) )
    allocate( slct( ncv ) )
    allocate( dr( ncv ) )
    allocate( di( ncv ) )
    allocate( xglobal( nglobal ) )

    ! set the which flag
    if( job .eq. 0 ) then
        which = 'SM'
    else
        which = 'LM'
    end if

    ! exact comparison on purpose: the C++ side passes a literal 0.0 for the
    ! unfolded operator, never a computed value that might round to it
    if( sigma .eq. 0.0D0 ) then
        tSign = 1.0D0
    else
        tSign = -1.0D0
    end if

    ! set up parameters
    iparam = 0
    iparam( 1 ) = 1     ! exact shifts
    iparam( 3 ) = maxit ! maximum number of restart iterations
    iparam( 7 ) = 1     ! mode 1, standard eigenvalue problem

    residuals = 0.0D0
    tTol      = tol

    ! begin reverse communication loop
    ido = 0

    do

        call pdnaupd( &
                tComm, &
                ido, &
                bmat, &
                nloc, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                nloc, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 1 ) )

        if (ido .eq. -1 .or. ido .eq. 1) then

            x => workd(ipntr(1):ipntr(1)+nloc-1)
            y => workd(ipntr(2):ipntr(2)+nloc-1)

            ! assemble the global input. This is the price of a row
            ! distribution: a local row can reference any column, so every
            ! rank needs the whole vector. It is also the scaling limit of
            ! this driver -- the message is O( nglobal ) per matvec.
            ! tError is deliberately not branched on inside the loop: MPI's
            ! default handler is MPI_ERRORS_ARE_FATAL, and a rank reacting to
            ! a failure on its own would deadlock the ranks that did not
            call MPI_ALLGATHERV( x, tMyRows, MPI_DOUBLE_PRECISION, &
                                 xglobal, tCounts, tDispls, MPI_DOUBLE_PRECISION, &
                                 tComm, tError )

            ! local rows against the assembled vector. The row sum goes into a
            ! scalar rather than into y( i ): x and y are both pointers into
            ! workd, so the compiler cannot prove they do not overlap and
            ! would reload y( i ) on every inner iteration
#ifdef BELFEM_OMP
!$omp parallel do private(i, j, acc) shared(nloc, pointers, values, indices, xglobal, x, y, sigma, tSign) schedule(guided)
#endif
            do i = 1, nloc
                acc = 0.0D0
                do j = pointers(i), pointers(i+1) - 1
                    acc = acc + values(j) * xglobal(indices(j))
                end do
                ! the identity term takes the LOCAL entry x( i ), not
                ! xglobal( i ) : row i of this rank is global row
                ! tDispls( rank + 1 ) + i, so xglobal( i ) would be the wrong
                ! component on every rank above zero. sigma = 0 collapses this
                ! to y( i ) = acc
                y(i) = sigma * x(i) + tSign * acc
            end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

        else if (ido .eq. 2) then
            ! B is the identity for bmat = 'I', so this reduces to a copy
            x => workd(ipntr(1):ipntr(1)+nloc-1)
            y => workd(ipntr(2):ipntr(2)+nloc-1)
            y = x
        else
            ! ido = 99 is the normal exit, anything else is unexpected
            if ( ido .ne. 99 ) info( 1 ) = 100
            exit
        end if
    end do

    info( 3 ) = iparam( 5 )  ! NCONV
    info( 4 ) = iparam( 3 )  ! MXITER
    info( 5 ) = iparam( 9 )  ! NUMOP
    info( 6 ) = iparam( 11 ) ! NUMREO
    info( 7 ) = ncv

    ! pdnaupd info 1 ( maximum iterations ) and 3 ( no shifts applied ) are
    ! warnings -- what converged is still extractable. This branch is uniform
    ! across ranks because pdnaupd returns the same flag everywhere
    if ( info( 1 ) .ge. 0 .and. info( 1 ) .ne. 100 ) then

        call pdneupd( &
                tComm, &
                rvec, &
                howmny, &
                slct, &
                dr, &
                di, &
                z, &
                ldz, &
                sigmar, &
                sigmai, &
                workev, &
                bmat, &
                nloc, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                nloc, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 2 ) )

        info( 3 ) = iparam( 5 )

        ! hand back what the caller has room for
        do i = 1, nev + 1
            lambdareal( i ) = dr( i )
            lambdaimag( i ) = di( i )
        end do

    end if

    deallocate( xglobal )
    deallocate( di )
    deallocate( dr )
    deallocate( slct )
    deallocate( workev )
    deallocate( z )
    deallocate( workl )
    deallocate( workd )
    deallocate( vectors )
    deallocate( residuals )
    deallocate( tDispls )
    deallocate( tCounts )

end subroutine parpack_standard_eigen


!-----------------------------------------------------------------------
!
! Symmetric counterpart of parpack_standard_eigen. Identical contract,
! identical info layout, identical row-distribution rules -- the difference is
! pdsaupd / pdseupd instead of pdnaupd / pdneupd. See arpack_symmetric_eigen
! in arpacktools.f90 for why the symmetric pair is the right one when the
! matrix is symmetric, and note that the ERROR TABLES differ: these flags are
! decoded by check_saupd / check_seupd, never by the nonsymmetric pair.
!
! The CALLER decides which driver to use. It is not detected here.
!
!-----------------------------------------------------------------------
subroutine parpack_symmetric_eigen( nloc, nglobal, nnz, values, indices, pointers, &
        job, nev, ncvmin, tol, maxit, sigma, lambdareal, lambdaimag, info ) bind( c )
#ifdef BELFEM_OMP
    use omp_lib
#endif
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

!-----------------------------------------------------------------------
!   input/output parameters
!-----------------------------------------------------------------------

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix data ( the LOCAL row block )
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: nloc     ! rows owned by this rank
    integer( int_t ), intent( in ) :: nglobal  ! rows of the whole matrix
    integer( int_t ), intent( in ) :: nnz      ! nonzeros in the local block
    real( c_double ), dimension( nnz ),      intent( in ) :: values    ! local values
    integer( int_t ), dimension( nnz ),      intent( in ) :: indices   ! GLOBAL column indices
    integer( int_t ), dimension( nloc + 1 ), intent( in ) :: pointers  ! LOCAL row pointers

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! settings
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: job      ! 0: smallest values, 1: largest values
    integer( int_t ), intent( in ) :: nev      ! number of values to compute
    integer( int_t ), intent( in ) :: ncvmin   ! floor for the Krylov subspace
    real( c_double ), intent( in ) :: tol      ! epsilon environment
    integer( int_t ), intent( in ) :: maxit    ! maximum number of restart iterations

    ! spectral fold, identical in meaning to the serial driver. sigma = 0 runs
    ! OP = A ; any other value runs OP = sigma*I - A, whose eigenvalues are
    ! sigma - lambda. Folding about a sigma just above max|lambda| turns the
    ! SMALL end of A into the LARGE end of OP, which is the end Arnoldi
    ! converges. Every rank must be given the SAME sigma -- it changes the
    ! operator, so a divergent value would put the ranks on different problems
    real( c_double ), intent( in ) :: sigma

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! output values ( identical on every rank -- pdseupd is collective )
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdareal
    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdaimag

    ! status array, laid out exactly as in arpacktools.f90 so one decoder
    ! serves both drivers
    !   info( 1 ) : pdsaupd error flag, or a driver code :
    !                 100 unexpected reverse communication request
    !                 101 the row blocks do not sum to nglobal, or nglobal
    !                     is not the same on every rank
    !                 102 at least one rank owns no rows
    !                 103 nnz does not match pointers( nloc + 1 ) - 1
    !                 104 job is neither 0 nor 1
    !               For 101-104 nothing was computed and info( 3:7 ) stay zero;
    !               for 100 the loop ran, so info( 3:7 ) hold the counters
    !   info( 2 ) : pdseupd error flag
    !   info( 3 ) : NCONV,  number of converged Ritz values ( iparam( 5 ) )
    !   info( 4 ) : MXITER, restart iterations taken        ( iparam( 3 ) )
    !   info( 5 ) : NUMOP,  number of OP*x operations       ( iparam( 9 ) )
    !   info( 6 ) : NUMREO, number of re-orthogonalizations ( iparam( 11 ) )
    !   info( 7 ) : NCV,    size of the Krylov subspace actually used
    integer( int_t ), intent( inout ), dimension( 7 ) :: info

!-----------------------------------------------------------------------
!   variables required for the pdsaupd loop
!-----------------------------------------------------------------------
    character, parameter :: bmat = 'I'     ! standard problem
    character( len=2 ) :: which
    integer( int_t ), dimension( 11 ) :: iparam
    ! ELEVEN, not the nonsymmetric driver's fourteen: the symmetric driver
    ! has no complex workspace to point at
    integer( int_t ), dimension( 11 ) :: ipntr

    integer( int_t ) :: ido

    ! PARPACK OVERWRITES tol when it is <= 0 : pdsaupd assigns
    ! tol = pdlamch10( comm, 'EpsMach' ). The dummy is intent( in ), and the C++
    ! side hands over a pointer to a const member, so the library is given this
    ! writable copy instead -- the prescribed tolerance cannot drift
    real( c_double ) :: tTol

    real( c_double ), dimension( : ),    allocatable         :: residuals
    real( c_double ), dimension( :, : ), allocatable         :: vectors

    real( c_double ), dimension( : ),    allocatable, target :: workd
    real( c_double ), dimension( : ),    allocatable         :: workl

    ! contiguous: both point into workd with unit stride. Saying so lets the
    ! matvec address them directly and keeps the MPI_ALLGATHERV send buffer
    ! free of a copy in / copy out through the implicit mpif.h interface
    real( c_double ), dimension( : ), pointer, contiguous :: x   ! local input slice
    real( c_double ), dimension( : ), pointer, contiguous :: y   ! local output slice

    ! the assembled input vector. A local row may reference any global column,
    ! so the full vector is needed on every rank
    real( c_double ), dimension( : ), allocatable :: xglobal

    integer( int_t ) :: ncv
    integer( int_t ) :: lworkl
    integer( int_t ) :: i, j
    real( c_double ) :: acc          ! row accumulator, see the matvec below

    ! +1 leaves the matrix term as it is ( OP = A ), -1 flips it for the fold
    ! ( OP = sigma*I - A ). Branched once, outside the loop
    real( c_double ) :: tSign

!-----------------------------------------------------------------------
!   variables required for the pdseupd
!-----------------------------------------------------------------------
    ! eigenVALUES only. rvec = .false. skips the Ritz vector computation, and
    ! pdseupd then does not reference Z at all ( "If RVEC = .FALSE. or
    ! HOWMNY = 'P', then Z is not referenced", and "In any case, LDZ >= 1" ),
    ! so z collapses to a placeholder instead of a second n x ncv block.
    ! To get eigenvectors back: rvec = .true., z( nloc, ncv ), ldz = nloc.
    ! The type is logical because pdseupd declares it so --
    ! this constant is internal and never crosses the C++ boundary, so the
    ! integer-across-the-boundary rule does not apply to it
    logical, parameter :: rvec = .false.
    integer( int_t ), parameter :: ldz = 1

    character, parameter :: howmny = 'A'
    real( c_double ), dimension( :, : ), allocatable :: z

    ! ONE real array: the spectrum of a symmetric matrix is real, so there is
    ! no imaginary companion and no complex pair that can straddle the
    ! nev / nev+1 boundary. Sized at NCV rather than NEV for the same reason
    ! the nonsymmetric driver oversizes dr/di -- it costs a few kB and removes
    ! a failure mode
    real( c_double ), dimension( : ), allocatable :: d

    ! pdseupd's OWN sigma: the shift of shift-invert mode. We run mode 1, so it
    ! is zero and unreferenced. Deliberately NOT named sigma -- confusing it
    ! with the fold argument would silently change the operator
    real( c_double ), parameter :: sigmashift = 0.0d0

    logical, dimension( : ), allocatable :: slct

!-----------------------------------------------------------------------
!   MPI bookkeeping
!-----------------------------------------------------------------------
    ! MPI counts are default INTEGER, which is NOT int_t when BELFEM is built
    ! with 64 bit indices -- keep the two apart and convert explicitly
    integer :: tComm, tNumProcs, tError
    integer :: tMyRows, tSumRows
    integer :: tLocalError, tGlobalError
    integer, dimension( : ), allocatable :: tCounts   ! rows per rank
    integer, dimension( : ), allocatable :: tDispls   ! offset of each rank
    integer :: p

    tComm = MPI_COMM_WORLD
    call MPI_COMM_SIZE( tComm, tNumProcs, tError )

    info = 0

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! row map. Gathered once, then reused for every matvec
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    allocate( tCounts( tNumProcs ) )
    allocate( tDispls( tNumProcs ) )

    tMyRows = int( nloc )
    call MPI_ALLGATHER( tMyRows, 1, MPI_INTEGER, &
                        tCounts, 1, MPI_INTEGER, tComm, tError )

    tDispls( 1 ) = 0
    tSumRows     = tCounts( 1 )
    do p = 2, tNumProcs
        tDispls( p ) = tDispls( p - 1 ) + tCounts( p - 1 )
        tSumRows     = tSumRows + tCounts( p )
    end do

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! preconditions. Each is tested LOCALLY, then reduced with MPI_MAX so the
    ! decision to bail is identical on every rank. A local branch here would
    ! be the very hazard it is meant to catch: one rank returning while the
    ! others walk into the collective pdsaupd hangs the job
    !
    ! nglobal is checked against the reduced row map, which also catches a
    ! caller that passes a different nglobal on different ranks -- that rank
    ! raises 101 and MPI_MAX hands it to everyone
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    tLocalError = 0

    if ( nloc .le. 0 ) then
        ! pdsaupd rejects n <= 0 and returns with ido = 99,
        ! so an empty rank would leave the reverse communication loop while
        ! the others are still inside collective calls
        tLocalError = 102
    else if ( pointers( nloc + 1 ) - 1 .ne. nnz ) then
        tLocalError = 103
    else if ( job .ne. 0 .and. job .ne. 1 ) then
        tLocalError = 104
    else if ( tSumRows .ne. int( nglobal ) ) then
        tLocalError = 101
    end if

    if ( tError .ne. MPI_SUCCESS ) tLocalError = 105

    call MPI_ALLREDUCE( tLocalError, tGlobalError, 1, MPI_INTEGER, MPI_MAX, &
                        tComm, tError )

    if ( tGlobalError .ne. 0 ) then
        info( 1 ) = tGlobalError
        deallocate( tDispls )
        deallocate( tCounts )
        return
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! size of the Krylov subspace. PARPACK only rejects ncv <= nev + 1, it
    ! does NOT bound ncv by the n it is handed -- and the n handed to pdsaupd
    ! is the LOCAL row count. ( The serial dsaupd DOES bound it where the
    ! parallel one does not, which is why the two drivers clamp differently. )
    ! The Krylov space lives in the GLOBAL space, so the clamp is against
    ! nglobal ; clamping against nloc would shrink the subspace as ranks are
    ! added, and nothing in the library would complain.
    ! ncvmin is the caller's floor -- raising it at fixed nev usually reduces
    ! the total OP*x count ( remark 4 of dsaupd ), at the price of a larger
    ! basis: vectors( nloc, ncv ) grows with it on every rank
    ! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -
    ncv = 2*nev + 1
    if ( ncv < ncvmin ) ncv = ncvmin
    if ( ncv > nglobal ) ncv = nglobal

    ! THE SYMMETRIC SIZE. pdsaupd wants ncv*( ncv + 8 ) against pdnaupd's
    ! ncv*( 3*ncv + 8 ). Using the nonsymmetric size here only wastes memory ;
    ! using this one there would corrupt the heap
    lworkl = ncv * ( ncv + 8 )

    ! every array is LOCAL apart from xglobal
    allocate( residuals( nloc ) )
    allocate( vectors( nloc, ncv ) )
    allocate( workd( 3 * nloc ) )
    allocate( workl( lworkl ) )
    allocate( z( ldz, 1 ) )         ! placeholder, not referenced for rvec = .false.
    allocate( slct( ncv ) )
    allocate( d( ncv ) )
    allocate( xglobal( nglobal ) )

    ! set the which flag
    if( job .eq. 0 ) then
        which = 'SM'
    else
        which = 'LM'
    end if

    ! exact comparison on purpose: the C++ side passes a literal 0.0 for the
    ! unfolded operator, never a computed value that might round to it
    if( sigma .eq. 0.0D0 ) then
        tSign = 1.0D0
    else
        tSign = -1.0D0
    end if

    ! set up parameters
    iparam = 0
    iparam( 1 ) = 1     ! exact shifts
    iparam( 3 ) = maxit ! maximum number of restart iterations
    iparam( 7 ) = 1     ! mode 1, standard eigenvalue problem

    residuals = 0.0D0
    tTol      = tol

    ! begin reverse communication loop
    ido = 0

    do

        call pdsaupd( &
                tComm, &
                ido, &
                bmat, &
                nloc, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                nloc, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 1 ) )

        if (ido .eq. -1 .or. ido .eq. 1) then

            x => workd(ipntr(1):ipntr(1)+nloc-1)
            y => workd(ipntr(2):ipntr(2)+nloc-1)

            ! assemble the global input. This is the price of a row
            ! distribution: a local row can reference any column, so every
            ! rank needs the whole vector. It is also the scaling limit of
            ! this driver -- the message is O( nglobal ) per matvec.
            ! tError is deliberately not branched on inside the loop: MPI's
            ! default handler is MPI_ERRORS_ARE_FATAL, and a rank reacting to
            ! a failure on its own would deadlock the ranks that did not
            call MPI_ALLGATHERV( x, tMyRows, MPI_DOUBLE_PRECISION, &
                                 xglobal, tCounts, tDispls, MPI_DOUBLE_PRECISION, &
                                 tComm, tError )

            ! local rows against the assembled vector. The row sum goes into a
            ! scalar rather than into y( i ): x and y are both pointers into
            ! workd, so the compiler cannot prove they do not overlap and
            ! would reload y( i ) on every inner iteration
#ifdef BELFEM_OMP
!$omp parallel do private(i, j, acc) shared(nloc, pointers, values, indices, xglobal, x, y, sigma, tSign) schedule(guided)
#endif
            do i = 1, nloc
                acc = 0.0D0
                do j = pointers(i), pointers(i+1) - 1
                    acc = acc + values(j) * xglobal(indices(j))
                end do
                ! the identity term takes the LOCAL entry x( i ), not
                ! xglobal( i ) : row i of this rank is global row
                ! tDispls( rank + 1 ) + i, so xglobal( i ) would be the wrong
                ! component on every rank above zero. sigma = 0 collapses this
                ! to y( i ) = acc
                y(i) = sigma * x(i) + tSign * acc
            end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

        else if (ido .eq. 2) then
            ! B is the identity for bmat = 'I', so this reduces to a copy
            x => workd(ipntr(1):ipntr(1)+nloc-1)
            y => workd(ipntr(2):ipntr(2)+nloc-1)
            y = x
        else
            ! ido = 99 is the normal exit, anything else is unexpected
            if ( ido .ne. 99 ) info( 1 ) = 100
            exit
        end if
    end do

    info( 3 ) = iparam( 5 )  ! NCONV
    info( 4 ) = iparam( 3 )  ! MXITER
    info( 5 ) = iparam( 9 )  ! NUMOP
    info( 6 ) = iparam( 11 ) ! NUMREO
    info( 7 ) = ncv

    ! pdsaupd info 1 ( maximum iterations ) and 3 ( no shifts applied ) are
    ! warnings -- what converged is still extractable. This branch is uniform
    ! across ranks because pdsaupd returns the same flag everywhere
    if ( info( 1 ) .ge. 0 .and. info( 1 ) .ne. 100 ) then

        call pdseupd( &
                tComm, &
                rvec, &
                howmny, &
                slct, &
                d, &
                z, &
                ldz, &
                sigmashift, &
                bmat, &
                nloc, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                nloc, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 2 ) )

        info( 3 ) = iparam( 5 )

        ! hand back what the caller has room for. The spectrum is real, so the
        ! imaginary buffer is zeroed rather than left as the caller found it
        do i = 1, nev
            lambdareal( i ) = d( i )
        end do
        lambdareal( nev + 1 ) = 0.0D0

        do i = 1, nev + 1
            lambdaimag( i ) = 0.0D0
        end do

    end if

    deallocate( xglobal )
    deallocate( d )
    deallocate( slct )
    deallocate( z )
    deallocate( workl )
    deallocate( workd )
    deallocate( vectors )
    deallocate( residuals )
    deallocate( tDispls )
    deallocate( tCounts )

end subroutine parpack_symmetric_eigen
