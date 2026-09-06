
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

subroutine arpack_standard_eigen( n, nnz, values, indices, pointers, job, nev, ncvmin, tol, maxit, &
        sigma, lambdareal, lambdaimag, info ) bind( c )
#ifdef BELFEM_OMP
    use omp_lib
#endif
    use, intrinsic :: iso_c_binding
    implicit none
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

!-----------------------------------------------------------------------
!   input/output parameters
!-----------------------------------------------------------------------

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix data
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: n    ! size of matrix
    integer( int_t ), intent( in ) :: nnz  ! number of nonzeros
    real( c_double ),   dimension( nnz ),   intent( in ) :: values   ! values of matrix
    integer( int_t ), dimension( nnz ),   intent( in ) :: indices    ! CSR: column indices, CSC: row indices
    integer( int_t ), dimension( n + 1 ), intent( in ) :: pointers   ! pointers of matrix

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! settings
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: job                            ! 0: smallest values, 1: largest values
    integer( int_t ), intent( in ) :: nev                            ! number of values to compute
    integer( int_t ), intent( in ) :: ncvmin                         ! floor for the Krylov subspace
    real( c_double ), intent( in ) :: tol                            ! epsilon environment
    integer( int_t ), intent( in ) :: maxit                          ! maximum number of restart iterations

    ! spectral fold. sigma = 0 runs the plain operator OP = A ; any other value
    ! runs OP = sigma*I - A, whose eigenvalues are sigma - lambda. Folding about
    ! sigma = max|lambda| turns the SMALL end of A into the LARGE end of OP,
    ! which is where Arnoldi actually converges -- the caller subtracts the
    ! SAME sigma back out, so its own accuracy cancels exactly
    real( c_double ), intent( in ) :: sigma

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! output values
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdareal ! real parts of eigenvalues
    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdaimag ! imaginary parts of eigenvalues

    ! status array. dnaupd and dneupd have DIFFERENT error tables, so they
    ! get separate slots -- a shared one cannot be decoded on the C++ side
    !   info( 1 ) : dnaupd error flag ( 100 = unexpected reverse comm. request )
    !   info( 2 ) : dneupd error flag
    !   info( 3 ) : NCONV,  number of converged Ritz values ( iparam( 5 ) )
    !   info( 4 ) : MXITER, restart iterations taken        ( iparam( 3 ) )
    !   info( 5 ) : NUMOP,  number of OP*x operations       ( iparam( 9 ) )
    !   info( 6 ) : NUMREO, number of re-orthogonalizations ( iparam( 11 ) )
    !   info( 7 ) : NCV,    size of the Krylov subspace actually used
    integer( int_t ), intent( inout ), dimension( 7 ) :: info

!-----------------------------------------------------------------------
!   variables required for the dnaupd loop
!-----------------------------------------------------------------------
    character, parameter :: bmat = 'I'     ! the job flag, using 'I' for the standard problem
    character( len=2 ) :: which            ! tells what value we are looking for
    integer( int_t ), dimension( 11 ) :: iparam ! parameters for dneupd
    integer( int_t ), dimension( 14 ) :: ipntr   ! pointers for the work array

    integer( int_t ) :: ido                ! reverse communication flag

    ! ARPACK OVERWRITES tol when it is <= 0 : dnaupd assigns
    ! tol = dlamch( 'EpsMach' ). The dummy is intent( in ), and the C++ side
    ! hands over a pointer to a const member, so the library is given this
    ! writable copy instead -- the prescribed tolerance cannot drift
    real( c_double ) :: tTol

    ! the work arrays scale with n and ncv, so they live on the heap. as
    ! automatic arrays they were allocated on the stack ( gfortran keeps them
    ! there, and -fopenmp forces it ), which caps ncv at a few units before
    ! the stack overflows on a production sized matrix
    real( c_double ), dimension( : ),    allocatable         :: residuals ! residual vector
    real( c_double ), dimension( :, : ), allocatable         :: vectors   ! the arnoldi vectors

    real( c_double ), dimension( : ),    allocatable, target :: workd ! first work vector
    real( c_double ), dimension( : ),    allocatable         :: workl ! second work vector
    real( c_double ), dimension( : ), pointer :: x ! left hand side for multiplication
    real( c_double ), dimension( : ), pointer :: y ! right hand side for multiplication

    integer( int_t ) :: ncv
    integer( int_t ) :: lworkl
    integer( int_t ) :: i, j ! iterators

    ! row accumulator. It is NOT an optimization only: x and y are both
    ! pointers into workd, so the compiler cannot prove they do not alias and
    ! would reload y( i ) on every inner iteration ( the distributed driver
    ! has always done it this way ). The row bounds are inlined into the loop
    ! header for a second reason -- as named locals they were SHARED across
    ! the OpenMP team, and every thread wrote them
    real( c_double ) :: acc

    ! +1 leaves the matrix term as it is ( OP = A ), -1 flips it for the fold
    ! ( OP = sigma*I - A ). Branching once here keeps the inner loop clean
    real( c_double ) :: tSign

!-----------------------------------------------------------------------
!   variables required for the dneupd
!-----------------------------------------------------------------------
    ! eigenVALUES only. rvec = .false. skips the Ritz vector computation, and
    ! dneupd then does not reference Z at all ( "If RVEC = .FALSE. or
    ! HOWMNY = 'P', then Z is not referenced", and "In any case, LDZ >= 1" ),
    ! so z collapses to a placeholder instead of a full n x ( nev + 1 ) block.
    ! To get eigenvectors back: rvec = .true., z( n, nev + 1 ), ldz = n.
    ! logical, not integer, because dneupd declares it so -- this constant is
    ! internal and never crosses the C++ boundary
    logical, parameter :: rvec = .false.
    integer( int_t ), parameter :: ldz = 1

    character, parameter :: howmny = 'A'
    real( c_double ), dimension( :, : ), allocatable :: z
    real( c_double ), dimension( : ),    allocatable :: workev

    real( c_double ), parameter :: sigmar = 0.0d0
    real( c_double ), parameter :: sigmai = 0.0d0

    logical, dimension( : ), allocatable :: slct

    ! size of the Krylov subspace. ARPACK requires 2 <= ncv - nev and
    ! ncv <= n ; remark 4 of dnaupd recommends ncv >= 2*nev + 1. That lower
    ! bound is far too tight in practice -- at nev = 1 it gives a subspace of
    ! dimension 3, which restarts constantly and discards the Krylov
    ! information every cycle. ncvmin is the caller's floor, and remark 4 is
    ! explicit that raising it at fixed nev USUALLY REDUCES the total number
    ! of OP*x operations even though each restart costs more
    ncv = 2*nev + 1
    if ( ncv < ncvmin ) ncv = ncvmin
    if ( ncv > n ) ncv = n

    ! ARPACK reports info = -3 if this cannot be satisfied, which the caller
    ! decodes -- do not silently repair a matrix that is too small to work on
    lworkl = ncv * ( 3 * ncv + 8 )

    ! every array below is dimensioned per the dnaupd / dneupd contract.
    ! note that slct, workev and z scale with NCV or NEV+1, NOT with NEV:
    ! sizing them with nev overruns them inside dneupd
    allocate( residuals( n ) )
    allocate( vectors( n, ncv ) )
    allocate( workd( 3 * n ) )
    allocate( workl( lworkl ) )
    allocate( z( ldz, 1 ) )         ! placeholder, not referenced for rvec = .false.
    allocate( workev( 3 * ncv ) )   ! dneupd: 3*NCV
    allocate( slct( ncv ) )         ! dneupd: NCV, workspace when howmny = 'A'

    ! set the which flag
    if( job .eq. 0 ) then
        which = 'SM'
    else
        which = 'LM'
    end if

    ! sigma is compared against an exact zero on purpose: the C++ side passes a
    ! literal 0.0 for the unfolded operator, never a computed value that might
    ! round to it
    if( sigma .eq. 0.0D0 ) then
        tSign = 1.0D0
    else
        tSign = -1.0D0
    end if

    ! set up parameters
    iparam = 0
    iparam( 1 ) = 1   ! exact shifts
    iparam( 3 ) = maxit ! maximum number of restart iterations
    iparam( 7 ) = 1   ! set to 1 for standard eigenvalue problem

    ! reset residuals
    residuals = 0.0D0
    tTol      = tol

    ! rest info flags
    info = 0

    ! begin reverse communication loop
    ido = 0

    do

        call dnaupd( &
                ido, &
                bmat, &
                n, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                n, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 1 ) )

        if (ido .eq. -1 .or. ido .eq. 1) then

! perform multiplication ( note, while this is the csr routine, it does not matter for eigenvalues)
            x => workd(ipntr(1):ipntr(1)+n-1)
            y => workd(ipntr(2):ipntr(2)+n-1)

! guided, not dynamic: bare dynamic means chunk 1, so every row costs an
! atomic dispatch for a handful of flops. The chunks also shrink towards the
! end of the loop, which is where the dense cut/lambda rows live.
! acc is private and y( i ) is written once, so the y = 0.0 prefill the
! accumulating form needed is gone with it
#ifdef BELFEM_OMP
!$omp parallel do private(i, j, acc) shared(n, pointers, values, indices, x, y, sigma, tSign) schedule(guided)
#endif
            do i = 1, n
                acc = 0.0D0
                do j = pointers(i), pointers(i+1) - 1
                    acc = acc + values(j) * x(indices(j))
                end do
                ! sigma = 0 collapses this to y( i ) = acc
                y(i) = sigma * x(i) + tSign * acc
            end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

        else if (ido .eq. 2) then
            ! B is the identity for bmat = 'I', so this reduces to a copy.
            ! dnaupd does not ask for it in that mode, but servicing it is
            ! cheaper than risking a spin if it ever does
            x => workd(ipntr(1):ipntr(1)+n-1)
            y => workd(ipntr(2):ipntr(2)+n-1)
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

    ! dnaupd info 1 ( maximum iterations reached ) and 3 ( no shifts applied )
    ! are warnings -- the Ritz values found so far are still extractable, so
    ! dneupd is called. A negative flag means there is nothing to extract
    if ( info( 1 ) .ge. 0 .and. info( 1 ) .ne. 100 ) then

        call dneupd( &
                rvec, &
                howmny, &
                slct, &
                lambdareal, &
                lambdaimag, &
                z, &
                ldz, &
                sigmar, &
                sigmai, &
                workev, &
                bmat, &
                n, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                n, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 2 ) )

        ! dneupd recomputes the converged count
        info( 3 ) = iparam( 5 )

    end if

    deallocate( slct )
    deallocate( workev )
    deallocate( z )
    deallocate( workl )
    deallocate( workd )
    deallocate( vectors )
    deallocate( residuals )

end subroutine arpack_standard_eigen

!-----------------------------------------------------------------------
!
! Symmetric counterpart of arpack_standard_eigen. Same contract, same info
! layout, same fold argument -- the difference is dsaupd / dseupd instead of
! dnaupd / dneupd, which is the right pair when the matrix is symmetric:
!
!   - Lanczos uses a THREE-TERM recurrence where Arnoldi carries a full
!     Hessenberg column, so the re-orthogonalization work per restart drops
!     from O( ncv^2 ) to O( ncv ) vectors, and workl from ncv*( 3*ncv + 8 )
!     to ncv*( ncv + 8 )
!   - the spectrum is real, so there is no complex pair to straddle the
!     nev / nev+1 boundary and lambdaimag comes back as exact zeros
!   - for a symmetric matrix the spectral ratio lambda_max / lambda_min IS
!     the 2-norm condition number, which the nonsymmetric driver cannot
!     promise
!
! The CALLER decides which one to use. It is not detected here: the thermal
! Jacobian is symmetric by construction and the magnetic h-phi one is not
! ( Christian, 2026-08-28 ).
!
!-----------------------------------------------------------------------
subroutine arpack_symmetric_eigen( n, nnz, values, indices, pointers, job, nev, ncvmin, tol, maxit, &
        sigma, lambdareal, lambdaimag, info ) bind( c )
#ifdef BELFEM_OMP
    use omp_lib
#endif
    use, intrinsic :: iso_c_binding
    implicit none
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

!-----------------------------------------------------------------------
!   input/output parameters
!-----------------------------------------------------------------------

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix data
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: n    ! size of matrix
    integer( int_t ), intent( in ) :: nnz  ! number of nonzeros
    real( c_double ),   dimension( nnz ),   intent( in ) :: values   ! values of matrix
    integer( int_t ), dimension( nnz ),   intent( in ) :: indices    ! CSR: column indices, CSC: row indices
    integer( int_t ), dimension( n + 1 ), intent( in ) :: pointers   ! pointers of matrix

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! settings
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    integer( int_t ), intent( in ) :: job                            ! 0: smallest values, 1: largest values
    integer( int_t ), intent( in ) :: nev                            ! number of values to compute
    integer( int_t ), intent( in ) :: ncvmin                         ! floor for the Krylov subspace
    real( c_double ), intent( in ) :: tol                            ! epsilon environment
    integer( int_t ), intent( in ) :: maxit                          ! maximum number of restart iterations

    ! spectral fold, exactly as in the nonsymmetric driver. sigma = 0 runs
    ! OP = A ; any other value runs OP = sigma*I - A. Note this is OUR shift,
    ! applied in the matvec -- it is NOT dseupd's sigma, which belongs to
    ! shift-invert mode ( iparam( 7 ) = 3 ) and stays zero here
    real( c_double ), intent( in ) :: sigma

! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -
! output values
! - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -

    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdareal ! eigenvalues
    real( c_double ), dimension( nev + 1 ), intent( inout ) :: lambdaimag ! zeroed: the spectrum is real

    ! same seven slots as the nonsymmetric driver, so the caller's info
    ! handling does not fork. The ERROR TABLES differ, though -- dsaupd and
    ! dnaupd do not share codes -- so the C++ side decodes these with
    ! check_saupd / check_seupd, never with the nonsymmetric pair
    !   info( 1 ) : dsaupd error flag ( 100 = unexpected reverse comm. request )
    !   info( 2 ) : dseupd error flag
    !   info( 3 ) : NCONV,  number of converged Ritz values ( iparam( 5 ) )
    !   info( 4 ) : MXITER, restart iterations taken        ( iparam( 3 ) )
    !   info( 5 ) : NUMOP,  number of OP*x operations       ( iparam( 9 ) )
    !   info( 6 ) : NUMREO, number of re-orthogonalizations ( iparam( 11 ) )
    !   info( 7 ) : NCV,    size of the Krylov subspace actually used
    integer( int_t ), intent( inout ), dimension( 7 ) :: info

!-----------------------------------------------------------------------
!   variables required for the dsaupd loop
!-----------------------------------------------------------------------
    character, parameter :: bmat = 'I'     ! standard problem
    character( len=2 ) :: which
    integer( int_t ), dimension( 11 ) :: iparam

    ! ELEVEN, not fourteen: the symmetric driver has no complex workspace to
    ! point at. Sizing this like the nonsymmetric one would be harmless, but
    ! sizing the nonsymmetric one like this would overrun
    integer( int_t ), dimension( 11 ) :: ipntr

    integer( int_t ) :: ido

    ! see the nonsymmetric driver: ARPACK overwrites tol when it is <= 0, and
    ! the C++ side hands over a pointer to a const member
    real( c_double ) :: tTol

    real( c_double ), dimension( : ),    allocatable         :: residuals
    real( c_double ), dimension( :, : ), allocatable         :: vectors

    real( c_double ), dimension( : ),    allocatable, target :: workd
    real( c_double ), dimension( : ),    allocatable         :: workl
    real( c_double ), dimension( : ), pointer :: x
    real( c_double ), dimension( : ), pointer :: y

    integer( int_t ) :: ncv
    integer( int_t ) :: lworkl
    integer( int_t ) :: i, j

    ! row accumulator, and the fold sign. Same reasoning as the nonsymmetric
    ! driver: x and y both point into workd, so y( i ) cannot be held in a
    ! register, and the row bounds stay in the loop header rather than in
    ! shared locals
    real( c_double ) :: acc
    real( c_double ) :: tSign

!-----------------------------------------------------------------------
!   variables required for the dseupd
!-----------------------------------------------------------------------
    logical, parameter :: rvec = .false.
    integer( int_t ), parameter :: ldz = 1

    character, parameter :: howmny = 'A'
    real( c_double ), dimension( :, : ), allocatable :: z

    ! dseupd's OWN sigma: the shift of shift-invert mode. We run mode 1, so it
    ! is zero and unreferenced. Deliberately NOT named sigma -- confusing it
    ! with the fold argument above would silently change the operator
    real( c_double ), parameter :: sigmashift = 0.0d0

    logical, dimension( : ), allocatable :: slct

    ! eigenvalues come back here. dseupd fills D( NEV ) ; the caller's buffer
    ! is nev + 1 long, so this local is copied out rather than aliased, which
    ! also keeps the nev + 1 slot from carrying stale values
    real( c_double ), dimension( : ), allocatable :: d

    ! ncv rule as in the nonsymmetric driver: ARPACK requires 2 <= ncv - nev
    ! and ncv <= n, remark 4 recommends 2*nev + 1, and ncvmin is the caller's
    ! floor
    ncv = 2*nev + 1
    if ( ncv < ncvmin ) ncv = ncvmin
    if ( ncv > n ) ncv = n

    ! THE SYMMETRIC SIZE. dsaupd wants ncv*( ncv + 8 ) where dnaupd wants
    ! ncv*( 3*ncv + 8 ) -- the Lanczos recurrence is short, so the workspace
    ! is smaller. Using the nonsymmetric size here would merely waste memory ;
    ! using this one there would corrupt the heap
    lworkl = ncv * ( ncv + 8 )

    allocate( residuals( n ) )
    allocate( vectors( n, ncv ) )
    allocate( workd( 3 * n ) )
    allocate( workl( lworkl ) )
    allocate( z( ldz, 1 ) )         ! placeholder, not referenced for rvec = .false.
    allocate( slct( ncv ) )         ! dseupd: NCV
    allocate( d( ncv ) )            ! dseupd fills NEV of these

    ! set the which flag. 'SM' / 'LM' keeps the job contract identical to the
    ! nonsymmetric driver so the two are interchangeable from the caller's
    ! side. For a positive definite spectrum these coincide with 'SA' / 'LA'
    if( job .eq. 0 ) then
        which = 'SM'
    else
        which = 'LM'
    end if

    ! exact comparison on purpose, see the nonsymmetric driver
    if( sigma .eq. 0.0D0 ) then
        tSign = 1.0D0
    else
        tSign = -1.0D0
    end if

    iparam = 0
    iparam( 1 ) = 1   ! exact shifts
    iparam( 3 ) = maxit
    iparam( 7 ) = 1   ! mode 1, standard eigenvalue problem

    residuals = 0.0D0
    tTol      = tol
    d         = 0.0D0

    info = 0

    ido = 0

    do

        call dsaupd( &
                ido, &
                bmat, &
                n, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                n, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 1 ) )

        if (ido .eq. -1 .or. ido .eq. 1) then

            x => workd(ipntr(1):ipntr(1)+n-1)
            y => workd(ipntr(2):ipntr(2)+n-1)

! the matrix is symmetric, so CSR and CSC address the same operator and this
! one loop serves both storage orders
#ifdef BELFEM_OMP
!$omp parallel do private(i, j, acc) shared(n, pointers, values, indices, x, y, sigma, tSign) schedule(guided)
#endif
            do i = 1, n
                acc = 0.0D0
                do j = pointers(i), pointers(i+1) - 1
                    acc = acc + values(j) * x(indices(j))
                end do
                ! sigma = 0 collapses this to y( i ) = acc
                y(i) = sigma * x(i) + tSign * acc
            end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

        else if (ido .eq. 2) then
            ! B is the identity for bmat = 'I'
            x => workd(ipntr(1):ipntr(1)+n-1)
            y => workd(ipntr(2):ipntr(2)+n-1)
            y = x
        else
            if ( ido .ne. 99 ) info( 1 ) = 100
            exit
        end if
    end do

    info( 3 ) = iparam( 5 )  ! NCONV
    info( 4 ) = iparam( 3 )  ! MXITER
    info( 5 ) = iparam( 9 )  ! NUMOP
    info( 6 ) = iparam( 11 ) ! NUMREO
    info( 7 ) = ncv

    ! dsaupd info 1 ( maximum iterations ) is a warning: the Ritz values found
    ! so far are still extractable. A negative flag means there is nothing to
    ! extract. dsaupd also documents info = 3 ( no shifts could be applied,
    ! raise NCV ); this branch admits it as a warning, and check_saupd on
    ! the C++ side decodes it as one
    if ( info( 1 ) .ge. 0 .and. info( 1 ) .ne. 100 ) then

        call dseupd( &
                rvec, &
                howmny, &
                slct, &
                d, &
                z, &
                ldz, &
                sigmashift, &
                bmat, &
                n, &
                which, &
                nev, &
                tTol, &
                residuals, &
                ncv, &
                vectors, &
                n, &
                iparam, &
                ipntr, &
                workd, &
                workl, &
                lworkl, &
                info( 2 ) )

        ! dseupd recomputes the converged count
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

    deallocate( d )
    deallocate( slct )
    deallocate( z )
    deallocate( workl )
    deallocate( workd )
    deallocate( vectors )
    deallocate( residuals )

end subroutine arpack_symmetric_eigen

!-----------------------------------------------------------------------
!
! SHIFT-INVERT reverse-communication shims.
!
! Unlike arpack_standard_eigen / arpack_symmetric_eigen, these do NOT own
! the iteration -- they are single steps, and the loop lives in C++. That is
! the point: mode 3 applies OP = ( A - sigma*I )^-1, and only the C++ side can
! call the sparse solver that inverts. So C++ owns every array, calls
! arpack_si_step until ido = 99, services each request with a linear solve,
! then calls arpack_si_extract once.
!
! WHY SHIMS RATHER THAN dsaupd_ DIRECTLY FROM C++: this file exposes only
! NUMERIC arguments across bind( c ). BMAT, WHICH, HOWMNY, RVEC and SELECT
! stay Fortran-side, so no character descriptor and no Fortran LOGICAL ever
! crosses the boundary -- neither has a portable C representation, and the
! tree contains no name-mangled extern anywhere ( verified 2026-08-29 ).
!
! CONTRACT, verified against arpack-ng 3.9.1 ( SRC/dsaupd.f, SRC/dseupd.f ),
! which is the version BELFEM links ( libarpack.so.2.1.0 ) :
!
!   iparam( 7 ) = 3, bmat = 'I', which = 'LM', sigma = 0 asks for the
!   eigenvalues of A closest to zero -- i.e. the SMALLEST ones -- because
!   they are the LARGEST of A^-1, which is the end Arnoldi converges.
!
!   ido = -1 : y = OP * x, x at workd( ipntr( 1 ) )
!   ido =  1 : y = OP * x, x at workd( ipntr( 3 ) )   <-- NOT ipntr( 1 ).
!              In modes 3/4/5 B*x is already there and must not be recomputed
!   ido =  2 : y = B * x, which for bmat = 'I' is a copy. Legal but NOT
!              expected on this path -- handle it, never rely on seeing it
!   ido = 99 : done
!   result   : always to workd( ipntr( 2 ) )
!
!   dseupd applies lambda = 1/theta + sigma ITSELF, so the caller takes the
!   values as returned. Applying a reciprocal on top would invert twice.
!
! ONE PROBLEM AT A TIME: dsaupd carries SAVEd internal state across its
! reverse-communication calls, so two interleaved eigenproblems would corrupt
! each other. That matches BELFEM's deliberate not-thread-safe posture, but it
! means these shims must never be driven from two loops at once.
!
!-----------------------------------------------------------------------
subroutine arpack_si_step( ido, n, nev, ncv, tol, resid, vectors, ldv, &
        iparam, ipntr, workd, workl, lworkl, info ) bind( c )
    use, intrinsic :: iso_c_binding
    implicit none
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

    integer( int_t ), intent( inout ) :: ido     ! reverse communication flag
    integer( int_t ), intent( in )    :: n       ! size of the matrix
    integer( int_t ), intent( in )    :: nev     ! number of values wanted
    integer( int_t ), intent( in )    :: ncv     ! Krylov subspace size
    real( c_double ), intent( inout ) :: tol     ! dsaupd may overwrite this
    integer( int_t ), intent( in )    :: ldv
    integer( int_t ), intent( in )    :: lworkl

    real( c_double ), dimension( n ),        intent( inout ) :: resid
    real( c_double ), dimension( ldv, ncv ), intent( inout ) :: vectors
    integer( int_t ), dimension( 11 ),       intent( inout ) :: iparam
    integer( int_t ), dimension( 11 ),       intent( inout ) :: ipntr
    real( c_double ), dimension( 3 * n ),    intent( inout ) :: workd
    real( c_double ), dimension( lworkl ),   intent( inout ) :: workl
    integer( int_t ), intent( inout ) :: info

    ! kept Fortran-side, see the header note
    character, parameter :: bmat = 'I'
    character( len=2 ), parameter :: which = 'LM'

    call dsaupd( ido, bmat, n, which, nev, tol, resid, ncv, vectors, ldv, &
                 iparam, ipntr, workd, workl, lworkl, info )

end subroutine arpack_si_step

!-----------------------------------------------------------------------
subroutine arpack_si_extract( sigma, n, nev, ncv, tol, resid, vectors, ldv, &
        iparam, ipntr, workd, workl, lworkl, lambdareal, info ) bind( c )
    use, intrinsic :: iso_c_binding
    implicit none
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

    real( c_double ), intent( in )    :: sigma   ! the shift dseupd transforms back with
    integer( int_t ), intent( in )    :: n
    integer( int_t ), intent( in )    :: nev
    integer( int_t ), intent( in )    :: ncv
    real( c_double ), intent( inout ) :: tol
    integer( int_t ), intent( in )    :: ldv
    integer( int_t ), intent( in )    :: lworkl

    ! dsaupd's state, UNTOUCHED since the last step call -- dseupd requires
    ! that explicitly, and the C++ caller owns these arrays so nothing else
    ! can reach them
    real( c_double ), dimension( n ),        intent( inout ) :: resid
    real( c_double ), dimension( ldv, ncv ), intent( inout ) :: vectors
    integer( int_t ), dimension( 11 ),       intent( inout ) :: iparam
    integer( int_t ), dimension( 11 ),       intent( inout ) :: ipntr
    real( c_double ), dimension( 3 * n ),    intent( inout ) :: workd
    real( c_double ), dimension( lworkl ),   intent( inout ) :: workl

    real( c_double ), dimension( nev ), intent( inout ) :: lambdareal
    integer( int_t ), intent( inout ) :: info

    character, parameter :: bmat = 'I'
    character( len=2 ), parameter :: which = 'LM'
    character, parameter :: howmny = 'A'

    ! values only. rvec = .false. leaves z unreferenced, but ldz >= 1 is still
    ! required and select is still part of the ABI, so both exist as
    ! placeholders rather than being omitted
    logical, parameter :: rvec = .false.
    integer( int_t ), parameter :: ldz = 1

    logical, dimension( : ),    allocatable :: slct
    real( c_double ), dimension( :, : ), allocatable :: z
    real( c_double ), dimension( : ),    allocatable :: d

    integer( int_t ) :: i

    allocate( slct( ncv ) )
    allocate( z( ldz, 1 ) )

    ! dseupd fills nev entries ; sized at ncv for the same reason the other
    ! drivers oversize theirs -- it costs a few kB and removes a failure mode
    allocate( d( ncv ) )
    d = 0.0D0

    call dseupd( rvec, howmny, slct, d, z, ldz, sigma, &
                 bmat, n, which, nev, tol, resid, ncv, vectors, ldv, &
                 iparam, ipntr, workd, workl, lworkl, info )

    ! already transformed back to eigenvalues of A by dseupd
    do i = 1, nev
        lambdareal( i ) = d( i )
    end do

    deallocate( d )
    deallocate( z )
    deallocate( slct )

end subroutine arpack_si_extract
