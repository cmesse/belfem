
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

!> see also https://software.intel.com/content/www/us/en/develop/articles/pardiso-parameter-table.html

!> MKL's own Fortran interface for pardiso ( sixteen arguments, assumed-size arrays ), included the
!> way mumpstools.f90 includes dmumps_struc.h: the vendor header is the source of truth, and every
!> call below is checked against it at compile time. It declares default INTEGER throughout, which
!> is the LP64 MKL this project links; Intel's ILP64 story for Fortran is -fdefault-integer-8, which
!> BELFEM does not pass, so the ILP64 combination is refused below instead of compiling silently.
include 'mkl_pardiso.f90'

#ifdef BELFEM_INT64
#error "MKL PARDISO with USE_MKL_64BIT_API is not supported by pardisotools.f90: the vendor interface declares default INTEGER and would need -fdefault-integer-8"
#endif

module pardisotools
    use, intrinsic :: iso_c_binding
    use mkl_pardiso
#ifdef OMP
    use omp_lib
#endif
    implicit none

    !> make all globals public
    public
#ifdef BELFEM_INT64
    integer, parameter :: int_t = C_INT64_T
#else
    integer, parameter :: int_t = C_INT32_T
#endif

    !> Parameter list as fed to PARDISO
    integer, dimension( 64 ) :: gParameters

    !> the matrix type
    !>   1 real and structurally symmetric
    !>   2 real and symmetric positive defnite
    !>  -2 real and symmetric indefinite
    !>   3 complex and structurally symmetric
    !>   4 complex and Hermitian positive dfnite
    !>  -4 complex and Hermitian indfnite
    !>   6 complex and symmetric
    !>  11 real and nonsymmetric
    !>  13 complex and nonsymmetric
    integer :: gMatrixType

    !> verbosity flag
    integer :: gInfoLevel

    !>  memory pointers ( MKL's opaque handle, one 64-bit word each )
    type( MKL_PARDISO_HANDLE ), dimension( 64 ) :: gMemoryPointers

    integer :: gMaxNumFactors
    integer :: gNumFactors

    integer :: gN
    integer :: gNRHS

end module pardisotools

!------------------------------------------------------------------------------

!> List of Parameters
!> 1: Transposed Flag  : 0 - CSR
!>                     : 1 - CSC
!>
!> 2: Indexing Base    : 0 - C++
!>                       1 - Fortran
!>
!> 3: Matrix Type      : 11 - Unsymmetric,
!>                        2 - PositiveDefiniteSymmetric,
!>                       -2 - GeneralSymmetric
!> 4: Info Level       :  0 - Silent
!> 5: Precon Exponent
!> 6: Max Number of
!>    Refinement steps : 0 - auto
!> 7: Flag that selects the reduce ordering
!>
!>
!> 8: reserved ( was "compute determinant"; MKL PARDISO returns no determinant and
!>    documents iparm( 33 ) as reserved, so the value is read and ignored )

function pardisotools_initialize_parameters( aParameters ) bind( c ) result( aStatus )
    use pardisotools
    implicit none
    integer( int_t ), intent( in ),    dimension( 8 )         :: aParameters
    integer( int_t ) :: aStatus

    ! tells if the matrix has to be transposed, e.g. if this is a nonsymmetric
    ! CSC Matrix
    integer( int_t ) :: tTransposedFlag

    ! the indexing base 0: C++ Style, 1: Fortran Style ( is inverted later )
    integer( int_t ) :: tIndexingBase

    ! exponent for the preconditioning
    integer( int_t ) :: tPreconExponent

    ! Reduce ordering scheme
    integer( int_t ) :: tReduceOrdering

    ! refinement step limit
    integer( int_t ) :: tMaxNumRefinementSteps

    ! kept for the C interface ( parameter 8 ); MKL PARDISO has no determinant
    integer( int_t ) :: tComputeDeterminant

    ! iterator
    integer( int_t ) :: k

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! USER DEFINED PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! set the transposed flag
    tTransposedFlag          = aParameters( 1 )

    ! set the indexing base
    tIndexingBase            = aParameters( 2 )

    ! set the symmetry mode
    gMatrixType              = aParameters( 3 )

    ! set the verbosity flag
    gInfoLevel               = aParameters( 4 )

    ! set the preconditioning exponent
    tPreconExponent          = aParameters( 5 )

    ! set refinement step limit
    tMaxNumRefinementSteps   = aParameters( 6 )

    ! select the reduce ordering flag
    tReduceOrdering          = aParameters( 7 )

    ! set determinant flag
    tComputeDeterminant      = aParameters( 8 )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   SOLVER PARAMETERS
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! reset all parameters
    forall( k = 1:64 ) gParameters( k ) = 0

    ! the user sets all values
    gParameters(  1 ) = 1

    ! set the reduce ordering flag
    gParameters(  2 ) = tReduceOrdering

! NOT a mistake: OMP, not BELFEM_OMP. This sets PARDISO's OWN thread count, a
! third-party budget. BELFEM_OMP gates only BELFEM's own kernels; gating
! this on it would silently serialise PARDISO on every default build.
#ifdef OMP
    gParameters(3) = OMP_GET_MAX_THREADS()
#else
    gParameters(3) = 1
#endif

    ! preconditioning
    if( tPreconExponent .eq. 0 ) then
        ! no preconditioning
        gParameters( 4 ) = 0
    else
        ! for general matrices
        if( gMatrixType .eq. 2 ) then
            ! special case for positive definite symmetric
            gParameters( 4 ) = 10*tPreconExponent + 2
        else
            ! general matrix or symmetric but not positive definite
            gParameters( 4 ) = 10*tPreconExponent + 1
        end if
    end if

    gParameters( 8 ) = tMaxNumRefinementSteps

    ! Pertubation stuff
    if ( gMatrixType .eq. 11 ) then
        gParameters( 10 ) = 13 ! perturbe the pivot elements with 1E-13 (unsymmetric default)
        gParameters( 11 ) = 1 ! use nonsymmetric permutation and scaling MPS
    else
        gParameters( 10 ) = 8 ! perturbe the pivot elements with 1E-8 (symmetric default)
        gParameters( 11 ) = 0 ! disable scaling
    end if

    ! set the transposed flag
    gParameters( 12 ) = tTransposedFlag

    ! information output
    if( gInfoLevel .eq. 1 ) then
        gParameters( 18 ) = -1 ! Output: number of nonzeros in the factor LU
        gParameters( 19 ) = -1 ! Output: Mflops for LU factorization
    end if

#ifdef DEBUG
    ! Matrix checker
    gParameters( 27 ) = 1
#endif

    ! 33: reserved in MKL PARDISO ( must stay 0 ); the Panua/Schenk determinant flag does not exist here
    gParameters( 33 ) = 0

    ! 1-based/0-based input data indexing
    gParameters( 35 ) = 1 - tIndexingBase

    aStatus = 0

end function pardisotools_initialize_parameters

!------------------------------------------------------------------------------

function pardisotools_symbolic_factorization( &
        aN,        & ! size of matrix
        aNNZ,      & ! number of nonzeros
        aNRHS,     & ! number of RHS columns
        aPointers, & ! pointers of CSR / CSC matrix
        aIndices,  & ! indoces of CSR / CSC matrix
        aValues    & ! values of matrix
        ) bind( c ) result( aStatus )
    use pardisotools
    implicit none
    integer( int_t ), intent( in )                            :: aN
    integer( int_t ), intent( in )                            :: aNNZ
    integer( int_t ), intent( in )                            :: aNRHS
    integer( int_t ), intent( in ),    dimension( aN+1 )      :: aPointers
    integer( int_t ), intent( in ),    dimension( aNNZ )      :: aIndices
    real( c_double ), intent( in ),    dimension( aNNZ )      :: aValues
    integer( int_t )                                          :: aStatus

    ! iteration index
    integer( int_t ) :: k

    ! the phase of the current call
    integer( int_t ) :: tPhase

    ! placeholders for the arguments this phase does not use; one object per dummy argument
    ! that MKL declares INTENT( INOUT ) or INTENT( OUT ), so no actual is bound twice
    integer :: tPerm( 1 ) = 0
    real*8  :: tB( 1 ) = 0.0d0
    real*8  :: tX( 1 ) = 0.0d0

    !  Reordering and Symbolic Factorization, This step also allocates
    ! all memory that is necessary for the factorization
    tPhase = 11

    ! initialize factors
    gMaxNumFactors = 1
    gNumFactors    = 1

    !  Initiliaze the internal solver memory pointer.
    forall( k = 1:64 ) gMemoryPointers( k )%DUMMY = 0

    call pardiso ( &
            gMemoryPointers, &
            gMaxNumFactors, &
            gNumFactors, &
            gMatrixType, &
            tPhase, &
            aN, &
            aValues, &
            aPointers, &
            aIndices, &
            tPerm, &
            aNRHS, &
            gParameters, &
            gInfoLevel, &
            tB, &
            tX, &
            aStatus )

end function pardisotools_symbolic_factorization

!------------------------------------------------------------------------------

function pardisotools_solve( &
        aN,        & ! size of matrix
        aNNZ,      & ! number of nonzeros
        aNRHS,     & ! number of RHS columns
        aPointers, & ! pointers of CSR / CSC matrix
        aIndices,  & ! indoces of CSR / CSC matrix
        aValues,   & ! values of matrix
        aLHS,      & ! Left hand side
        aRHS,      & ! Right hand side
        aInfo      & ! debug information
        ) bind( c ) result( aStatus )
    use pardisotools
    implicit none
    integer( int_t ), intent( in )                            :: aN
    integer( int_t ), intent( in )                            :: aNNZ
    integer( int_t ), intent( in )                            :: aNRHS
    integer( int_t ), intent( in ),    dimension( aN+1 )      :: aPointers
    integer( int_t ), intent( in ),    dimension( aNNZ )      :: aIndices
    real( c_double ), intent( in ),    dimension( aNNZ )      :: aValues
    real( c_double ), intent( inout ), dimension( aN, aNRHS ) :: aLHS
    ! MKL declares the right-hand side INTENT( INOUT ) ( it is overwritten when iparm( 6 ) = 1,
    ! which this wrapper never sets ), so the Fortran intent has to say inout as well
    real( c_double ), intent( inout ), dimension( aN, aNRHS ) :: aRHS
    integer( int_t ), intent( out ) ,  dimension( 8 )         :: aInfo
    integer( int_t )                                          :: aStatus

    ! iteration index
    integer( int_t ) :: k

    ! the phase of the current call
    integer( int_t ) :: tPhase

    ! placeholders for the arguments the factorization does not use, one object per dummy
    integer :: tPerm( 1 ) = 0
    real*8  :: tB( 1 ) = 0.0d0
    real*8  :: tX( 1 ) = 0.0d0

    ! local copy of parameters
    integer, dimension( 64 ) :: tParameters


    ! create a local copy of the parameters
    forall( k=1:64 ) tParameters( k ) = gParameters( k )

    ! remember size
    gN = aN
    gNRHS = aNRHS

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   NUMERIC FACTORIZATION
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    !  Factorization.
    tPhase = 22 ! only factorization
    call pardiso ( &
            gMemoryPointers, &
            gMaxNumFactors, &
            gNumFactors, &
            gMatrixType, &
            tPhase, &
            aN, &
            aValues, &
            aPointers, &
            aIndices, &
            tPerm, &
            aNRHS, &
            tParameters, &
            gInfoLevel, &
            tB, &
            tX, &
            aStatus )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   SOLVING
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! check of factorization was successful
    if( aStatus .eq. 0 ) then

        !  Back substitution and iterative refinement
        tPhase = 33

        call pardiso ( &
                gMemoryPointers, &
                gMaxNumFactors, &
                gNumFactors, &
                gMatrixType, &
                tPhase, &
                aN, &
                aValues, &
                aPointers, &
                aIndices, &
                tPerm, &
                aNRHS, &
                tParameters, &
                gInfoLevel, &
                aRHS, &
                aLHS, &
                aStatus )
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   COPY RETURN PARAMETERS INTO INFO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    aInfo( 1 ) = tPhase

    ! Number of performed iterative refinement steps
    aInfo( 2 ) = tParameters( 7 )

    ! Output: number of nonzeros in the factor LU
    aInfo( 3 ) = tParameters( 18 )

    ! Output: Mflops for LU factorization
    aInfo( 4 ) = tParameters( 19 )

    ! CGS diagnostic
    aInfo( 5 ) = tParameters( 20 )

    ! Number of positive eigenvalues
    aInfo( 6 ) = tParameters( 22 )

    ! Number of negative eigenvalues
    aInfo( 7 ) = tParameters( 23 )

    ! iparm( 33 ) is reserved in MKL PARDISO and stays 0; slot 8 is kept for the C interface
    aInfo( 8 ) = tParameters( 33 )

end function pardisotools_solve

!------------------------------------------------------------------------------

function pardisotools_free() bind( c ) result( aStatus )
    use pardisotools
    implicit none

    integer( int_t ) :: aStatus

    ! placeholders: the release phase reads none of them, but MKL declares perm, b and x as
    ! INTENT( INOUT ) / INTENT( OUT ), so each gets its own object
    real*8  :: tA( 1 ) = 0.0d0
    integer :: tIA( 1 ) = 0
    integer :: tJA( 1 ) = 0
    integer :: tPerm( 1 ) = 0
    real*8  :: tB( 1 ) = 0.0d0
    real*8  :: tX( 1 ) = 0.0d0
    integer :: tPhase = -1

    call pardiso ( &
            gMemoryPointers, &
            gMaxNumFactors, &
            gNumFactors, &
            gMatrixType, &
            tPhase, &
            gN, &
            tA, &
            tIA, &
            tJA, &
            tPerm, &
            gNRHS, &
            gParameters, &
            gInfoLevel, &
            tB, &
            tX, &
            aStatus )
    
end function pardisotools_free
