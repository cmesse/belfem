
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

!> see also https://software.intel.com/content/www/us/en/develop/articles/pardiso-parameter-table.html

module pardisotools
    use, intrinsic :: iso_c_binding
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

    !>  memory pointers
    integer*8, dimension( 64 ) :: gMemoryPointers

    !> output variables
    double precision, dimension( 64 ) :: gDPARM

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
!> 8: Compute Determinant : 0 - off
!>                          1 - on

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

    ! flag that tells if the determinant is to be computer
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

    ! 33: Determinant of a matrix.
    gParameters( 33 ) = tComputeDeterminant

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

    ! some dummy values
    integer( int_t ) :: tIntDummy
    real*8  :: tRealDummy

    !  Reordering and Symbolic Factorization, This step also allocates
    ! all memory that is necessary for the factorization
    tPhase = 11

    ! initialize factors
    gMaxNumFactors = 1
    gNumFactors    = 1

    !  Initiliaze the internal solver memory pointer.
    forall( k = 1:64 ) gMemoryPointers( k ) = 0

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
            tIntDummy, &
            aNRHS, &
            gParameters, &
            gInfoLevel, &
            tRealDummy, &
            tRealDummy, &
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
    real( c_double ), intent( in ),    dimension( aN, aNRHS ) :: aRHS
    integer( int_t ), intent( out ) ,  dimension( 8 )         :: aInfo
    integer( int_t )                                          :: aStatus

    ! iteration index
    integer( int_t ) :: k

    ! the phase of the current call
    integer( int_t ) :: tPhase

    ! some dummy values
    integer( int_t ) :: tIntDummy
    real*8  :: tRealDummy

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
            tIntDummy, &
            aNRHS, &
            tParameters, &
            gInfoLevel, &
            tRealDummy, &
            tRealDummy, &
            aStatus )

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !   SOLVING
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! check of factorization was successful
    if( aStatus .eq. 0 ) then

        ! initialize parameters
        forall( k=1:64 ) gDPARM( k ) = 0.0

        !  Back substitution and iterative refinement
        tPhase = 33 ! only factorization

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
                tIntDummy, &
                aNRHS, &
                tParameters, &
                gInfoLevel, &
                aRHS, &
                aLHS, &
                aStatus, &
                gDPARM )
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

    ! compute-determinant flag ( iparm 33 ); the determinant is in gDPARM( 33 )
    aInfo( 8 ) = tParameters( 33 )

end function pardisotools_solve

!------------------------------------------------------------------------------

function pardisotools_free() bind( c ) result( aStatus )
    use pardisotools
    implicit none

    integer( int_t ) :: aStatus

    ! dummy values
    integer( int_t ) :: tIntDummy
    real*8  :: tRealDummy
    integer :: tPhase = -1

    call pardiso ( &
            gMemoryPointers, &
            gMaxNumFactors, &
            gNumFactors, &
            gMatrixType, &
            tPhase, &
            gN, &
            tRealDummy, &
            tIntDummy, &
            tIntDummy, &
            tIntDummy, &
            gNRHS, &
            gParameters, &
            gInfoLevel, &
            tRealDummy, &
            tRealDummy, &
            aStatus )
    
end function pardisotools_free

!------------------------------------------------------------------------------

function pardisotools_get_determinant() bind( c ) result( aDet )
    use pardisotools
    implicit none
    real( c_double ) :: aDet

    aDet = gDPARM( 33 )

end function pardisotools_get_determinant