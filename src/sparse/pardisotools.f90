
!> see also https://software.intel.com/content/www/us/en/develop/articles/pardiso-parameter-table.html

module pardisotools
    use, intrinsic :: iso_c_binding
    implicit none

    !> make all globals public
    public

    !> Parameter list as fed to PARDISO
    integer, dimension( 64 ) :: gParameters

    !> the matrix type
    !>   1 real and structurally symmetric
    !>   2 real and symmetric positive denite
    !>  -2 real and symmetric indefinite
    !>   3 complex and structurally symmetric
    !>   4 complex and Hermitian positive denite
    !>  -4 complex and Hermitian indenite
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
    integer( c_int ), intent( in ),    dimension( 8 )         :: aParameters
    integer( c_int ) :: aStatus

    ! tells if the matrix has to be transposed, e.g. if this is a nonsymmetric
    ! CSC Matrix
    integer :: tTransposedFlag

    ! the indexing base 0: C++ Style, 1: Fortran Style ( is inverted later )
    integer :: tIndexingBase

    ! exponent for the preconditioning
    integer :: tPreconExponent

    ! Reduce ordering scheme
    integer :: tReduceOrdering

    ! refinement step limit
    integer :: tMaxNumRefinementSteps

    ! flag that tells if the determinant is to be computer
    integer :: tComputeDeterminant

    ! iterator
    integer :: k

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

#ifdef BELFEM_OPENMP
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

#ifdef BELFEM_DEBUG
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
    integer( c_int ), intent( in )                            :: aN
    integer( c_int ), intent( in )                            :: aNNZ
    integer( c_int ), intent( in )                            :: aNRHS
    integer( c_int ), intent( in ),    dimension( aN+1 )      :: aPointers
    integer( c_int ), intent( in ),    dimension( aNNZ )      :: aIndices
    real( c_double ), intent( in ),    dimension( aNNZ )      :: aValues
    integer( c_int )                                          :: aStatus

    ! iteration index
    integer :: k

    ! the phase of the current call
    integer :: tPhase

    ! some dummy values
    integer :: tIntDummy
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
    integer( c_int ), intent( in )                            :: aN
    integer( c_int ), intent( in )                            :: aNNZ
    integer( c_int ), intent( in )                            :: aNRHS
    integer( c_int ), intent( in ),    dimension( aN+1 )      :: aPointers
    integer( c_int ), intent( in ),    dimension( aNNZ )      :: aIndices
    real( c_double ), intent( in ),    dimension( aNNZ )      :: aValues
    real( c_double ), intent( inout ), dimension( aN, aNRHS ) :: aLHS
    real( c_double ), intent( in ),    dimension( aN, aNRHS ) :: aRHS
    integer( c_int ), intent( out ) ,  dimension( 8 )         :: aInfo
    integer( c_int )                                          :: aStatus

    ! iteration index
    integer :: k

    ! the phase of the current call
    integer :: tPhase

    ! some dummy values
    integer :: tIntDummy
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

    ! natural log of determinant ( if computed )
    aInfo( 8 ) = tParameters( 33 )

end function pardisotools_solve

!------------------------------------------------------------------------------

function pardisotools_free() bind( c ) result( aStatus )
    use pardisotools
    implicit none

    integer( c_int ) :: aStatus

    ! dummy values
    integer :: tIntDummy
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