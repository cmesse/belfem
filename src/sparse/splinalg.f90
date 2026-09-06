/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */
 
subroutine matvec_csr( n, m, nnz, values, indices, pointers, x, y, base ) bind ( c )
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
    integer( int_t ), intent( in ) :: n ! number of rows in matrix
    integer( int_t ), intent( in ) :: m ! number of columns in matrix
    integer( int_t ), intent( in ) :: nnz ! number of nonzeros
    real( c_double ), dimension( nnz ), intent( in ) :: values
    integer( int_t ), dimension( nnz ), intent( in  ) :: indices    ! column indices
    integer( int_t ), dimension( n + 1 ), intent( in  ) :: pointers ! row pointers
    real( c_double ), dimension( m ), intent( in ) :: x    ! left hand side
    real( c_double ), dimension( n ), intent( inout ) :: y ! right hand side

    ! indexing base of pointers and indices: 0 for C, 1 for Fortran. The
    ! caller passes pointers( 1 ) so the kernel works in either base and
    ! the matrix never has to be rewritten for a matvec.
    integer( int_t ), intent( in ) :: base

    integer( int_t ) i, j, a, b, off
    real( c_double ) :: yi

    ! shift from the caller's base to Fortran's 1-based array access.
    ! The slice upper bound needs off - 1, which is just -base, so it is
    ! spelled that way below rather than computing ( - 1 + off ).
    off = 1 - base

    ! reset values
    y = 0.0d0

#ifdef BELFEM_OMP
!$omp parallel do &
!$omp shared(n,m,nnz,x,y,pointers,values,indices,off,base) &
!$omp schedule(static) &
!$omp private(a,b,yi)
#endif
    do i = 1, n
        a = pointers(i) + off
        b = pointers(i+1) - base
        yi = 0.0d0
        do j = a, b
            yi = yi + values(j) * x(indices(j) + off)
        end do
        y(i) = yi
    end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

end subroutine matvec_csr

subroutine matvec_csc( n, m, nnz, values, indices, pointers, x, y, base ) bind ( c )
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
    integer( int_t ), intent( in ) :: n ! number of rows in matrix
    integer( int_t ), intent( in ) :: m ! number of columns in matrix
    integer( int_t ), intent( in ) :: nnz ! number of nonzeros
    real( c_double ), dimension( nnz ), intent( in ) :: values
    integer( int_t ), dimension( nnz ), intent( in  ) :: indices    ! row indices
    integer( int_t ), dimension( m + 1 ), intent( in  ) :: pointers ! column pointers
    real( c_double ), dimension( m ), intent( in ) :: x ! left hand side
    real( c_double ), dimension( n ), intent( inout ) :: y ! right hand side

    ! indexing base of pointers and indices: 0 for C, 1 for Fortran.
    ! See the note in matvec_csr.
    integer( int_t ), intent( in ) :: base

    integer( int_t ) i, j, a, b, off
    real( c_double ) :: xj

    ! shift from the caller's base to Fortran's 1-based array access.
    ! The slice upper bound needs off - 1, which is just -base, so it is
    ! spelled that way below rather than computing ( - 1 + off ).
    off = 1 - base

    ! reset values
    y = 0.0d0

    ! CSC is a scatter, so the column loop collides on y. An array
    ! reduction gives each thread a private y and combines once at the
    ! end, instead of one atomic per nonzero, which serialised the inner
    ! loop at any real thread count. The cost is n * 8 bytes per thread of
    ! scratch; if that ever becomes the binding constraint, the fallback
    ! is the transposed-CSR duality path, which gathers instead.
    !
    ! Those n * 8 bytes are STACK. y is an explicit-shape
    ! dummy, so gfortran materialises each thread's private copy in that
    ! thread's frame, and an OpenMP worker gets 2 MiB on Darwin -- so this
    ! crashed with SIGBUS above n = 262,144. The directives are now gated on
    ! BELFEM_OMP, which is OFF by default: this matvec runs master-only and is
    ! noise beside assembly and factorization, so the threading never paid for
    ! the exposure. Turning BELFEM_OMP back on re-arms the defect unless the
    ! private copy is moved off the stack first. OMP_STACKSIZE raises the
    ! worker stack if you need it. See doc/parallel_execution.md.
#ifdef BELFEM_OMP
!$omp parallel do &
!$omp shared(n,m,nnz,x,pointers,values,indices,off,base) &
!$omp private(a,b,xj) &
!$omp reduction(+:y) &
!$omp schedule(static)
#endif
    do j = 1, m
        a  = pointers(j) + off
        b =  pointers(j+1) - base
        xj = x(j)
        do i = a, b
            y(indices(i) + off) = y(indices(i) + off) + values(i) * xj
        end do
    end do
#ifdef BELFEM_OMP
!$omp end parallel do
#endif

end subroutine matvec_csc