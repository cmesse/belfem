
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

module intpoints
    use, intrinsic :: iso_c_binding
    implicit none
    private
!   the high percision datatype
    integer, parameter :: quadruple = selected_real_kind(r=60)
    integer, parameter :: QP = quadruple

    real(quadruple), private, parameter :: eps = epsilon( real(1D0,quadruple) )
    real(quadruple), private, parameter :: pi = 2.0d0 * acos( 0.0D0 )
    public :: intpoints_gauss
    public :: intpoints_gauss_quadruple
    public :: intpoints_lobatto
    public :: quadruple
    public :: QP
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

    !>  evaluates the legendre polynomial
    subroutine legendre(n, x, P, dPdx, d2Pdx2)
        integer, intent(in) :: n
        real(quadruple), intent(in) :: x
        real(quadruple), intent(out) :: P
        real(quadruple), intent(out), optional :: dPdx
        real(quadruple), intent(out), optional :: d2Pdx2

        integer :: i
        real(quadruple) :: f0, f1, f2, g0, g1, h0, c, k, m

        select case(n)
        case(1)
            f0 = 0._QP
            g0 = 0._QP
            h0 = 0._QP
        case(2)
            f0 = 1._QP
            g0 = 0._QP
            h0 = 0._QP
        case default
            f1 = x                      ! P0
            f0 = 0.5_QP*(3._QP*x**2-1._QP) ! P1

            g0 = 3._QP*x                 ! dP1dx
            h0 = 2._QP                   ! d2P1dx2

            if(n .gt. 2) then

                do i = 3, n
                    k = real(i, quadruple)
                    f2 = f1
                    f1 = f0
                    f0 = (((2._QP*k-1._QP)*x)*f1-(k-1._QP)*f2)/k
                end do


                if(present(dPdx)) then
                    m = real(n, quadruple)
                    c = 1._QP/(x**2-1._QP)
                    g1 = (m-1._QP)*(x*f1-f2)*c
                    g0 = m*(x*f0-f1)*c
                    if(present(d2Pdx2)) h0 = (m*f0+x*(m-2._QP)*g0-n*g1)*c
                end if
            end if
        end select
        P = f0
        if(present(dPdx)) then
            dPdx = g0
            if(present(d2Pdx2)) d2Pdx2 = h0
        end if
    end subroutine legendre

!-------------------------------------------------------------------------------

    subroutine intpoints_gauss(n, w, xi) bind( c )
        !> number of integration points
        integer( c_int ), intent(in) :: n

        !> weighting factors
        real( c_double ), dimension( n ), intent(inout) :: w

        !> integration points
        real( c_double ), dimension( n ), intent(inout) :: xi

        integer :: k, m
        real(quadruple) :: xi1, xi0, P, dPdxi

        w = 0D0
        xi = 0D0
        select case(n)
        case(1)
            w(1)  =  2D0
            xi(1) =  0D0
        case(2)
            w(1)  =  1D0
            w(2)  =  1D0
            xi(1) = -dsqrt(1D0/3D0)
            xi(2) =  dsqrt(1D0/3D0)
        case(3)
            w(1)  =  5D0/9D0
            w(2)  =  8D0/9D0
            w(3)  =  5D0/9D0
            xi(1) = -dsqrt(0.6D0)
            xi(2) =  0D0
            xi(3) =  dsqrt(0.6D0)
        case default
            m = n/2
            do k = 1, m
                xi1 = real(cos(pi*real(4*k-1,quadruple)/real(4*n+2,quadruple)),quadruple)
                xi0 = -2._QP
                do while(abs(xi0-xi1) .gt. eps)
                    call legendre(n, xi1, P, dPdxi)
                    xi0 = xi1
                    xi1 = xi1 - P/dPdxi
                end do
                xi(k) = -dble(xi1)
                w(k) = dble(2._QP/((1._QP-xi1**2)*dPdxi**2))
                xi(n-k+1) = dble(xi1)
                w(n-k+1) = w(k)
            end do
            if(mod(n, 2) .eq. 1) w(m+1) = 2D0-sum(w)
        end select
    end subroutine intpoints_gauss

!-------------------------------------------------------------------------------

    subroutine intpoints_gauss_quadruple(n, w, xi)
        !> number of integration points
        integer( c_int ), intent(in) :: n

        !> weighting factors
        real( quadruple ), dimension( n ), intent(inout) :: w

        !> integration points
        real( quadruple ), dimension( n ), intent(inout) :: xi

        integer :: k, m
        real(quadruple) :: xi1, xi0, P, dPdxi

        w = 0_QP
        xi = 0_QP
        select case(n)
        case(1)
            w(1)  =  2_QP
            xi(1) =  0_QP
        case(2)
            w(1)  =  1_QP
            w(2)  =  1_QP
            xi(1) = -sqrt(1._QP/3._QP)
            xi(2) =  sqrt(1._QP/3._QP)
        case(3)
            w(1)  =  5_QP/9_QP
            w(2)  =  8_QP/9_QP
            w(3)  =  5_QP/9_QP
            xi(1) = -sqrt(0.6_QP)
            xi(2) =  0_QP
            xi(3) =  sqrt(0.6_QP)
        case default
            m = n/2
            do k = 1, m
                xi1 = real(cos(pi*real(4*k-1,quadruple)/real(4*n+2,quadruple)),quadruple)
                xi0 = -2_QP
                do while(abs(xi0-xi1) .gt. eps)
                    call legendre(n, xi1, P, dPdxi)
                    xi0 = xi1
                    xi1 = xi1 - P/dPdxi
                end do
                xi(k) = -xi1
                w(k) = 2._QP/((1._QP-xi1**2)*dPdxi**2)
                xi(n-k+1) = xi1
                w(n-k+1) = w(k)
            end do
            if(mod(n, 2) .eq. 1) w(m+1) = 2._QP-sum(w)
        end select

    end subroutine intpoints_gauss_quadruple

!-------------------------------------------------------------------------------

    subroutine intpoints_lobatto( n, w, xi ) bind( c )
        !> number of integration points
        integer( c_int ), intent(in) :: n

        !> weighting factors
        real( c_double ), dimension( n ), intent(inout) :: w

        !> integration points
        real( c_double ), dimension( n ), intent(inout) :: xi

        integer :: k, m, i
        real(quadruple) :: xi0, xi1, z0, z1, P, dPdxi, d2Pdxi2

        w = 0D0
        xi = 0D0
        xi(1) = -1D0
        xi(n) = 1D0
        w(1) = 2D0/dble(n*(n-1))
        w(n) = w(1)

        if(n .gt. 3) then
            m = n/2
            z1 = real(-cos(pi*3D0/dble(4*n-2)),quadruple)
            do k = 2, m
                xi0 = -2._QP
                z0 = z1
                z1 = real(-cos(pi*dble(4*k-1)/dble(4*n-2)),quadruple)
                xi1 = 0.5_QP*(z0+z1)
                i = 0
                do while(abs(xi0-xi1) .gt. eps)
                    i = i + 1
                    call legendre(n-1, xi1, P, dPdxi, d2Pdxi2)
                    xi0 = xi1
                    xi1 = xi1 - dPdxi/d2Pdxi2
                end do
                xi(k) = dble(xi1)
                xi(n-k+1) = -xi(k)
                w(k) = 2D0/(dble(real(n*(n-1),quadruple)*P**2))
                w(n-k+1) = w(k)
            end do
            if(mod(n, 2) .eq. 1) w(m+1) = 2D0-sum(w)
        else
            w(2) = 4D0/3D0
        end if
    end subroutine intpoints_lobatto
!-------------------------------------------------------------------------------
end module intpoints