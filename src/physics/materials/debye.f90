
! BELFEM -- The Berkeley Lab Finite Element Framework
! Copyright (c) 2026, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any required
! approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Developers: Christian Messe, Gregory Giard
!
! See the top-level LICENSE file for the complete license and disclaimer.

module debye
    use, intrinsic :: iso_c_binding
    use intpoints, only: intpoints_gauss_quadruple, quadruple, QP
    implicit none
    private
    real( QP ), parameter :: pi = 3.141592653589793238462643383279502884_QP
    real( QP ), parameter :: kB = 1.380649e-23_QP
    real( QP ), parameter :: NA = 6.02214076e23_QP
    real( QP ), parameter :: h = 6.62607015e-34_QP
    real( QP ), parameter :: hbar = 0.5_QP*h/pi
    real( QP ), parameter :: el = 1.602176634e-19_QP
    real( QP ), parameter :: me = 9.1093837e-31_QP
    real( QP ), parameter :: eV =  1.602176634e-19_QP
    real( QP ), parameter :: c =  299792458.0


    ! integration points, scaled from 0 to 1
    real( QP ), dimension( : ), allocatable :: xi
    real( QP ), dimension( : ), allocatable :: eta

    ! integration weights, so that Σw = 1
    real( QP ), dimension( : ), allocatable :: w
    public :: debye_table, callaway_conductivity

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

    real( quadruple ) function debye_function( n, z ) result ( f )
        real( c_double ), intent( in ) :: n
        real( quadruple ), intent( in ) :: z
        real( quadruple ) :: exp_z

        exp_z = exp( z )

        f = z**n * exp_z / ( ( exp_z - 1_QP ) * ( exp_z - 1_QP ) )

    end function debye_function


    subroutine debye_derivatives( n, z, y, dy )
        real( c_double ), intent( in ) :: n
        real( quadruple ), intent( in )  :: z
        real( quadruple ), intent( out ) :: y, dy
        real( quadruple ) :: u, v, f, g
        real( quadruple ) :: du, dv, df, dg
        real( quadruple ) :: ddu, ddv, ddf, ddg

        u = z**n
        v = exp( z )
        f = u * v
        g = (v-1._QP)*(v-1._QP)

        du = n * u / z
        dv = v
        df = du * v + dv*u
        dg = 2._QP * v * (v - 1._QP)

        ddu = du * ( n - 1._QP ) / z
        ddv = v

        ddf = ddu * v + 2._QP * du * dv + ddv * u
        ddg = 2._QP*v*(v+v-1._QP)

        y = ( g * df - f * dg ) / ( g * g )
        dy = ( g * ( ddf- f * ddg - 2._QP * df * dg ) + 2._QP * f * dg * dg ) / ( g * g * g )

    end subroutine debye_derivatives

    real( quadruple ) function debye_zpeak(n) result ( z )
        real( c_double ), intent( in ) :: n
        real( quadruple ) :: f, df

        real( quadruple ) :: g, f0
        integer :: k

        ! initial guess

        if ( n < 3. ) then
            write( *,* ) "invalid value for n"
            stop
        end if

        z = 1.1279 * n - 0.7499

        k = 0
        f = 1._QP
        f0 = 0._QP

        do while ( abs(f-f0) > 1D-16 )
            f0 =  f
            call debye_derivatives( n, z, f, df )
            z = z - f/df
            k = k + 1
            if ( k == 1000 ) then
                write( *,* ) "error in debye_zpeak: too many iterations"
                stop
            end if
        end do
    end function debye_zpeak

    subroutine debye_table( n, m, p, q, z, y ) bind( c )

        !> exponent value for integral
        real( c_double ),  intent( in )  :: n ! set to 3, 4 or 5

        !> number of integration points per inteval
        integer( c_int ),  intent( in )  :: m

        !> number of points to peak value
        integer( c_int ),  intent( in )  :: p ! set to 20

        !> total number of points, n=3: 330,  n=4: 237, n=5: 196
        integer( c_int ),  intent( in )  :: q

        !> output table of abscissae
        real(  c_double ), dimension( q ), intent( out ) :: z

        !> output table of function values
        real(  c_double ), dimension( q ), intent( out ) :: y

        ! temporary value for zk
        real( quadruple ) :: zi, z0, z1, yk, dz, y1

        integer :: k, i

        ! compute maximum value
        zi = debye_zpeak( n )


        ! calculate stepwidth
        dz = zi / ( p - 1 )

        call allocate_intpoints( m )

        z1 = 0.0_QP
        y1 = 0.0_QP

        z( 1 ) = dble( z1 )
        y( 1 ) = dble( y1 )

        ! loop over all integrals
        do k = 2, q
            z0 = z1
            z1 = z1 + dz

            ! integrate
            yk = 0.0_QP
            do i = 1, m
                zi = xi( i ) * z0 + eta( i ) * z1
                yk = yk + w( i ) * debye_function( n, zi )
            end do
            y1 = y1 + dz * yk

            z( k ) = dble( z1 )
            y( k ) = dble( y1 )
        end do

    end subroutine debye_table

    real( quadruple ) function scatter_function( n, z, K ) result ( f )
        real( quadruple ), intent( in ) :: n
        real( quadruple ), intent( in ) :: z
        real( quadruple ), dimension( 5 ), intent( in ) :: K
        real( quadruple ) :: exp_z, omega, iTauU, iTauM, iTauB, iTauPhE, tau

        ! frequency
        omega   = K( 1 ) * z

        ! Umklapp scattering
        iTauU   = K( 2 ) * omega * omega

        ! mass difference impurity
        iTauM   = K( 3 ) * omega * omega  * omega * omega

        ! boundary scattering
        iTauB   = K( 4 )

        ! phonon-electron scattering
        iTauPhE = K( 5 ) * omega

        ! Matthiessen's rule
        tau = 1.0_QP/( iTauU + iTauM + iTauB + iTauPhE )

        exp_z = exp( z )

        f = tau * z**n * exp_z / ( ( exp_z - 1_QP ) * ( exp_z - 1_QP ) )

    end function scatter_function

    ! see J. Callaway, Physical Review 113, 1046 (1959) –
    ! "Model for Lattice Thermal Conductivity at Low Temperatures"

    ! tau computation see  Zou and Balandin, 2001,
    ! "Phonon heat conduction in a semiconductor nanowire"

    subroutine callaway_conductivity( params, result, status ) bind( c )
        real( c_double ), dimension(20), intent( in ) :: params
        real( c_double ), intent( out ) :: result
        integer, intent( out ) :: status

        integer, parameter :: numintpoints = 7
        integer, parameter :: numintervals_callaway_equation = 100
        integer, parameter :: numintervals_to_ered_equals_one = 10
        integer            :: numintervals_density_equation
        real( quadruple ), parameter :: E_cutoff = 10.0_QP

        ! user defined parameters, see description below
        real( quadruple ) :: T, theta, vg, rho, G, gruen, M, Gamma, Tc, L0
        real( quadruple ) :: n, b, d, eps, meeff, ne, Delta0, lambda_pe

        real( quadruple ) :: lambda_opt, omega_opt
        real( quadruple ) :: V0, omegaD, A, Nratio, f, U, V

        real( quadruple ) :: E, Ered, dE, Ered0, Ered1, dEred, Eredmax, Emax

        real( quadruple ) :: Z, Z0, Z1, dZ

        real( quadruple ) :: Y, yj

        ! these are constants we precompute for the integration
        real( quadruple ), dimension(5) :: K

        integer :: i, j
        logical :: compute_density_function = .false.

        ! - - - - - - - - - - - - - - - - - - - -
        ! unpack the parameters
        ! - - - - - - - - - - - - - - - - - - - -

        ! first the  variables computed by the material model
        ! they may or may not depend on the temperature
        T       = params( 1 )        ! temperature in K


        theta   = params( 2 )        ! debye temperature in K
        vg      = params( 3 )        ! group velocity in m/s

        rho     = params( 4 )        ! density in kg/m³, used for tau_ph-e

        G       = params( 5 )        ! Shear Modulus in Pa, used for tau_U
        gruen   = params( 6 )        ! Grüneisen parameter, dimensionless, used for tau_U


        ! fixed material properties based on composition
        M       = params(  7 )        ! molar mass in kg/mol
        Gamma   = params(  8 )        ! impurity parameter, used for tau_M
        Tc      = params(  9 )        ! critical temperature in K

        ! this is actually the thickness of the REBCO Layer
        L0      = params( 10 )       ! thickness of the Rebco layer in m, used for tau_B

        ! now the model parameters
        n       = params( 11 )       ! exponent for the function, usually n=4


        ! parameters for Umklapp scattering, ignored if b=0
        ! C++ index | Fortran index | meaning
        !    11     | params(12)    | b, Umklapp parameter
        !    12     | params(13)    | d, Umklapp temperature correction
        b       = params( 12 )
        d       = params( 13 )

        eps     = params( 14 )*eV    ! deformation potential, passed in eV, used for tau_ph-e
        meeff   = params( 15 )*me    ! factor for effective electron mass, used for tau_ph-e
        ne      = params( 16 )       !  conduction electrons concentration in 1/ m³

        Delta0  = params( 17 )*kB*Tc ! factor for d-wave coupling, set to zero if no SC reduction

        lambda_pe    = params( 18 )   ! value for acoustic electron-phonon coupling

        ! values for optical coupling
        ! C++ index | Fortran index | meaning
        !    18     | params(19)    | lambda_opt, optical e-ph coupling
        !    19     | params(20)    | omega_opt, Raman shift in 1/cm
        lambda_opt = params( 19 )   ! value for optical electron-phonon coupling

        omega_opt  = 100_QP*params( 20 )*h*c ! frequency for Oxygen band peak, value passed in 1/cm

        ! initialize if necessary
        call allocate_intpoints( numintpoints )

        status = 0

        ! we can't use both optical and acoustic coupling, this will result in an error
        if ( ( lambda_opt > 0.0_QP ) .and. ( lambda_pe > 0.0_QP ) ) then
            status = 1
            result = 0.0_QP
            return
        end if

        ! constant to compute the debye frequency
        K(1) = kB * T / hbar ! * z = omega

        ! constant for phonon-phonon scattering
        V0 = M / ( rho * NA )          ! volume per atom
        omegaD = theta * kB / hbar

        K(2) = 2_QP * gruen * gruen * kB * T / ( G * V0 * omegaD ) ! * omega^2 = 1/tau_U

        if ( b /= 0.0_QP ) K( 2 ) = K( 2 ) * exp( -theta / ( ( b * ( 1._QP + d * T / theta ) * T ) ) )

        ! constant for mass-difference impuriy scattering
        K(3) = V0 * Gamma / ( 4_QP * pi * vg * vg * vg ) ! * omega^4 = 1/tau_M

        ! constant for boundary scattering
        K(4) = vg / L0 ! = 1/tau_B

        K( 5 ) = 0.0_QP

        if ( lambda_opt > 0.0_QP ) then ! use optical model

            ! omega_opt already holds the phonon energy h*c*nu ( = hbar*omega ),
            ! see the conversion above; the earlier extra factor hbar made U
            ! ~ 1e-32 and switched the optical channel off
            U = omega_opt / ( kB * T )
            V = exp( U )

            A = ne * kB * T / ( rho * vg * vg ) * ( U * U * V ) /( V - 1.0_QP )

            K(5) = A * lambda_opt

            compute_density_function = .true.
        end if

        if ( lambda_pe > 0.0_QP ) then  ! use acoustic model

            ! base_prefactor
            A= ne * eps * eps / ( rho * vg * vg * kB * T )

            K(5) = A * lambda_pe

            compute_density_function = .true.
        end if

        if ( ( T < Tc ) .and. ( Delta0 > 0.0_QP ) .and. compute_density_function ) then
            Y = 0.0_QP

            dEred = 1.0_QP / numintervals_to_ered_equals_one

            dE = dEred  *  Delta0

            Ered1 = 0.0_QP

            Emax = E_cutoff * kB * T

            Eredmax = Emax / Delta0

            numintervals_density_equation = ceiling( Eredmax / dEred )

            do j = 1, numintervals_density_equation
                Ered0 = Ered1
                Ered1 = Ered1 + dEred

                yj = 0.0_QP
                do i = 1, numintpoints
                    ! energy at point i
                    Ered = xi( i ) * Ered0 + eta ( i ) * Ered1

                    ! unreduced energy
                    E = Ered * Delta0

                    ! adding  d-wave density of states ratio
                    if ( Ered < 1.0_QP ) then
                        ! Subgap: pure nodal linear behavior
                        Nratio = 2.0_QP/pi * Ered
                    else
                        ! Above gap: BCS-like with sqrt divergence
                        ! since we chose dEred wisely,
                        ! we know that dividing by the square root below
                        ! will never cause a numerical issue!
                        Nratio = Ered / sqrt( Ered * Ered - 1.0_QP )
                    end if

                    ! Fermi factor
                    f = 1.0_QP / (1.0_QP + exp( E / ( kB*T ) ) )

                    yj = yj + w( i ) * Nratio * f * ( 1._QP - f )
                end do
                Y = Y + yj * dE
            end do
            K( 5 ) = K( 5 ) * Y / ( kB * T )
        else if ( .not. compute_density_function ) then
            ! constant for phonon-electron scattering

            if ( ( ne > 0._QP) .and. (eps > 0._QP) .and. ( meeff > 0._QP)) then
                A = ne * eps * eps / ( rho * vg * vg * kB * T )
                U = meeff * vg * vg / ( 2._QP * kB * T )
                V = pi * U
                K(5) = A * sqrt( V ) / exp( U ) ! * omega = 1/tau_ph-e

                ! superconducting reduction
                if( Delta0 /= 0._QP ) K(5) = K(5)*2._QP / exp( Delta0 / ( kB * T ) )
            else
                K( 5 ) = 0._QP
            end if
        end if


        Z = theta / T
        dZ = Z / numintervals_callaway_equation

        Y = 0.0_QP
        Z1 = 0.0_QP

        do j = 1, numintervals_callaway_equation
            Z0 = Z1
            Z1 = Z1 + dZ
            yj = 0.0_QP
            do i = 1, numintpoints
                Z = xi( i ) * Z0 + eta( i ) * Z1
                f = scatter_function( n, Z, K )
                yj = yj + w(i) * f
            end do
            Y = Y + yj * dZ
        end do

        U = kB/(2._QP * pi * pi * vg )
        V = kB/hbar*T

        result = Y*U*V*V*V
    end subroutine callaway_conductivity

    subroutine allocate_intpoints( n )
        integer, intent( in ) :: n
        integer :: i
        logical :: init

        init = .false.


        if( allocated( xi ) ) then
            if( size( xi ) /= n ) then
                init = .true.
                deallocate( w )
                deallocate( xi )
                deallocate( eta )
            end if
        else
           init = .true.
        end if

        if ( init ) then

            ! allocate memory
            allocate( w( n ) )
            allocate( xi( n ) )
            allocate( eta( n ) )

            ! populate data
            call intpoints_gauss_quadruple( n, w, xi )

            ! scale weights to 1
            forall ( i = 1:n ) eta( i ) = 0.5_QP * xi( i ) + 0.5_QP
            forall ( i = 1:n ) xi( i ) = 1.0_QP - eta( i )
            forall( i = 1:n )  w( i ) = 0.5_QP * w ( i )
        end if

    end subroutine allocate_intpoints
end module debye
