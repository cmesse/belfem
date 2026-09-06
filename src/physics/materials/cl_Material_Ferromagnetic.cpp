/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_Material_Ferromagnetic.hpp"

namespace belfem
{
    namespace material
    {
        Ferromagnetic::Ferromagnetic(
            const string & aLabel,
            const bool aBuildTables ) :
            Metal( aLabel, MaterialType::PureMetal, aBuildTables )
        {
            this->set_have( MaterialProperty::mu, false );
        }

        Ferromagnetic::~Ferromagnetic()
        {

        }

        void
        Ferromagnetic::set_angular_momentum( const real S )
        {
            mBrillouinS  = S ;
            mBrillouinQ  = 0.5/S ;
            mBrillouinQ2 = mBrillouinQ*mBrillouinQ;
            mBrillouinP  = mBrillouinQ*(S+S+1.);
            mBrillouinP2 =mBrillouinP*mBrillouinP;
            mBrillouinR = (S+S+S)/(S+1.);
        }

        real
        Ferromagnetic::compute_mred( real T ) const
        {
            real Tc = this->constant_property( MaterialProperty::Tcurie );
            if (T >= Tc) return 0.0;

            real t = T / Tc;

            // Initial guess
            real m = 0.5 ;

            uint tCount = 0;
            real f = 1 ;
            while ( std::abs(f) > BELFEM_EPSILON )
            {
                real y = mBrillouinR * m / t;     // argument to Brillouin function

                real B, dBdy;
                this->brillouin(y, B, dBdy);

                f = B - m;

                real dfdm = dBdy * mBrillouinR / t - 1.0;
                m -= f / dfdm;                    // ← Newton update on m

                BELFEM_ERROR( tCount++ < 100, "compute_mred did not converge for T = %g", T);
            }



            return m;
        }

        real
        Ferromagnetic::rho_custom( const real T ) const
        {
            // resistive resistivity
            real rho_0 = this->constant_property( MaterialProperty::rho_0 );

            if ( T < BELFEM_EPSILON ) return rho_0 ;

            return rho_0 + this->rho_i_custom( T ) + this->rho_mag( T );
        }

        inline real
        Ferromagnetic::compute_debye_from_rho( const real T, const real rho, const real theta_guess )
        {
            return Metal::compute_debye_from_rho( T, rho-this->rho_mag( T ), theta_guess);
        }

        real
        Ferromagnetic::rho_mag( const real T ) const
        {
            // magnetization
            real m = this->compute_mred( T );

            // magnetic weight
            real w = 1. - m * m ;

            // blend of low and high temperatue value
            return  std::max( m * m * this->constant_property( MaterialProperty::A_electron_magnon) * T * T
                  + w * w * this->constant_property( MaterialProperty::A_spin_disorder ), 0.0 );
        }

        void
        Ferromagnetic::brillouin(
           real x,
           real & B,
           real & dBdx ) const
        {
            // minimizing the usage of expensive exp function
            real fP = std::exp( mBrillouinP * x );
            real gP = 1./fP;

            real sinhP = 0.5*(fP-gP);
            real cothP = (fP+gP)/(fP-gP);

            real fQ = std::exp( mBrillouinQ * x );
            real gQ = 1./fQ;
            real sinhQ = 0.5*(fQ-gQ);
            real cothQ = (fQ+gQ)/(fQ-gQ);

            B = mBrillouinP * cothP - mBrillouinQ * cothQ ;
            dBdx = mBrillouinQ2/(sinhQ*sinhQ) - mBrillouinP2/(sinhP*sinhP);
        }

    }
}