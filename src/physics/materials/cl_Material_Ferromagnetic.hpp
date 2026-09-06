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

#ifndef BELFEM_CL_MATERIAL_FERROMAGNETIC_HPP
#define BELFEM_CL_MATERIAL_FERROMAGNETIC_HPP

#include "cl_Material_Metal.hpp"

namespace belfem
{
    namespace material
    {
        class Ferromagnetic : public Metal
        {
            real mBrillouinS = BELFEM_QUIET_NAN ;
            real mBrillouinP = BELFEM_QUIET_NAN ;
            real mBrillouinQ = BELFEM_QUIET_NAN ;
            real mBrillouinR = BELFEM_QUIET_NAN ;
            real mBrillouinP2 = BELFEM_QUIET_NAN ;
            real mBrillouinQ2 = BELFEM_QUIET_NAN ;
        public:

            Ferromagnetic( const string & aLabel,
                  const bool aBuildTables = true );

            ~Ferromagnetic() override;


            real
            compute_debye_from_rho( const real T, const real rho, const real theta_guess=BELFEM_QUIET_NAN ) override ;


            // to be moved into private
            real
            rho_mag( const real T ) const ;

            void
            set_angular_momentum( const real S );

        protected:

            real
            rho_custom( const real T ) const override;


            virtual real
            compute_mred( real T ) const ;


        private:
            void
            brillouin(
             real x,
             real & B,
             real & dBdx ) const ;
        };

    }
}

#endif //BELFEM_CL_MATERIAL_FERROMAGNETIC_HPP
