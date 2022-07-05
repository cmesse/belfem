//
// Created by Christian Messe on 02.12.19.
//

#ifndef BELFEM_CL_IWG_TATCAD_HPP
#define BELFEM_CL_IWG_TATCAD_HPP


#include "cl_IWG.hpp"
#include "en_SolverEnums.hpp"

#include "cl_IWG_TransientHeatConduction.hpp"

namespace belfem
{
    class Gas;

    namespace fem
    {
        /**
         * Transient Aerothermodynamic Convection and Diffusion
         */
        class IWG_TATCAD : public IWG_TransientHeatConduction
        {
            Gas * mGas = nullptr;

            // reference conditions
            real mT     = BELFEM_QUIET_NAN;
            real mP     = BELFEM_QUIET_NAN;
            real mU     = BELFEM_QUIET_NAN;
            real mRho   = BELFEM_QUIET_NAN;
            real mH     = BELFEM_QUIET_NAN;
            real mHt    = BELFEM_QUIET_NAN;

            real mTpow4 = BELFEM_QUIET_NAN;
            real mRhoU  = BELFEM_QUIET_NAN;
            real mRhoU2 = BELFEM_QUIET_NAN;

            // surface emissivity
            real mEpsilon =  0.83;

            Vector< index_t > mNodesOnHotgasSiteSets ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_TATCAD( const uint aNumberOfDimensions );

//------------------------------------------------------------------------------

            ~IWG_TATCAD() = default;

//------------------------------------------------------------------------------
// User Settings
//------------------------------------------------------------------------------

            /**
             * the reference gas for the convective boundary condition
             * @param aGas
             */
            void
            set_gas( Gas * aGas );

//------------------------------------------------------------------------------

            /**
             * freestream flow conditions
             *
             * @param aT  Temperature in K
             * @param aP  Pressure in Pa
             * @param aU  Velocity in m/s
             */
            void
            set_reference_conditions(
                    real & aT,
                    real & aP,
                    real & aU );

//------------------------------------------------------------------------------

            /**
             * sets the surface emmisivity to a constant factor
             */
            void
            set_surface_emmissivity( const real & aEpsilon );

//------------------------------------------------------------------------------

            /**
             * set hotgas side
             */
             void
             set_wetted_sidesets_hotcold(
                     const Vector< id_t > & aHotgasSidesets,
                     const Vector< id_t > & aColdgasSidesets );


//------------------------------------------------------------------------------
// kernel relevant functions
//------------------------------------------------------------------------------

            // this function must be called before the convection terms
            // are computed
            //fixme: make parallel
            void
            compute_heatloads();

//------------------------------------------------------------------------------

            // linearizes the enthalpy function for the wall to accellerate the
            // the computation fixme: make parallel
            void
            compute_surface_enthalpy();

//------------------------------------------------------------------------------
        };

    }
}

#endif //BELFEM_CL_IWG_TATCAD_HPP
