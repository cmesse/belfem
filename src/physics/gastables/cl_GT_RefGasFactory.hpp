//
// Created by Christian Messe on 30.08.19.
//

#ifndef BELFEM_CL_GT_REFGASFACTORY_HPP
#define BELFEM_CL_GT_REFGASFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_GT_RefGas.hpp"
#include "fn_GT_data_path.hpp"

namespace belfem
{
    class SpMatrix;

    namespace gastables
    {
        class InputThermo;
        class InputTransport;
        class InputData;
        class InputAlpha;

//------------------------------------------------------------------------------

        class RefGasFactory
        {
            InputThermo    * mCryoThermo;
            InputThermo    * mThermo;
            InputTransport * mCryoTransport;
            InputTransport * mTransport;
            InputData      * mData;
            InputAlpha     * mAlpha;

            Vector< real >   mTemperatures;
            SpMatrix         mHelpMatrix;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            RefGasFactory( const string aDataPath = gastables::data_path() );

//------------------------------------------------------------------------------

            ~RefGasFactory();

//------------------------------------------------------------------------------

            RefGas *
            create_refgas( const string & aLabel );

//------------------------------------------------------------------------------

            void
            create_temperature_steps( Vector< real > & aTemperatureSteps );

//------------------------------------------------------------------------------

            void
            create_helpmatrix( SpMatrix & aHelpMatrix );

//------------------------------------------------------------------------------

            bool
            interaction_viscosity_exists(
                    const string & aA,
                    const string & aB );

//------------------------------------------------------------------------------

            RefGas *
            create_interaction_viscosity( const string & aA, const string & aB );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_REFGASFACTORY_HPP
