//
// Created by Christian Messe on 15.09.19.
//

#ifndef BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP
#define BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP

#include "cl_GM_EoS_AlphaFunction.hpp"

namespace belfem
{
    namespace gastables
    {
        class GasData;
    }

    namespace gasmodels
    {

//----------------------------------------------------------------------------
        class AlphaFunctionFactory
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

             AlphaFunctionFactory() = default;

            ~AlphaFunctionFactory() = default;

//----------------------------------------------------------------------------

            AlphaFunction *
            create_empty();

//----------------------------------------------------------------------------

            AlphaFunction *
            create_srk( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pr78( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_ccr_pr( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_ccr_mc_srk( const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pm_srk(  const gastables::GasData * aData );

//----------------------------------------------------------------------------

            AlphaFunction *
            create_pm_pr(  const gastables::GasData * aData );

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GM_ALPHAFUNCTIONFACTORY_HPP
