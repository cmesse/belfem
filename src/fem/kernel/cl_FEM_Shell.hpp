//
// Created by christian on 3/29/22.
//

#ifndef BELFEM_CL_FEM_SHELL_HPP
#define BELFEM_CL_FEM_SHELL_HPP

#include "cl_FEM_SideSet.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class DofManager ;

        class Shell :  public SideSet
        {
            // these integration data are needed for the integration
            // along thickness
            IntegrationData * mThinShellIntegrationData ;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Shell( DofManager * aParent,
                 const id_t aID,  Cell< mesh::Facet * > & aFacets );

//------------------------------------------------------------------------------

            ~Shell();

//------------------------------------------------------------------------------

            const IntegrationData *
            thinshell_integration() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Shell::thinshell_integration() const
        {
            return mThinShellIntegrationData ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_SHELL_HPP
