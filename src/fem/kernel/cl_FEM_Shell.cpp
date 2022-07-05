//
// Created by christian on 3/29/22.
//


#include "cl_FEM_Group.hpp"
#include "cl_FEM_Shell.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Shell::Shell(
                DofManager * aParent,
                const id_t aID,
                Cell< mesh::Facet * > & aFacets ) :
                SideSet(  aParent, aID, aFacets, GroupType::SHELL )
        {
            this->set_domain_type( DomainType::ThinShell );

            // get the element type
            switch( mesh::interpolation_order_numeric( mElementType )  )
            {
                case( 1 ) :
                {
                    mThinShellIntegrationData = new IntegrationData( ElementType::LINE2 );
                    break ;
                }
                case( 2 ) :
                {
                    mThinShellIntegrationData = new IntegrationData( ElementType::LINE3 );
                    break ;
                }
                case( 3 ) :
                {
                    mThinShellIntegrationData = new IntegrationData( ElementType::LINE4 );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid element order for sideset %lu",
                                  ( long unsigned int ) this->id() );
                }
            }


            mThinShellIntegrationData->populate( aParent->sideset_integration_order(),
                                        aParent->integration_scheme() );

        }

//------------------------------------------------------------------------------

        Shell::~Shell()
        {
            delete mThinShellIntegrationData ;
        }

//------------------------------------------------------------------------------
    }
}