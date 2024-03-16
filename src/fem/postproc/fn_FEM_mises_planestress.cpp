//
// Created by Christian Messe on 18.11.19.
//

#include "typedefs.hpp"
#include "commtools.hpp"

#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"

#include "fn_FEM_mises_planestress.hpp"

#include "cl_FEM_Group.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        void
        mises_planestress( DofManagerBase * aField )
        {
            if( comm_rank() == aField->parent()->master() )
            {
                // start timer
                Timer tTimer;

                // get pointer to mesh
                Mesh * tMesh = aField->mesh();

                // add a new field
                Vector< real > & tSigmaXX = tMesh->create_field( "sigma_11", EntityType::NODE );
                Vector< real > & tSigmaYY = tMesh->create_field( "sigma_22" , EntityType::NODE );
                Vector< real > & tSigmaXY = tMesh->create_field( "tau_12" , EntityType::NODE );
                Vector< real > & tMises   = tMesh->create_field( "mises" , EntityType::NODE );

                // allocate a C-Matrix
                Matrix< real > tC( 3, 3 );
                Vector< real > tEpsilon( 3 );
                Vector< real > tSigma( 3 );

                Vector< real > & tDUDX = tMesh->field_data( "grad_uxx");
                Vector< real > & tDUDY = tMesh->field_data( "grad_uxy");

                Vector< real > & tDVDX = tMesh->field_data( "grad_uyx");
                Vector< real > & tDVDY = tMesh->field_data( "grad_uyy");

                // loop over all blocks of field
                for ( id_t tBlockID : aField->iwg()->selected_blocks() )
                {
                    Block * tBlock = aField->block( tBlockID );

                    // get material of current block
                    const Material * tMat = tBlock->material();

                    BELFEM_ASSERT( tMat != nullptr, "no material was assigned to block %u",
                                  ( unsigned int ) tBlock->id() );

                    // unflag all nodes on mesh
                    tMesh->unflag_all_nodes();

                    // flag all nodes on block
                    tBlock->block()->flag_nodes();

                    // loop over all nodes on mesh
                    for( mesh::Node * tNode : tMesh->nodes() )
                    {
                        // test if node is flaged
                        if( tNode->is_flagged() )
                        {
                            index_t k = tNode->index();

                            // populate displacements
                            tEpsilon( 0 ) = tDUDX( k );
                            tEpsilon( 1 ) = tDVDY( k );
                            tEpsilon( 2 ) = tDUDY( k ) + tDVDX( k );

                            // populate elasticity
                            tMat->C_ps( tC );

                            // apply Hooke's law
                            tSigma = tC * tEpsilon;

                            // extract individual stresses
                            tSigmaXX( k ) = tSigma( 0 );
                            tSigmaYY( k ) = tSigma( 1 );
                            tSigmaXY( k ) = tSigma( 2 );

                            // compute van mises stress
                            tMises( k ) = std::sqrt( tSigma( 0 ) * tSigma( 0 )
                                                  +  tSigma( 1 ) * tSigma( 1 )
                                                  -  tSigma( 0 ) * tSigma( 1 )
                                                  +  3.0 * tSigma( 2 ) ) ; // in MPa
                        }
                    }

                    message( 4, "    ... time for computing Stresses             : %u ms\n",
                             ( unsigned int ) tTimer.stop() );
                }
            }
        }

//------------------------------------------------------------------------------
    }
}