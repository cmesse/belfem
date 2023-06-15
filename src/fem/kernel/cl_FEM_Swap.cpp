//
// Created by christian on 6/9/23.
//
#include "cl_IWG.hpp"
#include "cl_FEM_Swap.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Swap::Swap( Group * aGroup ) :
            mGroup( aGroup ),
            mMesh( aGroup->parent()->mesh() )
        {

        }

//------------------------------------------------------------------------------

        void
        Swap::allocate()
        {
            // get the number of nodes per element
            uint tNumNodes = mesh::number_of_nodes( mGroup->element_type() );

            uint tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->element_type() );

            // get the number of dimensions
            uint tNumDimensions = mMesh->number_of_dimensions() ;

            // allocate vector for coordinates
            mX.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

            if ( mHaveThermal )
            {
                mT.set_size( tNumNodes, BELFEM_QUIET_NAN  );
            }

            if( mHaveMechanical )
            {
                mux.set_size( tNumNodes, BELFEM_QUIET_NAN  );
                muy.set_size( tNumNodes, BELFEM_QUIET_NAN  );
                if ( tNumDimensions == 3 )
                {
                    muz.set_size( tNumNodes, BELFEM_QUIET_NAN  );
                }
            }

            if( mHaveMaxwellPhi )
            {
                mphi.set_size( tNumNodes, BELFEM_QUIET_NAN );
            }

            if( mHaveMaxwellA )
            {
                ma.set_size( tNumNodes, BELFEM_QUIET_NAN  );
            }

            if( mHaveMaxwellH )
            {
                mh.set_size( tNumNedelecDofs, BELFEM_QUIET_NAN );
            }

            if( mHaveLagrange )
            {
                mlambda.set_size( tNumNodes, BELFEM_QUIET_NAN );
            }

            // check if we are allocating master and slave elements
            if( mGroup->master_type() != ElementType::UNDEFINED )
            {
                tNumNodes = mesh::number_of_nodes( mGroup->master_type() );
                tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->master_type() );

                mXm.set_size( tNumNodes, tNumDimensions+1 );

                if( mHaveMaxwellPhi )
                {
                    mphim.set_size( tNumNodes, BELFEM_QUIET_NAN  );
                }
                if( mHaveMaxwellH )
                {
                    mhm.set_size( tNumNedelecDofs, BELFEM_QUIET_NAN );
                }
                if( mHaveMaxwellA )
                {
                    mam.set_size( tNumDimensions==2 ? tNumNodes : tNumNedelecDofs,
                                  BELFEM_QUIET_NAN);
                }
                if( mHaveLagrange )
                {
                    mlambdam.set_size(
                            mesh::number_of_nedelec_dofs( mGroup->element_type() ),
                            BELFEM_QUIET_NAN );
                }
                mlambda.set_size( tNumNodes, BELFEM_QUIET_NAN );
            }

            if( mGroup->slave_type() != ElementType::UNDEFINED )
            {
                tNumNodes = mesh::number_of_nodes( mGroup->master_type() );
                tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->slave_type() );
                mXs.set_size( tNumNodes, tNumDimensions+1 );

                if( mHaveMaxwellPhi )
                {
                    mphis.set_size( tNumNodes, BELFEM_QUIET_NAN  );
                }

                if( mHaveMaxwellH )
                {
                    mhs.set_size( tNumNedelecDofs, BELFEM_QUIET_NAN );
                }
                if( mHaveMaxwellA )
                {
                    mas.set_size( tNumDimensions==2 ? tNumNodes : tNumNedelecDofs,
                                  BELFEM_QUIET_NAN );
                }
                if( mHaveLagrange )
                {
                    mlambdas.set_size(
                            mesh::number_of_nedelec_dofs( mGroup->element_type() ),
                            BELFEM_QUIET_NAN );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Swap::set_switches()
        {
            // make sure that field is set
            if( mGroup->parent()->iwg() != nullptr )
            {
                // grab the fields
                const Cell< string > & tFields =  mGroup->parent()->iwg()->all_fields() ;

                // loop over all fields
                // ( we need to use the wrapped std::vector here )
                for ( const string & tField : tFields.vector_data() )
                {
                    if( tField == "T" ) // temperature field
                    {
                        mHaveThermal = true ;
                    }
                    else if ( tField == "ux" ) // mechanical field
                    {
                        mHaveMechanical = true ;
                    }
                    else if ( tField == "phi" ) // magnetic scalar potential
                    {
                        mHaveMaxwellPhi = true ;
                    }
                    else if ( tField == "az" or tField == "edge_a") // magnetic vector potential
                    {
                        mHaveMaxwellA = true ;
                    }
                    else if ( tField == "edge_h" ) // h-field
                    {
                        mHaveMaxwellH = true ;
                    }
                    else if ( tField == "jz" ) // b-field
                    {
                        mHaveMaxwellJ = true ;
                    }
                    else if ( tField == "jz" ) // b-field
                    {
                        mHaveMaxwellJ = true ;
                    }
                    else if ( tField.substr(0,6) == "lambda" ) // lagrange
                    {
                        mHaveLagrange = true ;
                    }
                }
            }


        }

//------------------------------------------------------------------------------
    }
}