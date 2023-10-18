//
// Created by christian on 2/9/23.
//
#include "commtools.hpp"
#include "cl_IWG_Maxwell_Thermal2D.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_Thermal2D::IWG_Maxwell_Thermal2D() :
                IWG_TimestepOld( IwgType::TransientHeatConduction,
                     IwgMode::Iterative,
                     SymmetryMode::PositiveDefiniteSymmetric,
                     DofMode::AllBlocksEqual ),
                     mMyRank( comm_rank() )
        {
            mNumberOfDofsPerNode = 1 ;
            mNumberOfDerivativeDimensions = 1 ;
            mNumberOfRhsCols = 1 ;

            // set the names for the fields
            mDofFields = { "T" };
            mFluxFields = { "ej" };
            mOtherFields = { "b", "T0" };

            this->initialize() ;
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_Thermal2D::~IWG_Maxwell_Thermal2D()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::set_material_table(
                Cell< Material * > & aTapeMaterials )
        {
            mMaterials.set_size( aTapeMaterials.size(), nullptr );

            for ( uint l=0; l<aTapeMaterials.size(); ++l )
            {
                mMaterials( l ) = aTapeMaterials( l );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::compute_geometry_data(
                Mesh * aMaxwellMesh,
                Mesh * aThermalMesh,
                const Vector< real > & aTapeThicknesses )
        {
            mTapeThickness = aTapeThicknesses ;

            mLayer.set_size( aThermalMesh->number_of_nodes(), BELFEM_UINT_MAX );
            mLength.set_size( aThermalMesh->number_of_nodes(), BELFEM_QUIET_NAN );

            if( mMyRank == 0 )
            {
                uint l = 0;

                for ( id_t b: aMaxwellMesh->ghost_block_ids())
                {
                    Cell< mesh::Element * > & tElements
                            = aMaxwellMesh->block( b )->elements();

                    for ( mesh::Element * tElement: tElements )
                    {
                        // get index of corresponding node on thermal mesh
                        index_t tIndex = aThermalMesh->node( tElement->id())->index();

                        mLayer( tIndex ) = l;

                        // compute length of element
                        real tdX = tElement->node( 1 )->x() - tElement->node( 0 )->x();
                        real tdY = tElement->node( 1 )->y() - tElement->node( 0 )->y();

                        // write element length into container
                        mLength( tIndex ) = std::sqrt( tdX * tdX + tdY * tdY );
                    }
                    ++l;
                }



                // create comm list
                proc_t tNumProcs = comm_size() ;

                // cancel if this is just one proc
                if( tNumProcs < 2 ) return;

                Vector< proc_t > tCommList( tNumProcs-1 );
                index_t tCount = 0 ;
                for( proc_t p=0; p<tNumProcs; ++p )
                {
                    if( p != mMyRank )
                    {
                        tCommList( tCount++ ) = p ;
                    }
                }

                comm_barrier() ;

                // request list
                Cell< Vector< id_t > > tAllRequest( tNumProcs-1, {} ) ;

                // length list
                Cell< Vector< real > > tAllLength( tNumProcs-1, {} ) ;

                // layer list
                Cell< Vector< uint > > tAllLayer( tNumProcs-1, {} ) ;

                receive( tCommList, tAllRequest );

                for( proc_t p=0; p<tNumProcs-1; ++p )
                {
                    Vector< id_t > & tRequest = tAllRequest( p );
                    Vector< real > & tLength  = tAllLength( p );
                    Vector< uint > & tLayer   = tAllLayer( p );

                    // reset counter
                    tCount = 0 ;

                    tLength.set_size( tRequest.length() );
                    tLayer.set_size( tRequest.length() );

                    for( id_t tID : tRequest )
                    {
                        index_t tIndex = aThermalMesh->node( tID )->index();
                        tLength( tCount ) = mLength( tIndex );
                        tLayer( tCount++ ) = mLayer( tIndex );
                    }
                }

                send( tCommList, tAllLength );
                send( tCommList, tAllLayer );
            }
            else
            {
                Cell< mesh::Node * > & tNodes = aThermalMesh->nodes() ;

                // create request list
                Vector< id_t > tRequest( aThermalMesh->number_of_nodes() );

                index_t tCount = 0 ;

                for( mesh::Node * tNode : tNodes )
                {
                    tRequest( tCount++ ) = tNode->id() ;
                }

                comm_barrier() ;

                send( 0, tRequest );
                receive( 0, mLength );
                receive( 0, mLayer );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // mass matrix
            Matrix< real > & tM  = aJacobian ;

            // stiffness matrix
            Matrix< real > & tK  = mGroup->work_K() ;

            // degrees of freedom
            Vector< real > & tT  = mGroup->work_phi() ;
            this->collect_node_data( aElement, "T", tT );

            // degrees of freedom at last timestep
            Vector< real > & tT0 = mGroup->work_psi() ;
            this->collect_node_data( aElement, "T0", tT0 );

            // magnetic flux density
            Vector< real > & tB  = mGroup->work_tau() ;
            this->collect_node_data( aElement, "b", tB );

            // heat load
            Vector< real > & tQ  = mGroup->work_chi() ;
            this->collect_node_data( aElement, "ej", tQ );

            // grab nodes of element
            mesh::Node * tNodeA = aElement->element()->node( 0 ) ;
            mesh::Node * tNodeB = aElement->element()->node( 1 );

            // layer IDs
            uint tLayerA = mLayer( tNodeA->index() ) ;
            uint tLayerB = mLayer( tNodeB->index() );

            // grab materials
            Material * tMatA = mMaterials( tLayerA );
            Material * tMatB = mMaterials( tLayerB );

           // std::cout << "dof : " << aElement->dof( 0 )->value() << " " << aElement->dof( 1 )->value() << std::endl ;

            // lumped mass part for node 0
            if( ( ! tNodeA->is_flagged() ) && tNodeA->owner() == mMyRank )
            {
                // volume
                real tV = mTapeThickness( tLayerA ) * mLength( tNodeA->index() );


                // contribution to lumped mass
                tM( 0, 0 ) = tV * tMatA->rho() * tMatA->c( tT( 0 ) ) ;

                //std::cout << "n : " << tNodeA->id() <<
                //          " V: " << tV <<
                //          " T: " << tT( 0 ) <<
                //          " rho: " <<  tMatA->rho() <<
                //          " cp:"  << tMatA->c( tT( 0 ) ) << std::endl ;

                // contribution to heat flux
                aRHS( 0 ) = tV * tQ( 0 );

                tNodeA->flag() ;
            }
            else
            {
                tM( 0, 0 ) = 0.0 ;
                aRHS( 0 ) = 0.0 ;
            }

            tM( 1, 0 ) = 0.0 ;
            tM( 0, 1 ) = 0.0 ;

            // lumped mass part for node 1
            if( ( ! tNodeB->is_flagged() ) && tNodeB->owner() == mMyRank )
            {
                // volume
                real tV = mTapeThickness( tLayerB ) * mLength( tNodeB->index() );

                // contribution to lumped mass
                tM( 1, 1 ) = tV * tMatB->rho() * tMatB->c( tT( 1 ) ) ;

                // contribution to heat flux
                aRHS( 1 ) = tV * tQ( 1 );

                tNodeB->flag() ;


                //std::cout << "n : " << tNodeB->id() <<
                //          " V: " << tV <<
                //          " T: " << tT( 1 ) <<
                //          " rho: " <<  tMatB->rho() <<
                //          " cp:"  << tMatB->c( tT( 1 ) ) << std::endl ;

            }
            else
            {
                tM( 1, 1 ) = 0.0 ;
                aRHS ( 1 ) = 0.0 ;
            }

            real tLa ;
            real tLb ;
            real tS ;
            //real tKa = std::max( tMatA->lambda( tT( 0 ), tB( 0 ) ), BELFEM_EPSILON );
            //real tKb = std::max( tMatB->lambda( tT( 1 ), tB( 1 )  ), BELFEM_EPSILON );
            real tKa = tMatA->lambda( tT( 0 ), tB( 0 ) ) ;
            real tKb = tMatB->lambda( tT( 1 ), tB( 1 )  ) ;

            //std::cout << "el : " << aElement->id() <<
            //    " T: " << tT( 0 ) << " " << tT( 1 ) <<
            //    " B: " <<  tB( 0 )  << " " << tB( 1 ) <<
            //    " k:"  << tKa << " " << tKb << " " << std::endl ;

            if( tLayerA == tLayerB )
            {
                // horizontal connector
                tLa = 0.5 * mLength( tNodeA->index() );
                tLb = 0.5 * mLength( tNodeB->index() );
                tS = mTapeThickness( tLayerA );
            }
            else
            {
                // vertical connector
                tLa = 0.5 * mTapeThickness( tLayerA );
                tLb = 0.5 * mTapeThickness( tLayerB );
                tS = mLength( tNodeA->index() );
            }

            real tValue = ( tKa * tKb * tS ) / ( tKa * tLb + tKb * tLa );
            tK( 0,0 )  =  tValue ;
            tK( 1, 0 ) = -tValue ;
            tK( 0, 1 ) = -tValue ;
            tK( 1, 1 ) =  tValue ;


            aRHS *= mDeltaTime ;
            aRHS += tM * tT0 ;
            aJacobian += mDeltaTime * tK ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::allocate_work_matrices( Group * aGroup )
        {
            aGroup->work_K().set_size( 2, 2 );
            aGroup->work_phi().set_size( 2 );
            aGroup->work_psi().set_size( 2 );
            aGroup->work_tau().set_size( 2 );
            aGroup->work_chi().set_size( 2 );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::shift_fields()
        {
            mMesh->field_data( "T0" ) = mMesh->field_data( "T");
        }

//------------------------------------------------------------------------------

    }
}