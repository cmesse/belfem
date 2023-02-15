//
// Created by christian on 2/9/23.
//
#include "cl_IWG_Maxwell_Thermal2D.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_Thermal2D::IWG_Maxwell_Thermal2D() :
                IWG_Timestep( IwgType::TransientHeatConduction,
                     IwgMode::Iterative,
                     SymmetryMode::PositiveDefiniteSymmetric,
                     DofMode::AllBlocksEqual )
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

            uint l = 0 ;

            for( id_t b : aMaxwellMesh->ghost_block_ids() )
            {
                Cell< mesh::Element * > & tElements
                    = aMaxwellMesh->block( b )->elements() ;

                for( mesh::Element * tElement : tElements )
                {
                    // get index of corresponding node on thermal mesh
                    index_t tIndex = aThermalMesh->node( tElement->id() )->index() ;

                    mLayer( tIndex ) = l ;

                    // compute length of element
                    real tdX = tElement->node( 1 )->x() - tElement->node( 0 )->x() ;
                    real tdY = tElement->node( 1 )->y() - tElement->node( 0 )->y() ;

                    // write element length into container
                    mLength( tIndex ) = std::sqrt( tdX * tdX + tdY * tdY );
                }
                ++l ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->print_dofs( aElement );

            Matrix< real > & tM  = aJacobian ;
            Matrix< real > & tK  = mGroup->work_K() ;
            Vector< real > & tT  = mGroup->work_phi() ;
            Vector< real > & tT0 = mGroup->work_psi() ;
            Vector< real > & tB  = mGroup->work_tau() ;
            Vector< real > & tQ  = mGroup->work_chi() ;

            mesh::Node * tNodeA = aElement->element()->node( 0 ) ;
            mesh::Node * tNodeB = aElement->element()->node( 1 );

            // layer IDs
            uint tLayerA = mLayer( tNodeA->index() ) ;
            uint tLayerB = mLayer( tNodeB->index() );

            // grab materials
            Material * tMatA = mMaterials( tLayerA );
            Material * tMatB = mMaterials( tLayerB );

            this->collect_node_data( aElement, "T", tT );
            this->collect_node_data( aElement, "T0", tT0 );
            this->collect_node_data( aElement, "b", tB );
            this->collect_node_data( aElement, "ej", tQ );

            std::cout << tNodeA->is_flagged() << " " << tNodeA->owner() << std::endl ;

            // lumped mass part for node 0
            if( ( ! tNodeA->is_flagged() ) && tNodeA->owner() == gComm.rank() )
            {
                // volume
                real tV = mTapeThickness( tLayerA ) * mLength( tNodeA->index() );

                // contribution to lumped mass
                tM( 0, 0 ) = tV * tMatA->rho() * tMatA->c( tT( 0 ) ) ;

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
            if( ( ! tNodeB->is_flagged() ) && tNodeB->owner() == gComm.rank() )
            {
                // volume
                real tV = mTapeThickness( tLayerA ) * mLength( tNodeA->index() );

                // contribution to lumped mass
                tM( 1, 1 ) = tV * tMatB->rho() * tMatB->c( tT( 1 ) ) ;

                // contribution to heat flux
                aRHS( 1 ) = tV * tQ( 1 );

                tNodeB->flag() ;
            }
            else
            {
                tM( 1, 1 ) = 0.0 ;
                aRHS ( 0 ) = 0.0 ;
            }

            real tLa ;
            real tLb ;
            real tS ;
            real tKa = tMatA->lambda( tT( 0 ), tB( 0 ) );
            real tKb = tMatB->lambda( tT( 1 ), tB( 1 ) );

            if( tLayerA == tLayerB )
            {
                tLa = 0.5 * mLength( tNodeA->index() );
                tLb = 0.5 * mLength( tNodeB->index() );
                tS = mTapeThickness( tLayerA );

            }
            else
            {
                tLa = 0.5 * mTapeThickness( tLayerA );
                tLb = 0.5 * mTapeThickness( tLayerB );
                tS = mLength( tNodeA->index() );
            }

            real tValue = ( tKa * tKb * tS ) / ( tKa * tLb + tKb * tLa );
            tK( 0,0 )  =  tValue ;
            tK( 1, 0 ) = -tValue ;
            tK( 0, 1 ) = -tValue ;
            tK( 1, 1 ) =  tValue ;

            tM.print("M");
            tK.print("K");

            aRHS *= mDeltaTime ;
            aRHS += tM * tT0 ;
            aJacobian += mDeltaTime * tK ;

            aJacobian.print("J");
            aRHS.print("RHS");

            exit( 0 );
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

    }
}