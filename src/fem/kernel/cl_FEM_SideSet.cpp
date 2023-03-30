//
// Created by Christian Messe on 03.11.19.
//
#include "cl_Mesh.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_IWG.hpp"
#include "commtools.hpp"
#include "meshtools.hpp"

#include "cl_FEM_Element.hpp"

#include "fn_IF_initialize_integration_points.hpp"
#include "fn_IF_initialize_shape_function.hpp"
#include "fn_intpoints_auto_integration_order.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        SideSet::SideSet(
                Field * aParent,
                const id_t     aID,
                Cell< mesh::Facet * > & aFacets ) :
                Group( aParent, GroupType::SIDESET,
                       aFacets.size() > 0
                       ? ( aParent->enforce_linear_interpolation() ?
                           mesh::linear_element_type( aFacets( 0 )->element()->type()  ) :
                           aFacets( 0 )->element()->type() )
                       : ( ElementType::UNDEFINED ),
                       aID, aFacets.size() )
        {
            // allocate set-wise boundary conditions
            mBcValues.set_size( aParent->number_of_dofs_per_node(), 0.0 );
            mBcTypes.set_size( aParent->number_of_dofs_per_node(), BoundaryConditionImposing::Free );

            if( mElementType != ElementType::UNDEFINED )
            {
                // set the number of nodes per element
                mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType );

                if( mParent->iwg() == NULL )
                {
                    mIntegrationData = new IntegrationData( mElementType,
                                                            InterpolationType::LAGRANGE );
                }
                else
                {
                    mIntegrationData = new IntegrationData( mElementType,
                                                            mParent->iwg()->interpolation_type());
                }

                mIntegrationData->populate( aParent->sideset_integration_order(),
                                            aParent->integration_scheme() );

                this->assume_isogeometry();
                this->initialize_elements( aFacets );
                this->collect_nodes( aFacets );
                this->create_element_map();
            }
        }

//------------------------------------------------------------------------------

        SideSet::SideSet( DofManager * aParent,
                          const id_t aID,
                          Cell< mesh::Facet * > & aFacets,
                          const GroupType aGroupType  ) :
                Group( aParent, aGroupType,
                       aFacets.size() > 0
                       ? ( aParent->enforce_linear_interpolation() ?
                           mesh::linear_element_type( aFacets( 0 )->element()->type()  ) :
                           aFacets( 0 )->element()->type() )
                       : ( ElementType::UNDEFINED ),
                       aID, aFacets.size() ),
                       mMasterType( aFacets.size() == 0 || aGroupType == GroupType::CUT ? ElementType::EMPTY :
                         aFacets( 0 )->master()->type() ),
                       mSlaveType( aFacets.size() == 0 || aGroupType == GroupType::CUT ? ElementType::EMPTY :
                       ( aFacets( 0 )->has_slave() ? aFacets( 0 )->slave()->type() : ElementType::EMPTY ) )
        {
            // get the number of dofs per node
            uint tNumDofTypes =  aParent->iwg()->dof_entity_types().length() ;

            // allocate set-wise boundary conditions
            mBcValues.set_size( tNumDofTypes, 0.0 );
            mBcTypes.set_size( tNumDofTypes, BoundaryConditionImposing::Free );

            // catch special case
            if( aGroupType == GroupType::CUT )
            {
                mNumberOfNodesPerElement = 2 ;

                this->initialize_elements( aFacets );
                this->collect_nodes( aFacets );
                this->create_element_map();
            }
            else if( mElementType != ElementType::UNDEFINED )
            {
                // set the number of nodes per element
                mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType );


                this->assume_isogeometry();
                this->initialize_elements( aFacets );
                this->collect_nodes( aFacets );
                this->create_element_map();

                uint tOrder = 0 ;
                if( aParent != nullptr )
                {
                    tOrder = aParent->sideset_integration_order() ;
                }
                this->initialize_lookup_tables( tOrder );
                this->initialize_work_matrices();
            }
        }

//------------------------------------------------------------------------------

        SideSet::SideSet( const ElementType aElementType,
                 const ElementType aMasterType,
                 const ElementType aSlaveType,
                 const GroupType aGroupType  ) :
                Group( nullptr, aGroupType, aElementType, 0, 0 ),
                mMasterType( aMasterType ),
                mSlaveType( aSlaveType )
        {
            // set the number of nodes per element
            mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType );


            // initialize lookup tables for integration of master and slave
            this->initialize_lookup_tables( 0 ) ;


            this->initialize_work_matrices();
        }

//------------------------------------------------------------------------------

        SideSet::~SideSet()
        {
            if( mElementType != ElementType::UNDEFINED )
            {
                for ( IntegrationData * tData : mMasterIntegration )
                {
                    delete tData ;
                }

                for ( IntegrationData * tData : mSlaveIntegration )
                {
                    delete tData ;
                }

                this->delete_pointers();
           }

        }

//------------------------------------------------------------------------------

        void
        SideSet::initialize_elements( Cell< mesh::Facet * > & aFacets  )
        {
            // get size of container
            index_t tNumberOfFacets = aFacets.size();

            // allocate memory
            mElements.set_size( tNumberOfFacets, nullptr );

            // legacy support
            switch ( mParent->type() )
            {
                case( DofManagerType::OLD ) :
                {
                    Field * tParent = reinterpret_cast< Field * >( mParent );

                    // create elements
                    for( index_t e=0; e<tNumberOfFacets; ++e )
                    {
                        mElements( e ) = new Element( this, tParent,
                                                      aFacets( e )->element() );
                    }
                    break;
                }
                case( DofManagerType::NEW ) :
                {

                    DofManager * tParent = reinterpret_cast< DofManager * >( mParent );

                    switch( mType )
                    {
                        case( GroupType::CUT ) :
                        {
                            index_t tCount = 0 ;

                            // loop over all facets
                            for ( mesh::Facet * tFacet: aFacets )
                            {
                                mElements( tCount++ ) = new Element(
                                        this,
                                        tParent,
                                        tFacet,
                                        SideSetDofLinkMode::Cut,
                                        tFacet->element()->block_id(), // block id contains plus block, needed for dof types
                                        0 );
                            }
                            break ;
                        }
                        case( GroupType::SHELL ) :
                        {
                            // counter for new elements
                            index_t tCount = 0 ;

                            // get number of layers
                            index_t tNumberOfLayers = mParent->mesh()->number_of_thin_shell_layers() ;


                            // allocate facet container
                            Cell< mesh::Facet * > tFacets( tNumberOfLayers, nullptr );

                            for ( mesh::Facet * tFacet: aFacets )
                            {
                                // collect facets
                                for ( index_t l = 0; l < tNumberOfLayers; ++l )
                                {
                                    tFacets( l ) = mParent->mesh()->ghost_facet( tFacet->id(), l );
                                }

                                // create the new element
                                mElements( tCount++ ) = new Element(
                                        this,
                                        tParent,
                                        tFacet,
                                        tFacets,
                                        tFacet->master()->block_id(), // block id contains plus block, needed for dof types
                                        tFacet->slave()->block_id() );
                            }

                            break ;
                        }
                        default :
                        {
                            // flag relevant elements on master block
                            tParent->mesh()->unflag_all_elements() ;

                            for ( mesh::Facet * tFacet: aFacets )
                            {
                                BELFEM_ASSERT( tFacet->master() != nullptr,
                                               "Master of facet %lu must not be null",
                                               ( long unsigned int ) tFacet->id());

                                tFacet->master()->flag();

                                if ( tFacet->slave() != nullptr )
                                {
                                    tFacet->slave()->flag();
                                }
                            }

                            // initialize counter
                            index_t tCount = 0;

                            // the block map tells which element sits on which block
                            Map< id_t, id_t > tBlockMap;
                            for ( mesh::Block * tBlock: tParent->mesh()->blocks())
                            {
                                for ( mesh::Element * tElement: tBlock->elements())
                                {
                                    if ( tElement->is_flagged())
                                    {
                                        tBlockMap[ tElement->id() ] = tBlock->id();
                                    }
                                }
                            }

                            // get the linking mode
                            SideSetDofLinkMode tMode = mParent->iwg()->sideset_dof_link_mode();

                            // loop over all facets
                            for ( mesh::Facet * tFacet: aFacets )
                            {
                                mElements( tCount++ ) = new Element(
                                        this,
                                        tParent,
                                        tFacet,
                                        tMode,
                                        tBlockMap( tFacet->master()->id()),
                                        tFacet->slave() == nullptr ?
                                        0 :
                                        tBlockMap( tFacet->slave()->id()));
                            }
                            break ;
                        }
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid DOF Manager type");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::collect_nodes( Cell< mesh::Facet * > & aFacets )
        {
            // reset node container
            mParent->mesh()->unflag_all_nodes();

            if ( mParent->enforce_linear_interpolation() )
            {
                // flag all corner nodes that belong to this set
                for( mesh::Facet * tFacet : aFacets )
                {
                    tFacet->flag_corner_nodes();
                }
            }
            else
            {
                // flag all nodes that belong to this set
                for ( mesh::Facet * tFacet : aFacets )
                {
                    tFacet->flag_nodes();
                }
            }

            // count flagged nodes
            index_t tCount = 0;
            for( mesh::Node * tNode : mParent->mesh()->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            // allocate container
            mNodes.set_size( tCount, nullptr );

            // reset counter
            tCount = 0;

            // collect nodes
            for( mesh::Node * tNode : mParent->mesh()->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    // write node into container
                    mNodes( tCount ++ ) = tNode;

                    // tidy up
                    tNode->unflag();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::impose_dirichlet( const real aValue, const uint aDofType )
        {
            // set value
            mBcValues( aDofType ) = aValue ;

            // set type
            mBcTypes( aDofType ) = BoundaryConditionImposing::Dirichlet ;

            // loop over all nodes of this mesh
            for ( mesh::Node * tNode : mNodes )
            {
                // grab dof and fix value
                mParent->dof( mParent->calculate_dof_id( tNode, aDofType ) )->fix( aValue );
            }

        }

//------------------------------------------------------------------------------

        void
        SideSet::impose_neumann( const real aValue, const uint aDofType )
        {
            // set value
            mBcValues( aDofType ) = aValue ;

            mBcTypes( aDofType ) = BoundaryConditionImposing::Neumann ;
        }

//------------------------------------------------------------------------------

        void
        SideSet::impose_alpha( const real aAlpha, const real aTinf )
        {
            // remember value of alpha
            mBcValues( 0 ) = aAlpha ;

            mBcTypes( 0 ) = BoundaryConditionImposing::Alpha ;
            mTinf = std::abs( aTinf ) >= 0 ? aTinf : BELFEM_TREF ;

            // create the fields if they don't exist already
            if( ! mParent->mesh()->field_exists( "alpha") )
            {
                Vector< real > & tAlpha = mParent->mesh()->create_field( "alpha" ) ;
                tAlpha.fill( 0.0 );
            }

            if( ! mParent->mesh()->field_exists( "Tinf") )
            {
                Vector< real > & tTinf = mParent->mesh()->create_field( "Tinf" ) ;
                tTinf.fill( mTinf );
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::free()
        {
            uint n = mParent->number_of_dofs_per_node();

            for( uint k=0; k<n; ++k )
            {
                mBcTypes( k ) = BoundaryConditionImposing::Free ;
            }

            // loop over all nodes of this set
            for ( mesh::Node * tNode : mNodes )
            {
                for( uint k=0; k<n; ++k )
                {
                    // grab dof and free value
                    mParent->dof( mParent->calculate_dof_id( tNode, k ) )->free();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::set_boundary_conditions()
        {
            uint tNumDofsPerNode = mBcTypes.size();

            for( uint k=0; k<tNumDofsPerNode; ++k )
            {
                switch( mBcTypes( k ) )
                {
                    case( BoundaryConditionImposing::Free ) :
                    {
                        /* do nothing */
                        break ;
                    }
                    case( BoundaryConditionImposing::Dirichlet ) :
                    {
                        /* do nothing */
                        break ;
                    }
                    case( BoundaryConditionImposing::Neumann ) :
                    {
                        // grab surface field
                        Vector< real > & tField
                                = mParent->field_data( mParent->iwg()->field(
                                        tNumDofsPerNode + k ) );

                        // loop over all nodes of this mesh
                        for ( mesh::Node * tNode : mNodes )
                        {
                            tField( tNode->index() ) = mBcValues( k );
                        }
                        break ;
                    }
                    case( BoundaryConditionImposing::Alpha ) :
                    {
                        // grab surface field
                        Vector< real > & tAlpha  = mParent->field_data( "alpha" );
                        Vector< real > & tTinf   = mParent->field_data( "Tinf" );

                        // loop over all nodes of this mesh
                        for ( mesh::Node * tNode : mNodes )
                        {
                            tAlpha( tNode->index() ) = mBcValues( 0 );
                            tTinf( tNode->index()  ) = mTinf ;
                        }
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "don't know what to do here");
                        break ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::initialize_lookup_tables( const uint aIntegrationOrder )

        {
            // determine integration order
            uint tIntegrationOrder = aIntegrationOrder ;
            InterpolationType tType = mParent == nullptr ? InterpolationType::LAGRANGE : mParent->iwg()->interpolation_type();
            IntegrationScheme tScheme = mParent == nullptr ? IntegrationScheme::GAUSSCLASSIC : mParent->integration_scheme() ;


            bool tHaveMaster =  mMasterType != ElementType::UNDEFINED && mMasterType != ElementType::EMPTY ;
            bool tHaveSlave =   mSlaveType != ElementType::UNDEFINED && mSlaveType != ElementType::EMPTY ;

            // check if we use auto setting
            if( tIntegrationOrder == 0 )
            {
                // auto define integration order
                tIntegrationOrder = auto_integration_order( mElementType );

                if( tHaveMaster )
                {
                    tIntegrationOrder = std::max(
                            tIntegrationOrder,
                            auto_integration_order( mMasterType ) );
                }

                if( tHaveSlave )
                {
                    tIntegrationOrder = std::max(
                            tIntegrationOrder,
                            auto_integration_order( mSlaveType ) );
                }
            }

            // note that the integration data is not populated!
            mIntegrationData = new IntegrationData( mElementType );

            // populate main integration data
            mIntegrationData->populate( tIntegrationOrder, tScheme );

            if( tHaveMaster )
            {
                uint tNumFacets = mesh::number_of_facets( mMasterType );
                mMasterIntegration.set_size( tNumFacets, nullptr );
                for( uint f=0; f<tNumFacets; ++f )
                {
                    mMasterIntegration( f ) = new IntegrationData( mMasterType,
                                                                   tType );
                    mMasterIntegration( f )->populate_for_master(
                            f , tIntegrationOrder, tScheme );
                }
            }
            if( tHaveSlave )
            {
                // count number of permutations
                uint tNumFacets = mesh::number_of_facets( mSlaveType );

                uint tNumPermutations = 0 ;
                for( uint f=0; f<tNumFacets; ++f )
                {
                    tNumPermutations += mesh::number_of_orientations( mSlaveType, f );
                }

                mSlaveIntegration.set_size( tNumPermutations, nullptr );

                uint tCount = 0 ;

                for( uint f=0; f<tNumFacets; ++f )
                {
                    for( uint o=0; o<mesh::number_of_orientations( mSlaveType, f );  ++o )
                    {
                        mSlaveIntegration( tCount ) = new IntegrationData( mSlaveType,  tType );
                        mSlaveIntegration( tCount++ )->populate_for_slave( f, o, tIntegrationOrder, tScheme );
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        SideSet::initialize_work_matrices()
        {
            uint tNumDim = 0 ;

            if(    mMasterType != ElementType::UNDEFINED
                && mMasterType != ElementType::EMPTY )
            {
                tNumDim = mesh::dimension( mMasterType );

                this->work_Xm().set_size( mesh::number_of_nodes( mMasterType ), tNumDim );
                this->work_normal().set_size( tNumDim );

            }

            if(       mSlaveType != ElementType::UNDEFINED
                   && mSlaveType != ElementType::EMPTY )
            {
                if( tNumDim == 0 )
                {
                    tNumDim = mesh::dimension( mSlaveType );
                }
                this->work_Xs().set_size( mesh::number_of_nodes( mSlaveType ), tNumDim );

            }

            if(          mElementType != ElementType::UNDEFINED
                      && mElementType != ElementType::EMPTY )
            {
                if( tNumDim == 0 )
                {
                    tNumDim = mesh::dimension( mElementType ) + 1 ;
                }
                this->work_X().set_size( mesh::number_of_nodes( mElementType ), tNumDim );
                this->work_J().set_size( tNumDim, tNumDim );
            }

        }

 //------------------------------------------------------------------------------

        void
        SideSet::set_thin_shells( Cell< Material* > & aMaterials,
                         const Vector< real > & aLayerThicknesses )
        {
            BELFEM_ERROR( aMaterials.size() == aLayerThicknesses.length() ,
                          "Number of material labels must match number of layer thicknesses ( %u vs %u )",
                          ( unsigned int ) aMaterials.size(), ( unsigned  int ) aLayerThicknesses.length() );

            // determine interpolation order
            uint tInterpolationOrder = mesh::interpolation_order_numeric( this->element_type() );

            // get number of layers
            mNumberOfThinShellLayers = aMaterials.size() ;
            mNumberOfGhostSideSets = tInterpolationOrder * mNumberOfThinShellLayers + 1 ;

            // allocate memory
            mThinShellMaterials.set_size( mNumberOfThinShellLayers, nullptr );

            // get pointer to kernel
            Kernel * tKernel = this->parent()->parent() ;

            mThinShellThickness = 0.0 ;

            for( uint k=0; k<mNumberOfThinShellLayers; ++k )
            {
                // get material
                Material * tMaterial = aMaterials( k );

                // check if material has been added already
                if( ! tMaterial->is_flagged() )
                {
                    // add material to kernel
                    tKernel->add_material( tMaterial );

                    // flag material
                    tMaterial->flag() ;
                }
                mThinShellMaterials( k ) = tMaterial ;
                mThinShellThickness += aLayerThicknesses( k );
            }
            mThinShellThicknesses = aLayerThicknesses ;
        }

//------------------------------------------------------------------------------
    }
}