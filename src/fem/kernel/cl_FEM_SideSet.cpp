/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */
#include <cmath>

#include "cl_Logger.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_IWG.hpp"
#include "commtools.hpp"
#include "meshtools.hpp"

#include "cl_FEM_Element.hpp"

#include "fn_IF_initialize_integration_points.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        bool
        pin_dirichlet_dof(
                const DofManagerBase * aParent,
                      Dof            * aDof,
                const real             aValue,
                const id_t             aGroupID,
                const id_t             aNodeID,
                      index_t        & aNumFirstFlips )
        {
            // only the single-source unit-weight condensation shape is
            // safely reroutable ( see header ); everything else is the
            // caller's business
            if ( ! ( aDof->is_hanging() && aDof->number_of_sources() == 1 ) )
            {
                return false ;
            }

            const real tWeight = aDof->weight( 0 );

            // the weight of a reroutable pair is a stored literal 1.0 with
            // no arithmetic on it, so the tight machine tolerance is safe;
            // anything off unit is a different condensation shape and must
            // not be pinned
            BELFEM_ERROR( std::isfinite( tWeight )
                          && std::abs( tWeight - 1.0 ) <= BELFEM_EPSILON,
                "Dirichlet group %lu ( node %lu ): the dof hangs on one source with weight %g != 1 and cannot be pinned",
                ( long unsigned int ) aGroupID,
                ( long unsigned int ) aNodeID,
                ( double ) tWeight );

            Dof * tSource = aDof->source( 0 );

            BELFEM_ERROR( tSource != nullptr && ! tSource->is_hanging(),
                "Dirichlet group %lu ( node %lu ): the source dof is null or itself hanging - condensation chains are flattened by construction, this is an internal error",
                ( long unsigned int ) aGroupID,
                ( long unsigned int ) aNodeID );

            if ( ! tSource->is_fixed() )
            {
                // flipping a free dof to fixed after initialize() has frozen
                // the free/fixed split would corrupt a classification the
                // solver containers are sized for. Factory-time impositions
                // land before the freeze; a first flip arriving later is a
                // wiring error and must be loud
                BELFEM_ERROR( ! aParent->is_initialized(),
                    "Dirichlet group %lu ( node %lu ): rerouting to the source dof would flip a free dof to fixed after initialization - impose this condition once before initialize()",
                    ( long unsigned int ) aGroupID,
                    ( long unsigned int ) aNodeID );

                ++aNumFirstFlips ;
            }

            // the alias reconstructs as weight * source; divide by the
            // ( unit-tolerance ) weight to land exactly on aValue
            tSource->fix( aValue / tWeight );

            return true ;
        }

//------------------------------------------------------------------------------

        SideSet::SideSet( DofManager * aParent,
                          const id_t aID,
                          Cell< mesh::Facet * > & aFacets,
                          const GroupType aGroupType  ) :
                Group( aParent, aGroupType,
                       aFacets.size() > 0
                       ? aFacets( 0 )->element()->type()
                       : ( ElementType::UNDEFINED ),
                       aID, aFacets.size() ),
                       mMasterType( aFacets.size() == 0 ? ElementType::EMPTY :
                         aFacets( 0 )->master()->type() ),
                       mSlaveType( aFacets.size() == 0 ? ElementType::EMPTY :
                       ( aFacets( 0 )->has_slave() ? aFacets( 0 )->slave()->type() : ElementType::EMPTY ) )
        {
            // inherit the sideset type if the sideset exists on the mesh
            if ( aParent->mesh()->sideset_exists( aID ) )
            {
                this->set_domain_type( aParent->mesh()->sideset( aID )->domain_type() );
            }

            // get the number of dofs per node
            uint tNumDofTypes =  aParent->iwg()->dof_entity_types().length() ;

            // allocate set-wise boundary conditions
            mBcValues.set_size( tNumDofTypes, 0.0 );
            mBcTypes.set_size( tNumDofTypes, BoundaryConditionImposing::Free );

            if( mElementType != ElementType::UNDEFINED )
            {
                // set the number of nodes per element
                mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType );

                this->initialize_elements( aFacets );
                this->collect_nodes( aFacets );
                this->create_element_map();

                mIntegrationOrder = 0 ;
                if( aParent != nullptr )
                {
                    mIntegrationOrder = aParent->sideset_integration_order() ;
                }

                this->initialize_lookup_tables( mIntegrationOrder );

                if( mCalc != nullptr )
                {
                    if ( mParent->iwg() == nullptr )
                    {
                        mCalc->initialize_integration( mElementType, InterpolationType::LAGRANGE );
                    }
                    else
                    {
                        mCalc->initialize_integration( mElementType, aParent->iwg()->interpolation_type());
                    }

                    mCalc->set_integration_order( mIntegrationOrder );
                    mCalc->link( this );
                }
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
            if ( mParent != nullptr )
            {
                this->initialize_lookup_tables( 0 ) ;
            }
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

                for ( IntegrationData * tData : mEnrichmentData )
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

            DofManager * tParent = reinterpret_cast< DofManager * >( mParent );

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
                // get master ID
                id_t tMasterBlockID = 0 ;
                id_t tSlaveBlockID = 0 ;

                if( tBlockMap.key_exists( tFacet->master()->id() ) )
                {
                    tMasterBlockID = tBlockMap[ tFacet->master()->id() ] ;
                }

                if( tFacet->has_slave() )
                {
                    // check if slave block is selected
                    if( tBlockMap.key_exists( tFacet->slave()->id() ) )
                    {
                        tSlaveBlockID = tBlockMap( tFacet->slave()->id() );
                    }
                }

                Element * tMaster = nullptr ;
                Element * tSlave = nullptr ;

                if( ! tParent->block_exists( tMasterBlockID ) )
                {
                    tMasterBlockID = 0 ;
                }
                else if ( tParent->block( tMasterBlockID )->element_exists( tFacet->master()->id() ) )
                {
                    tMaster = tParent->block( tMasterBlockID )->element( tFacet->master()->id() );
                }

                if( ! tParent->block_exists( tSlaveBlockID ) )
                {
                    tSlaveBlockID = 0 ;
                }
                else if ( tParent->block( tSlaveBlockID )->element_exists( tFacet->slave()->id() ) )
                {
                    tSlave = tParent->block( tSlaveBlockID )->element( tFacet->slave()->id() );
                }

                mElements( tCount++ ) = new Element(
                        this,
                        tParent,
                        tFacet,
                        tMode,
                        tMaster,
                        tSlave );
            }
        }

//------------------------------------------------------------------------------

        void
        SideSet::collect_nodes( Cell< mesh::Facet * > & aFacets )
        {
            // reset node container
            mParent->mesh()->unflag_all_nodes();

            // flag all nodes that belong to this set
            for ( mesh::Facet * tFacet : aFacets )
            {
                tFacet->flag_nodes();
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

            index_t tNumFirstFlips = 0 ;

            // loop over all nodes of this mesh
            for ( mesh::Node * tNode : mNodes )
            {
                Dof * tDof = mParent->dof(
                        mParent->calculate_dof_id( tNode, aDofType ) );

                // A condensed ( hanging ) dof is eliminated from the system;
                // fix() on it is historically ineffective for elimination —
                // the T-matrix reconstructs its value from the sources. For a
                // NON-DUPLICATE node with a single-source unit-weight
                // condensation the source sits outside this sideset and the
                // constraint would be lost silently — the shared
                // helper pins the source instead ( same rule as the Maxwell
                // factory pre-fix; Bearing::impose_dirichlet is the
                // template ). Everything else keeps the historical
                // behaviour: a multi-source condensation ( e.g. a refinement
                // hanging node, whose sources this loop typically fixes
                // itself ) is indistinguishable from other provenances at
                // the Dof level, and duplicate-pair condensations ( cuts,
                // thin shells ) can carry an inhomogeneous relation that a
                // plain source pin would violate
                if ( tNode->is_duplicate()
                     || ! pin_dirichlet_dof( mParent, tDof, aValue,
                                             this->id(), tNode->id(),
                                             tNumFirstFlips ) )
                {
                    tDof->fix( aValue );
                }
            }

            // announce once per sideset, not per node; only first flips
            // count, so the per-timestep Dirichlet path stays quiet after
            // its first imposition
            if ( tNumFirstFlips > 0 && comm_rank() == 0 )
            {
                message( InfoLevel::Default,
                    "    SideSet %lu: %lu node(s) are condensed onto sources outside the free set; pinning those source dofs instead",
                    ( long unsigned int ) this->id(),
                    ( long unsigned int ) tNumFirstFlips );
            }

            // set domain type of this shell
            this->set_domain_type( DomainType::Dirichlet );

            // Dirichlet sidesets need calculators for geometry but no DOF allocation
            this->set_activation_mode( GroupActivationMode::GeometryOnly );
        }

//------------------------------------------------------------------------------

        void
        SideSet::impose_neumann( const real aValue, const uint aDofType )
        {
            // set value
            mBcValues( aDofType ) = aValue ;

            mBcTypes( aDofType ) = BoundaryConditionImposing::Neumann ;

            this->set_domain_type( DomainType::Neumann );

            this->activate( true );
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

            // create field for reference temperature
            if( ! mParent->mesh()->field_exists( "Tinf") )
            {
                Vector< real > & tTinf = mParent->mesh()->create_field( "Tinf" ) ;
                tTinf.fill( mTinf );
            }

            // set domaon type of this shell
            this->set_domain_type( DomainType::ThermalAlpha );

            this->activate( true );
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

            this->set_domain_type( DomainType::Default );

            this->activate( false );
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

                if ( mParent != nullptr )
                {
                    if( mParent->iwg()->enrich_sidesets() )
                    {
                        InterpolationFunctionFactory tFactory ;

                        for( IntegrationData * tData : mEnrichmentData )
                        {
                            delete tData ;
                        }
                        mEnrichmentData.set_size( tNumFacets, nullptr );

                        for( uint f=0; f<tNumFacets; ++f )
                        {
                            InterpolationFunction * tBubble = tFactory.create_bubble_function( mMasterType, f );

                            mEnrichmentData( f ) = new IntegrationData( mMasterType, tBubble, true );
                            mEnrichmentData( f )->populate_for_master(
                                f , tIntegrationOrder, tScheme );
                        }
                    }
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
    }
}
