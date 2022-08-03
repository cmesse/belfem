//
// Created by Christian Messe on 26.10.19.
//

#include "cl_FEM_Element.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_DofManager.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element::Element(  Group * aParent, Field * aField, mesh::Element * aElement  ) :
            mParent( aParent ),
            mElement( aElement )
        {
            this->link_dofs( aField );

            // link the extended Jacobian
            this->link_second_derivative_functions( aElement->type() );

        }

//------------------------------------------------------------------------------

        Element::Element(
                Block         * aParent,
                DofManager    * aDofManager,
                mesh::Element * aElement ) :
                mParent( aParent ),
                mElement( aElement )
        {
            // grab ID from group
            this->link_dofs( aDofManager, aParent->id() );

            // link the extended Jacobian
            this->link_second_derivative_functions( aElement->type() );
        }

//------------------------------------------------------------------------------

        Element::Element(
                SideSet       * aParent,
                DofManager    * aDofManager,
                mesh::Facet   * aFacet,
                const SideSetDofLinkMode aSideSetDofLinkMode,
                const id_t      aMasterBlockID,
                const id_t      aSlaveBlockID ) :
                mParent( aParent ),
                mMaster( aFacet->has_master() ?
                    aDofManager->block( aMasterBlockID )->element( aFacet->master()->id()  ) : nullptr ),
                mSlave( aFacet->has_slave() ?
                    aDofManager->block( aSlaveBlockID )->element( aFacet->slave()->id() )  : nullptr ),
                mFacet( aFacet ),
                mElement( aFacet->element() )
        {
            this->link_dofs( aDofManager,
                             aSideSetDofLinkMode,
                             aParent->id(),
                             aMasterBlockID,
                             aSlaveBlockID );

            // link the extended Jacobian
            this->link_second_derivative_functions( aFacet->element()->type() );
        }

//------------------------------------------------------------------------------

        Element::Element( SideSet       * aParent,
                        DofManager    * aDofManager,
                                mesh::Facet   * aFacet,
                        Cell< mesh::Facet * > & aLayers,
                        const id_t      aMasterBlockID,
                        const id_t      aSlaveBlockID ) :
                mParent( aParent ),
                mMaster( aFacet->has_master() ?
                         aDofManager->block( aMasterBlockID )->element( aFacet->master()->id()  ) : nullptr ),
                mSlave( aFacet->has_slave() ?
                        aDofManager->block( aSlaveBlockID )->element( aFacet->slave()->id() )  : nullptr ),
                mFacet( aFacet ),
                mElement( aFacet->element() )
        {

            this->link_dofs_thin_shell( aDofManager,
                                        aLayers,
                                        aDofManager->iwg()->sideset_dof_link_mode(),
                                        aMasterBlockID,
                                        aSlaveBlockID );

        }

//------------------------------------------------------------------------------

        Element::~Element()
        {
            // free rotation data if they have been set
            if( mRotationData != nullptr )
            {
                free( mRotationData );
            }

            // free dof container
            if( mNumberOfDofs > 0 )
            {
                free( mDOFs );
            }
        }

//------------------------------------------------------------------------------

        // deprecated function
        void
        Element::link_dofs( Field * aField )
        {
            // get number of dofs per node and edge
            index_t tNumberOfDofsPerNode = aField->number_of_dofs_per_node();
            index_t tNumberOfNodesPerElement =
                    mParent->parent()->enforce_linear_interpolation() ?
                    mElement->number_of_corner_nodes() :
                    mElement->number_of_nodes();

            index_t tNumberOfDofsPerEdge = aField->number_of_dofs_per_edge();
            index_t tNumberOfEdgesPerElement = mElement->number_of_edges();

            // compute numner of dofs
            mNumberOfDofs  = tNumberOfDofsPerNode * tNumberOfNodesPerElement
                           + tNumberOfDofsPerEdge * tNumberOfEdgesPerElement ;


            // allocate dof container
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // initialize counter
            index_t tCount = 0;

            // loop over all edges
            if ( mParent->parent()->mesh()->edges_exist() && tNumberOfEdgesPerElement > 0 )
            {
                this->compute_edge_directions();

                // link element with dofs
                if ( tNumberOfDofsPerEdge > 0 && tNumberOfEdgesPerElement > 0 )
                {
                    // link element to edge dofs
                    for ( uint k = 0; k < tNumberOfEdgesPerElement; ++k )
                    {
                        if ( mEdgeDirections.test( k ) )
                        {
                            for ( uint i = 0; i < tNumberOfDofsPerEdge; ++i )
                            {
                                mDOFs[ tCount++ ] = aField->dof(
                                        aField->calculate_dof_id( mElement->edge( k ), i ));
                            }
                        }
                        else
                        {
                            for ( int i = tNumberOfDofsPerEdge - 1; i >= 0; i-- )
                            {
                                mDOFs[ tCount++ ] = aField->dof(
                                        aField->calculate_dof_id( mElement->edge( k ), i ));
                            }
                        }
                    }
                }
            }

            // loop over all nodes
            for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
            {
                for ( uint i = 0; i < tNumberOfDofsPerNode; ++i )
                {
                    mDOFs[ tCount++ ] = aField->dof(
                            aField->calculate_dof_id( mElement->node( k ), i ));
                }
            }

            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element with DOFs" );
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs( DofManager * aDofManager, const id_t aBlockID )
        {
            // node dofs
            const Vector< index_t > & tNodeDofTypes =
                    aDofManager->iwg()->dofs_per_node( aBlockID );

            // edge dofs
            const Vector< index_t > & tEdgeDofTypes =
                    aDofManager->iwg()->dofs_per_edge( aBlockID );

            // face dofs
            const Vector< index_t > & tFaceDofTypes =
                    aDofManager->iwg()->dofs_per_face( aBlockID );

            // get number of dofs per node and edge
            index_t tNumberOfDofsPerNode = tNodeDofTypes.length();

            index_t tNumberOfNodesPerElement =
                    aDofManager->enforce_linear_interpolation() ?
                    mElement->number_of_corner_nodes() :
                    mElement->number_of_nodes();

            index_t tNumberOfDofsPerEdge =  tEdgeDofTypes.length() * aDofManager->iwg()->edge_multiplicity();

            index_t tNumberOfDofsPerFace = tFaceDofTypes.length() * aDofManager->iwg()->face_multiplicity();

            index_t tNumberOfEdgesPerElement = mElement->number_of_edges();
            index_t tNumberOfFacesPerElement = mElement->number_of_faces() ;

            mNumberOfDofs = tNumberOfDofsPerNode * tNumberOfNodesPerElement
                          + tNumberOfDofsPerEdge * tNumberOfEdgesPerElement
                          + tNumberOfDofsPerFace * tNumberOfFacesPerElement;

            // allocate dof container
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // initialize counter
            index_t tCount = 0;

            BELFEM_ASSERT( tEdgeDofTypes.length() < 2 && tFaceDofTypes.length() < 2 , "Too many Edge and Face dofs" );

            if ( mElement->has_edges() && tNumberOfEdgesPerElement > 0 )
            {
                this->compute_edge_directions();

                // link element with dofs
                if ( tNumberOfDofsPerEdge > 0 && tNumberOfEdgesPerElement > 0 )
                {
                    // link element to edge dofs
                    for ( uint k = 0; k < tNumberOfEdgesPerElement; ++k )
                    {
                        if ( mEdgeDirections.test( k ) > 0 )
                        {
                            for ( uint i = 0; i < tNumberOfDofsPerEdge; ++i )
                            {
                                mDOFs[ tCount++ ] = aDofManager->dof(
                                        aDofManager->calculate_dof_id( mElement->edge( k ),
                                                                       tEdgeDofTypes( 0 ) + i));
                            }
                        }
                        else
                        {
                            for ( int i = tNumberOfDofsPerEdge - 1; i >= 0; i-- )
                            {
                                mDOFs[ tCount++ ] = aDofManager->dof(
                                        aDofManager->calculate_dof_id( mElement->edge( k ), tEdgeDofTypes( 0 ) + i ) );
                            }
                        }
                    }
                }
            }

            if ( mElement->has_faces() && tNumberOfFacesPerElement > 0 )
            {
                // face orientation is handled by interpolation function
                for ( uint k = 0; k < tNumberOfFacesPerElement; ++k )
                {
                    for ( uint i = 0; i < tNumberOfDofsPerFace; ++i )
                    {
                        mDOFs[ tCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id( mElement->face( k ),
                                                               tFaceDofTypes( 0 ) + i));
                    }
                }
            }

            // here is where we would link #celldof

            // check if node dofs exist
            if( tNumberOfDofsPerNode > 0 )
            {
                // loop over all nodes
                for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
                {

                    for ( index_t tDofType: tNodeDofTypes )
                    {
                        mDOFs[ tCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id( mElement->node( k ), tDofType ));
                    }
                }
            }
            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element %lu with DOFs ( counted %u dofs, expect %u )",
                          ( long unsigned int ) mElement->id(),
                          ( unsigned int ) tCount,
                          ( unsigned int ) mNumberOfDofs );

        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs(
                DofManager                * aDofManager,
                const SideSetDofLinkMode    aMode,
                const id_t                  aSideSetID,
                const id_t                  aMasterBlockID,
                const id_t                  aSlaveBlockID  )
        {
            // get equation object
            IWG * tIWG = aDofManager->iwg() ;

            // if facet has no slave, there are no dofs
            // in that case, we link against an empty vector
            Vector< index_t > tEmpty ;

            // node dofs
            const Vector< index_t > & tMasterNodeDofTypes =
                    aMasterBlockID == 0 ? tEmpty :
                    aDofManager->iwg()->dofs_per_node( aMasterBlockID );

            const Vector< index_t > & tSlaveNodeDofTypes =
                    aSlaveBlockID == 0 ? tEmpty :
                    tIWG->dofs_per_node( aSlaveBlockID );

            // edge dofs
            const Vector< index_t > & tMasterEdgeDofTypes =
                    aMasterBlockID == 0 ? tEmpty :
                    tIWG->dofs_per_edge( aMasterBlockID );

            const Vector< index_t > & tSlaveEdgeDofTypes =
                    aSlaveBlockID == 0 ? tEmpty :
                    tIWG->dofs_per_edge( aSlaveBlockID );

            BELFEM_ASSERT( ! (  tMasterEdgeDofTypes.length() > 0 && tSlaveEdgeDofTypes.length() > 0 ),
                          "can't have edge dofs on both master and slave" );

            // face dofs
            const Vector< index_t > & tMasterFaceDofTypes =
                    aMasterBlockID == 0 ? tEmpty :
                    tIWG->dofs_per_face( aMasterBlockID );

            const Vector< index_t > & tSlaveFaceDofTypes =
                    aSlaveBlockID == 0 ? tEmpty :
                    tIWG->dofs_per_face( aSlaveBlockID );

            BELFEM_ASSERT( ! (  tMasterFaceDofTypes.length() > 0 && tSlaveFaceDofTypes.length() > 0 ),
                          "can't have face dofs on both master and slave" );

            // lagrangian dofs
            const Vector<  index_t > & tLambdaDofTypes = aDofManager->iwg()->lambda_dofs( aSideSetID );

            // catch case for boundaries
            SideSetDofLinkMode    tMode = aMode ;
            if ( aMode == SideSetDofLinkMode::MasterAndSlave && aSlaveBlockID == 0 )
            {
                tMode = aMasterBlockID == 0 ?
                          SideSetDofLinkMode::FacetOnly
                        : SideSetDofLinkMode::FacetAndMaster ;
            }

            // need to do something fancy here for #Facedof
            switch ( tMode )
            {
                case( SideSetDofLinkMode::Cut ) :
                {
                    this->link_dofs_facet_only(
                            aDofManager,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tMasterEdgeDofTypes,
                            tSlaveEdgeDofTypes,
                            tMasterFaceDofTypes,
                            tSlaveFaceDofTypes,
                            tLambdaDofTypes );
                    break ;
                }
                case( SideSetDofLinkMode::FacetOnly ) :
                {
                    this->link_dofs_facet_only(
                            aDofManager,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tMasterEdgeDofTypes,
                            tSlaveEdgeDofTypes,
                            tMasterFaceDofTypes,
                            tSlaveFaceDofTypes,
                            tLambdaDofTypes );

                    break ;
                }
                case( SideSetDofLinkMode::FacetAndMaster ) :
                {

                    this->link_dofs_facet_and_master(
                            aDofManager,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tMasterEdgeDofTypes,
                            tSlaveEdgeDofTypes,
                            tMasterFaceDofTypes,
                            tSlaveFaceDofTypes,
                            tLambdaDofTypes );
                    break ;
                }
                case( SideSetDofLinkMode::FacetAndSlave ) :
                {
                    this->link_dofs_facet_and_slave(
                            aDofManager,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tMasterEdgeDofTypes,
                            tSlaveEdgeDofTypes,
                            tMasterFaceDofTypes,
                            tSlaveFaceDofTypes,
                            tLambdaDofTypes );
                    break ;
                }
                case( SideSetDofLinkMode::MasterAndSlave ) :
                {
                    this->link_dofs_master_and_slave(
                            aDofManager,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tMasterEdgeDofTypes,
                            tSlaveEdgeDofTypes,
                            tMasterFaceDofTypes,
                            tSlaveFaceDofTypes,
                            tLambdaDofTypes );
                    break;
                }
                case( SideSetDofLinkMode::Shell ) :
                {
                    // pass, we have to do this directly over the factory
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Linking Mode");
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_facet_only(
                DofManager  * aDofManager,
                const Vector< index_t > & aMasterNodeDofTypes,
                const Vector< index_t > & aSlaveNodeDofTypes,
                const Vector< index_t > & aMasterEdgeDofTypes,
                const Vector< index_t > & aSlaveEdgeDofTypes,
                const Vector< index_t > & aMasterFaceDofTypes,
                const Vector< index_t > & aSlaveFaceDofTypes,
                const Vector< index_t > & aLambdaDofTypes )
        {
            // get number of dofs per node
            index_t tNumberOfDofsPerNode =
                      aMasterNodeDofTypes.length()
                    + aSlaveNodeDofTypes.length() ;

            // get number of dofs per edge
            index_t tNumberOfDofsPerEdge =
                    ( aMasterEdgeDofTypes.length()
                    + aSlaveEdgeDofTypes.length() ) * aDofManager->iwg()->edge_multiplicity() ;

            // get number of dofs per face
            index_t tNumberOfDofsPerFace =
                    (  aMasterFaceDofTypes.length()
                     + aSlaveFaceDofTypes.length() ) * aDofManager->iwg()->face_multiplicity() ;

            // get the number of nodes on this facet
            index_t tNumberOfNodesPerElement =
                    aDofManager->enforce_linear_interpolation() ?
                    mElement->number_of_corner_nodes() :
                    mElement->number_of_nodes();

            index_t tNumberOfEdgesPerElement = mElement->number_of_edges() ;
            index_t tNumberOfFacesPerElement = mElement->number_of_faces() ;

            // compute the number of DOFs on this sideset
            mNumberOfDofs =
                      tNumberOfDofsPerNode * tNumberOfNodesPerElement
                    + tNumberOfDofsPerEdge * tNumberOfEdgesPerElement
                    + tNumberOfDofsPerFace * tNumberOfFacesPerElement
                    + aLambdaDofTypes.length() ;

            // allocate dof container
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // reset the counter
            index_t tCount = 0 ;

            // link edge directions
            if( aMasterEdgeDofTypes.length() > 0 )
            {
                this->grab_edge_directions_for_facet( mMaster );
            }
            else if ( aSlaveEdgeDofTypes.length() > 0 )
            {
                this->grab_edge_directions_for_facet( mSlave );
            }

            // link to edge dofs on master
            if( tNumberOfDofsPerEdge > 0 )
            {
                // link dofs on master
                this->link_edge_dofs(
                        aDofManager,
                        aMasterEdgeDofTypes,
                        mElement,
                        tNumberOfEdgesPerElement,
                        tCount );

                // link dofs on slave
                this->link_edge_dofs(
                        aDofManager,
                        aSlaveEdgeDofTypes,
                        mElement,
                        tNumberOfEdgesPerElement,
                        tCount );
            }

            // link face dofs on master
            if( tNumberOfDofsPerFace > 0 )
            {
                // link dofs on master
                this->link_face_dofs(
                        aDofManager,
                        aMasterFaceDofTypes,
                        mElement,
                        tNumberOfFacesPerElement,
                        tCount );

                // link dofs on slave
                this->link_face_dofs(
                        aDofManager,
                        aSlaveFaceDofTypes,
                        mElement,
                        tNumberOfFacesPerElement,
                        tCount );
            }

            // link to node dofs on master
            this->link_node_dofs(
                    aDofManager,
                    aMasterNodeDofTypes,
                    mElement,
                    tNumberOfNodesPerElement,
                    tCount );

            // link to node dofs on slave
            this->link_node_dofs(
                    aDofManager,
                    aSlaveNodeDofTypes,
                    mElement,
                    tNumberOfNodesPerElement,
                    tCount );

            this->link_lambda_dofs( aDofManager, aLambdaDofTypes, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element with DOFs" );
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_facet_and_master(
                DofManager  * aDofManager,
                const Vector< index_t > & aMasterNodeDofTypes,
                const Vector< index_t > & aSlaveNodeDofTypes,
                const Vector< index_t > & aMasterEdgeDofTypes,
                const Vector< index_t > & aSlaveEdgeDofTypes,
                const Vector< index_t > & aMasterFaceDofTypes,
                const Vector< index_t > & aSlaveFaceDofTypes,
                const Vector< index_t > & aLambdaDofTypes )
        {
            BELFEM_ASSERT( mMaster != nullptr, "Master element not linked" );

            index_t tNumberOfNodesOnMaster =
                    aDofManager->enforce_linear_interpolation() ?
                    mMaster->element()->number_of_corner_nodes() :
                    mMaster->element()->number_of_nodes();

            // get number of dofs on the facet
            index_t  tNumberOfNodesOnFacet =
                    aDofManager->enforce_linear_interpolation() ?
                    mElement->number_of_corner_nodes() :
                    mElement->number_of_nodes();


            index_t tNumberOfEdgesOnMaster =
                    mMaster->element()->number_of_edges() ;

            index_t tNumberOfEdgesOnFacet =
                    mElement->number_of_edges();

            index_t tNumberOfFacesOnMaster = mMaster->element()->number_of_faces() ;

            index_t tNumberOfFacesOnFacet = mElement->number_of_faces() ;

            // compute memory for dofs
            mNumberOfDofs =
                      tNumberOfNodesOnMaster * aMasterNodeDofTypes.length()
                    + tNumberOfNodesOnFacet  * aSlaveNodeDofTypes.length()
                    + aDofManager->iwg()->edge_multiplicity() *
                                          ( tNumberOfEdgesOnMaster * aMasterEdgeDofTypes.length()
                                          + tNumberOfEdgesOnFacet  * aSlaveEdgeDofTypes.length() )
                    + aDofManager->iwg()->face_multiplicity() *
                                          ( tNumberOfFacesOnMaster * aMasterFaceDofTypes.length()
                                          + tNumberOfFacesOnFacet  * aSlaveFaceDofTypes.length() )
                    + aLambdaDofTypes.length() ;

            // link edge directions
            if( aMasterEdgeDofTypes.length() > 0 )
            {
                mMaster->edge_directions( mEdgeDirections );
            }
            else if ( aSlaveEdgeDofTypes.length() > 0 )
            {
                this->grab_edge_directions_for_facet( mSlave );
            }

            // allocate memory
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // reset the counter
            index_t tCount = 0 ;

            if(  aDofManager->iwg()->edge_multiplicity() > 0 )
            {
                // link to edge dofs on master
                this->link_edge_dofs( aDofManager,
                                      aMasterEdgeDofTypes,
                                      mMaster->element(),
                                      tNumberOfEdgesOnMaster,
                                      tCount );

                // link to edge dofs on facet
                this->link_edge_dofs( aDofManager,
                                      aSlaveEdgeDofTypes,
                                      mElement,
                                      tNumberOfEdgesOnFacet,
                                      tCount );
            }

            if( aDofManager->iwg()->face_multiplicity() > 0 )
            {
                // link to face dofs on master
                this->link_face_dofs( aDofManager,
                                      aMasterFaceDofTypes,
                                      mMaster->element(),
                                      tNumberOfFacesOnMaster,
                                      tCount );

                // link to face dofs on facet
                this->link_face_dofs( aDofManager,
                                      aSlaveFaceDofTypes,
                                      mElement,
                                      tNumberOfFacesOnFacet,
                                      tCount );
            }

            // link to node dofs on master
            this->link_node_dofs( aDofManager,
                                  aMasterNodeDofTypes,
                                  mMaster->element(),
                                  tNumberOfNodesOnMaster,
                                  tCount ) ;

            // link to node dofs on facet
            this->link_node_dofs( aDofManager,
                                  aSlaveNodeDofTypes,
                                  mElement,
                                  tNumberOfNodesOnFacet,
                                  tCount );

            this->link_lambda_dofs( aDofManager, aLambdaDofTypes, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element with DOFs" );
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_facet_and_slave(
                DofManager  * aDofManager,
                const Vector< index_t > & aMasterNodeDofTypes,
                const Vector< index_t > & aSlaveNodeDofTypes,
                const Vector< index_t > & aMasterEdgeDofTypes,
                const Vector< index_t > & aSlaveEdgeDofTypes,
                const Vector< index_t > & aMasterFaceDofTypes,
                const Vector< index_t > & aSlaveFaceDofTypes,
                const Vector< index_t > & aLambdaDofTypes )
        {
            BELFEM_ASSERT( mSlave != nullptr, "Slave element not linked" );

            // get number of dofs on the facet
            index_t  tNumberOfNodesOnFacet =
                    aDofManager->enforce_linear_interpolation() ?
                    mElement->number_of_corner_nodes() :
                    mElement->number_of_nodes();

            index_t tNumberOfNodesOnSlave =
                    aDofManager->enforce_linear_interpolation() ?
                    mSlave->element()->number_of_corner_nodes() :
                    mSlave->element()->number_of_nodes();

            index_t tNumberOfEdgesOnFacet =
                    mElement->number_of_edges();

            index_t tNumberOfEdgesOnSlave =
                    mSlave->element()->number_of_edges() ;

            index_t tNumberOfFacesOnFacet =
                    mElement->number_of_faces();

            index_t tNumberOfFacesOnSlave =
                    mSlave->element()->number_of_faces() ;

            // compute memory for dofs
            mNumberOfDofs =
                      tNumberOfNodesOnFacet  * aMasterNodeDofTypes.length()
                    + tNumberOfNodesOnSlave  * aSlaveNodeDofTypes.length()
                    + aDofManager->iwg()->edge_multiplicity() *
                        (   tNumberOfEdgesOnFacet  * aMasterEdgeDofTypes.length()
                          + tNumberOfEdgesOnSlave  * aSlaveEdgeDofTypes.length() )
                        + aDofManager->iwg()->face_multiplicity() *
                        (   tNumberOfFacesOnFacet  * aMasterFaceDofTypes.length()
                          + tNumberOfFacesOnSlave  * aSlaveFaceDofTypes.length() )
                    + aLambdaDofTypes.length() ;

            // link edge directions
            if( aMasterEdgeDofTypes.length() > 0 )
            {
                this->grab_edge_directions_for_facet( mMaster );
            }
            else if ( aSlaveEdgeDofTypes.length() > 0 )
            {
                mSlave->edge_directions( mEdgeDirections );
            }

            // allocate memory
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // reset the counter
            index_t tCount = 0 ;

            if( aDofManager->iwg()->edge_multiplicity()  > 0 )
            {
                // link to edge dofs on facet
                this->link_edge_dofs( aDofManager,
                                      aMasterEdgeDofTypes,
                                      mElement,
                                      tNumberOfEdgesOnFacet,
                                      tCount );

                // link to edge dofs on master
                this->link_edge_dofs( aDofManager,
                                      aSlaveEdgeDofTypes,
                                      mSlave->element(),
                                      tNumberOfEdgesOnSlave,
                                      tCount );
            }

            if( aDofManager->iwg()->face_multiplicity() > 0 )
            {
                // link to face dofs on facet
                this->link_face_dofs( aDofManager,
                                      aMasterFaceDofTypes,
                                      mElement,
                                      tNumberOfFacesOnFacet,
                                      tCount );

                // link to face dofs on master
                this->link_face_dofs( aDofManager,
                                      aSlaveFaceDofTypes,
                                      mSlave->element(),
                                      tNumberOfFacesOnSlave,
                                      tCount );
            }

            // link to node dofs on facet
            this->link_node_dofs( aDofManager,
                                  aMasterNodeDofTypes,
                                  mElement,
                                  tNumberOfNodesOnFacet,
                                  tCount );

            // link to node dofs on slave
            this->link_node_dofs( aDofManager,
                                  aSlaveNodeDofTypes,
                                  mSlave->element(),
                                  tNumberOfNodesOnSlave,
                                  tCount ) ;


            this->link_lambda_dofs( aDofManager, aLambdaDofTypes, tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element with DOFs" );
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_master_and_slave(
                DofManager  * aDofManager,
                const Vector< index_t > & aMasterNodeDofTypes,
                const Vector< index_t > & aSlaveNodeDofTypes,
                const Vector< index_t > & aMasterEdgeDofTypes,
                const Vector< index_t > & aSlaveEdgeDofTypes,
                const Vector< index_t > & aMasterFaceDofTypes,
                const Vector< index_t > & aSlaveFaceDofTypes,
                const Vector< index_t > & aLambdaDofTypes )
        {
            BELFEM_ASSERT( mMaster != nullptr, "Master element not linked" );
            BELFEM_ASSERT( mSlave != nullptr, "Slave element not linked" );

            // get number of dofs on the facet
            index_t  tNumberOfNodesOnMaster =
                    aDofManager->enforce_linear_interpolation() ?
                    mMaster->element()->number_of_corner_nodes() :
                    mMaster->element()->number_of_nodes();

            index_t tNumberOfNodesOnSlave =
                    aDofManager->enforce_linear_interpolation() ?
                    mSlave->element()->number_of_corner_nodes() :
                    mSlave->element()->number_of_nodes();

            index_t tNumberOfEdgesOnMaster =
                    mMaster->element()->number_of_edges();

            index_t tNumberOfEdgesOnSlave =
                    mSlave->element()->number_of_edges() ;

            index_t tNumberOfFacesOnMaster =
                    mMaster->element()->number_of_faces();

            index_t tNumberOfFacesOnSlave =
                    mSlave->element()->number_of_faces() ;

            // compute memory for dofs
            mNumberOfDofs =
                      tNumberOfNodesOnMaster * aMasterNodeDofTypes.length()
                    + tNumberOfNodesOnSlave  * aSlaveNodeDofTypes.length()
                    + aDofManager->iwg()->edge_multiplicity() * (
                          tNumberOfEdgesOnMaster * aMasterEdgeDofTypes.length()
                        + tNumberOfEdgesOnSlave  * aSlaveEdgeDofTypes.length() )
                    + aDofManager->iwg()->face_multiplicity() * (
                          tNumberOfFacesOnMaster * aMasterFaceDofTypes.length()
                        + tNumberOfFacesOnSlave  * aSlaveFaceDofTypes.length() )
                    + aLambdaDofTypes.length() ;

            // link edge directions
            if( aMasterEdgeDofTypes.length() > 0 )
            {
                mMaster->edge_directions( mEdgeDirections );
            }
            else if ( aSlaveEdgeDofTypes.length() > 0 )
            {
                mSlave->edge_directions( mEdgeDirections );
            }

            // allocate memory
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // reset the counter
            index_t tCount = 0 ;

            if( aDofManager->iwg()->edge_multiplicity() > 0 )
            {
                // link to edge dofs on master
                this->link_edge_dofs( aDofManager,
                                      aMasterEdgeDofTypes,
                                      mMaster->element(),
                                      tNumberOfEdgesOnMaster,
                                      tCount );

                // link to edge dofs on slave
                this->link_edge_dofs( aDofManager,
                                      aSlaveEdgeDofTypes,
                                      mSlave->element(),
                                      tNumberOfEdgesOnSlave,
                                      tCount );
            }

            if( aDofManager->iwg()->face_multiplicity() > 0 )
            {
                // link to face dofs on master
                this->link_face_dofs( aDofManager,
                                      aMasterFaceDofTypes,
                                      mMaster->element(),
                                      tNumberOfFacesOnMaster,
                                      tCount );

                // link to face dofs on skave
                this->link_face_dofs( aDofManager,
                                      aSlaveFaceDofTypes,
                                      mSlave->element(),
                                      tNumberOfFacesOnSlave,
                                      tCount );
            }

            // link to node dofs on facet
            this->link_node_dofs( aDofManager,
                                  aMasterNodeDofTypes,
                                  mMaster->element(),
                                  tNumberOfNodesOnMaster,
                                  tCount );

            // link to node dofs on slave
            this->link_node_dofs( aDofManager,
                                  aSlaveNodeDofTypes,
                                  mSlave->element(),
                                  tNumberOfNodesOnSlave,
                                  tCount ) ;

            // for contact etc
            this->link_lambda_dofs( aDofManager,
                    aLambdaDofTypes,
                    tCount );

            BELFEM_ASSERT( tCount == mNumberOfDofs,
                          "Something went wrong while linking element with DOFs" );

        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_thin_shell(
                DofManager  * aDofManager,
                Cell< mesh::Facet * > & aLayers ,
                const SideSetDofLinkMode aMode,
                const id_t aMasterBlockID,
                const id_t aSlaveBlockID )
        {
            // help vector
            Vector< index_t > tEmpty ;

            // for this method, we assume that only node dofs exist on
            // master and slave
            const Vector< index_t > & tMasterNodeDofTypes =
                    aMasterBlockID == 0 ? tEmpty :
                    aDofManager->iwg()->dofs_per_node( aMasterBlockID );

            const Vector< index_t > & tSlaveNodeDofTypes =
                    aSlaveBlockID == 0 ? tEmpty :
                    aDofManager->iwg()->dofs_per_node( aSlaveBlockID );

            const Vector< index_t > & tShellNodeDofTypes = aDofManager->iwg()->dofs_per_node_on_sideset( mParent->id() );
            const Vector< index_t > & tShellEdgeDofTypes = aDofManager->iwg()->dofs_per_edge_on_sideset( mParent->id() );
            const Vector< index_t > & tShellFaceDofTypes = aDofManager->iwg()->dofs_per_face_on_sideset( mParent->id() );
            const Vector< index_t > & tLambdaDofTypes    = aDofManager->iwg()->lambda_dofs( mParent->id() );

            if( tShellEdgeDofTypes.length() > 0 )
            {
                this->compute_edge_directions() ;
            }
            switch( aMode )
            {
                case ( SideSetDofLinkMode::MasterAndSlave ) :
                {
                    this->link_dofs_thin_shell_master_and_slave(
                            aDofManager,
                            aLayers,
                            tMasterNodeDofTypes,
                            tSlaveNodeDofTypes,
                            tShellEdgeDofTypes,
                            tShellFaceDofTypes,
                            tLambdaDofTypes );
                    break ;
                }
                case( SideSetDofLinkMode::FacetOnly ) :
                {
                    this->link_dofs_thin_shell_facet_only(
                            aDofManager,
                            aLayers,
                            tShellNodeDofTypes,
                            tShellEdgeDofTypes,
                            tShellFaceDofTypes );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Linking mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_thin_shell_master_and_slave(
                DofManager  * aDofManager,
                Cell< mesh::Facet * > & aLayers ,
                const Vector< index_t > & aMasterNodeDofTypes,
                const Vector< index_t > & aSlaveNodeDofTypes,
                const Vector< index_t > & aThinShellEdgeDofTypes,
                const Vector< index_t > & aThinShellFaceDofTypes,
                const Vector< index_t > & aLambdaDofTypes )
        {
            index_t  tNumberOfNodesOnMaster =
                    aDofManager->enforce_linear_interpolation() ?
                    mMaster->element()->number_of_corner_nodes() :
                    mMaster->element()->number_of_nodes();

            index_t tNumberOfNodesOnSlave =
                    aDofManager->enforce_linear_interpolation() ?
                    mSlave->element()->number_of_corner_nodes() :
                    mSlave->element()->number_of_nodes();

            index_t tNumberOfEdges = mElement->number_of_edges() ;

            index_t tNumberOfFaces = mElement->number_of_faces() ;

            uint tNumberOfLayers = aLayers.size() ;

            // compute memory for dofs
            mNumberOfDofs =
                      tNumberOfNodesOnMaster * aMasterNodeDofTypes.length()
                    + tNumberOfNodesOnSlave  * aSlaveNodeDofTypes.length()
                    + tNumberOfLayers * (
                    + tNumberOfEdges * aDofManager->iwg()->edge_multiplicity() * aThinShellEdgeDofTypes.length()
                    + tNumberOfFaces * aDofManager->iwg()->face_multiplicity() * aThinShellFaceDofTypes.length() )
                    + aLambdaDofTypes.length() * aDofManager->iwg()->lambda_multiplicity() ;

            // allocate memory
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // initialize the counter
            index_t tCount = 0 ;

            // link node dofs on master
            this->link_node_dofs( aDofManager, aMasterNodeDofTypes, mMaster->element(), tNumberOfNodesOnMaster, tCount );

            // link node dofs on slave
            this->link_node_dofs( aDofManager, aSlaveNodeDofTypes, mSlave->element(), tNumberOfNodesOnSlave, tCount );

            // link dofs per layer
            for( mesh::Facet * tLayer : aLayers )
            {
                // link edge dofs
                for( uint k=0; k<tLayer->number_of_edges(); ++k )
                {
                    for( uint d: aThinShellEdgeDofTypes )
                    {
                        mDOFs[ tCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id(
                                        tLayer->edge( k ), d ) );
                    }
                }

                // link face dofs
                for( uint d: aThinShellFaceDofTypes )
                {
                    mDOFs[ tCount++ ] = aDofManager->dof(
                            aDofManager->calculate_dof_id(
                                    tLayer->face(), d ) );
                }

            }

            // link lambda dofs
            this->link_lambda_dofs( aDofManager, aLambdaDofTypes, tCount );

            BELFEM_ASSERT( mNumberOfDofs == tCount, "number of dofs for element %lu does not match (is %u, expect %u)",
                ( long unsigned int ) this->id(), tCount, mNumberOfDofs );
        }

//------------------------------------------------------------------------------

        void
        Element::link_dofs_thin_shell_facet_only(
                DofManager  * aDofManager,
                Cell< mesh::Facet * > & aLayers,
                const Vector< index_t > & aThinShellNodeDofTypes,
                const Vector< index_t > & aThinShellEdgeDofTypes,
                const Vector< index_t > & aThinShellFaceDofTypes )
        {
            uint tNumberOfLayers = aLayers.size() ;

            index_t tNumberOfNodes = mElement->number_of_nodes() ;

            index_t tNumberOfEdges = mElement->number_of_edges() ;

            index_t tNumberOfFaces = mElement->number_of_faces() ;

            mNumberOfDofs =
                    tNumberOfLayers * (
                            + tNumberOfNodes * aThinShellNodeDofTypes.length()
                            + aDofManager->iwg()->edge_multiplicity() * aThinShellEdgeDofTypes.length() * tNumberOfEdges
                            + aDofManager->iwg()->face_multiplicity() * aThinShellFaceDofTypes.length() * tNumberOfFaces ) ;

            // allocate memory
            mDOFs = ( Dof ** ) malloc( mNumberOfDofs * sizeof( Dof ) );

            // initialize the counter
            index_t tCount = 0 ;

            // link dofs per layer
            for( mesh::Facet * tLayer : aLayers )
            {
                // link the node dofs
                this->link_node_dofs( aDofManager, aThinShellNodeDofTypes, tLayer->element(), tNumberOfNodes, tCount );

                // link edge dofs
                for( uint k=0; k<tLayer->number_of_edges(); ++k )
                {
                    for( uint d: aThinShellEdgeDofTypes )
                    {
                        std::cout << "check layer dof " << tLayer->edge( k )->id() << " " << d << " " << aDofManager->calculate_dof_id(
                                tLayer->edge( k ), d ) << std::endl ;

                        mDOFs[ tCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id(
                                        tLayer->edge( k ), d ) );
                    }
                }

                // link face dofs
                for( uint d: aThinShellFaceDofTypes )
                {
                    mDOFs[ tCount++ ] = aDofManager->dof(
                            aDofManager->calculate_dof_id(
                                    tLayer->face(), d ) );
                }

            }

            BELFEM_ASSERT( mNumberOfDofs == tCount, "number of dofs for element %lu does not match (is %u, expect %u)",
                           ( long unsigned int ) this->id(), tCount, mNumberOfDofs );
        }

//------------------------------------------------------------------------------

        void
        Element::grab_edge_directions_for_facet( Element * aReferenceElement )
        {

            // we create a temporary map to assign the signs
            Map< id_t, bool > tSigns;
            for( uint e=0; e<aReferenceElement->element()->number_of_edges(); ++e )
            {
                tSigns[ aReferenceElement->element()->edge( e )->id() ] = aReferenceElement->edge_direction( e ) == 1.0 ;
            }


            for( uint e=0; e<mElement->number_of_edges(); ++e )
            {
                if(  tSigns( mElement->edge( e )->id() ) > 0 )
                {
                    mEdgeDirections.set( e );
                }
                else
                {
                    mEdgeDirections.reset( e );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Element::link_node_dofs(
                DofManager              * aDofManager,
                const Vector< index_t > & aDofTypes,
                mesh::Element           * aReferenceElement,
                const                uint aNumberOfNodesPerElement,
                index_t                 & aCount )
        {
            // do nothing if there are no dofs to link
            if ( aDofTypes.length() == 0 )
            {
                return;
            }

            for ( index_t tDofType : aDofTypes )
            {
                for ( uint k = 0; k < aNumberOfNodesPerElement; ++k )
                {
                    mDOFs[ aCount++ ] = aDofManager->dof(
                            aDofManager->calculate_dof_id(
                                    aReferenceElement->node( k ), tDofType ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::link_edge_dofs(
                       DofManager              * aDofManager,
                        const Vector< index_t > & aDofTypes,
                        mesh::Element           * aReferenceElement,
                        const                uint aNumberOfEdgesPerElement,
                        index_t                 & aCount )
        {
            // do nothing if there are no dofs to link
            if ( aDofTypes.length() == 0 )
            {
                return;
            }

            BELFEM_ASSERT( aDofTypes.length() == 1, "can only have one edge dof type" );

            // get the number of dof types, we must do this the old fashioned way
            // since the edge might in backwards direction
            int tNumberOfDofsPerEdge = aDofManager->iwg()->edge_multiplicity() ;

            // get multiplicity of edges
            uint tOff = aDofTypes( 0 ) ;

            for ( uint k = 0; k < aNumberOfEdgesPerElement; ++k )
            {
                if ( mEdgeDirections.test( k ) )
                {
                    for ( int i = 0; i < tNumberOfDofsPerEdge; ++i )
                    {
                        mDOFs[ aCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id( aReferenceElement->edge( k ), tOff + i ) );
                    }
                }
                else
                {
                    for ( int i = tNumberOfDofsPerEdge - 1; i >= 0; i-- )
                    {
                        mDOFs[ aCount++ ] = aDofManager->dof(
                                aDofManager->calculate_dof_id( aReferenceElement->edge( k ), tOff + i  ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::link_face_dofs(
                DofManager  * aDofManager,
                const Vector< index_t > & aDofTypes,
                mesh::Element * aReferenceElement,
                const uint aNumberOfFacesPerElement,
                index_t & aCount )
        {
            // do nothing if there are no dofs to link
            if ( aDofTypes.length() == 0 )
            {
                return;
            }

            BELFEM_ASSERT( aDofTypes.length() == 1, "can only have one face dof type" );

            uint tNumberOfDofsPerFace = aDofManager->iwg()->face_multiplicity() ;
            uint tOff = aDofTypes( 0 ) ;
            for ( uint k = 0; k < aNumberOfFacesPerElement; ++k )
            {
                for ( uint i = 0; i < tNumberOfDofsPerFace; ++i )
                {
                    mDOFs[ aCount++ ] = aDofManager->dof(
                            aDofManager->calculate_dof_id( aReferenceElement->face( k ), tOff + i ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::link_lambda_dofs(   DofManager  * aDofManager,
                            const Vector< index_t > & aDofTypes,
                            index_t & aCount )
        {
            // do nothing if there are no dofs to link
            if ( aDofTypes.length() == 0 )
            {
                return;
            }

            BELFEM_ASSERT( mFacet != nullptr, "No Facet connected to this element" );

            for ( index_t tDofType : aDofTypes )
            {
                mDOFs[ aCount++ ] = aDofManager->dof(
                        aDofManager->calculate_dof_id(
                                mFacet, tDofType ) );
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Element::N( const uint aPointIndex ) const
        {
            return mParent->N( aPointIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Element::G( const uint aPointIndex ) const
        {
            return mParent->G( aPointIndex );
        }

//------------------------------------------------------------------------------

        Matrix< real > &
        Element::J( const uint aPointIndex )
        {
            // remember index
            mParent->work_index() = aPointIndex;

            // collect node coordinates
            this->get_node_coors( mParent->node_coords() );

            // compute jacobian
            mParent->work_J() = mParent->dGdXi( aPointIndex )
                    * mParent->node_coords() ;

            // return matrix
            return mParent->work_J();
        }

//------------------------------------------------------------------------------

        Matrix< real > &
        Element::K( const uint aPointIndex )
        {
            BELFEM_ASSERT( mParent->work_index() == aPointIndex,
             "J-Matrix has not been computed for this integration point ( is %u but expect %u )",
             ( unsigned int ) aPointIndex,
             ( unsigned int ) mParent->work_index() );

            // collect node coordinates
            this->get_node_coors( mParent->node_coords() );

            // compute K-Matrix
            mParent->work_K() = mParent->d2GdXi2( aPointIndex ) * mParent->node_coords();

            // return matrix
            return mParent->work_K();
        }
//------------------------------------------------------------------------------

        Matrix< real > &
        Element::L( const uint aPointIndex )
        {
            BELFEM_ASSERT( mParent->work_index() == aPointIndex,
                          "J-Matrix has not been computed for this integration point ( is %u but expect %u )",
                                  ( unsigned int ) aPointIndex,
                          ( unsigned int ) mParent->work_index() );

            // compute K-Matrix
            ( this->*mL )( mParent->work_J(), mParent->work_L() );

            // return matrix
            return mParent->work_L();
        }

//------------------------------------------------------------------------------

        void
        Element::L1D( Matrix< real > & aJ, Matrix< real > & aL )
        {
            aL( 0, 0 ) = aJ( 0, 0 ) * aJ( 0, 0 );
        }

//------------------------------------------------------------------------------

        void
        Element::L2D( Matrix< real > & aJ, Matrix< real > & aL )
        {
            aL( 0, 0 ) = aJ( 0, 0 ) * aJ( 0, 0 );
            aL( 1, 0 ) = aJ( 1, 0 ) * aJ( 1, 0 );
            aL( 2, 0 ) = aJ( 0, 0 ) * aJ( 1, 0 );

            aL( 0, 1 ) = aJ( 0, 1 ) * aJ( 0, 1 );
            aL( 1, 1 ) = aJ( 1, 1 ) * aJ( 1, 1 );
            aL( 2, 1 ) = aJ( 0, 1 ) * aJ( 1, 1 );

            aL( 0, 2 ) = 2.0 * aJ( 0, 0 ) * aJ( 0, 1 );
            aL( 1, 2 ) = 2.0 * aJ( 1, 0 ) * aJ( 1, 1 );
            aL( 2, 2 ) = aJ( 0, 0 ) * aJ( 1, 1 ) + aJ( 1, 0 ) * aJ( 0, 1 );
        }

//------------------------------------------------------------------------------

        void
        Element::L3D( Matrix< real > & aJ, Matrix< real > & aL )
        {
            aL( 0, 0 ) = aJ( 0, 0 ) * aJ( 0, 0 );
            aL( 1, 0 ) = aJ( 1, 0 ) * aJ( 1, 0 );
            aL( 2, 0 ) = aJ( 2, 0 ) * aJ( 2, 0 );
            aL( 3, 0 ) = aJ( 1, 0 ) * aJ( 2, 0 );
            aL( 4, 0 ) = aJ( 0, 0 ) * aJ( 2, 0 );
            aL( 5, 0 ) = aJ( 0, 0 ) * aJ( 1, 0 );

            aL( 0, 1 ) = aJ( 0, 1 ) * aJ( 0, 1 );
            aL( 1, 1 ) = aJ( 1, 1 ) * aJ( 1, 1 );
            aL( 2, 1 ) = aJ( 2, 1 ) * aJ( 2, 1 );
            aL( 3, 1 ) = aJ( 1, 1 ) * aJ( 2, 1 );
            aL( 4, 1 ) = aJ( 0, 1 ) * aJ( 2, 1 );
            aL( 5, 1 ) = aJ( 0, 1 ) * aJ( 1, 1 );

            aL( 0, 2 ) = aJ( 0, 2 ) * aJ( 0, 2 );
            aL( 1, 2 ) = aJ( 1, 2 ) * aJ( 1, 2 );
            aL( 2, 2 ) = aJ( 2, 2 ) * aJ( 2, 2 );
            aL( 3, 2 ) = aJ( 1, 2 ) * aJ( 2, 2 );
            aL( 4, 2 ) = aJ( 0, 2 ) * aJ( 2, 2 );
            aL( 5, 2 ) = aJ( 0, 2 ) * aJ( 1, 2 );

            aL( 0, 3 ) = 2.0 * aJ( 0, 1 ) * aJ( 0, 2 );
            aL( 1, 3 ) = 2.0 * aJ( 1, 1 ) * aJ( 1, 2 );
            aL( 2, 3 ) = 2.0 * aJ( 2, 1 ) * aJ( 2, 2 );
            aL( 3, 3 ) = aJ( 1, 1 ) * aJ( 2, 2 ) + aJ( 2, 1 ) * aJ( 1, 2 );
            aL( 4, 3 ) = aJ( 0, 1 ) * aJ( 2, 2 ) + aJ( 2, 1 ) * aJ( 0, 2 );
            aL( 5, 3 ) = aJ( 0, 1 ) * aJ( 1, 2 ) + aJ( 1, 1 ) * aJ( 0, 2 );

            aL( 0, 4 ) = 2.0 * aJ( 0, 0 ) * aJ( 0, 2 );
            aL( 1, 4 ) = 2.0 * aJ( 1, 0 ) * aJ( 1, 2 );
            aL( 2, 4 ) = 2.0 * aJ( 2, 0 ) * aJ( 2, 2 );
            aL( 3, 4 ) = aJ( 1, 0 ) * aJ( 2, 2 ) + aJ( 2, 0 ) * aJ( 1, 2 );
            aL( 4, 4 ) = aJ( 0, 0 ) * aJ( 2, 2 ) + aJ( 2, 0 ) * aJ( 0, 2 );
            aL( 5, 4 ) = aJ( 0, 0 ) * aJ( 1, 2 ) + aJ( 1, 0 ) * aJ( 0, 2 );

            aL( 0, 5 ) = 2.0 * aJ( 0, 0 ) * aJ( 0, 1 );
            aL( 1, 5 ) = 2.0 * aJ( 1, 0 ) * aJ( 1, 1 );
            aL( 2, 5 ) = 2.0 * aJ( 2, 0 ) * aJ( 2, 1 );
            aL( 3, 5 ) = aJ( 1, 0 ) * aJ( 2, 1 ) + aJ( 2, 0 ) * aJ( 1, 1 );
            aL( 4, 5 ) = aJ( 0, 0 ) * aJ( 2, 1 ) + aJ( 2, 0 ) * aJ( 0, 1 );
            aL( 5, 5 ) = aJ( 0, 0 ) * aJ( 1, 1 ) + aJ( 1, 0 ) * aJ( 0, 1 );
        }

//------------------------------------------------------------------------------

        void
        Element::flag_dofs()
        {
            for( uint d=0; d<mNumberOfDofs; ++d )
            {
                mDOFs[ d ]->flag();
            }
        }

//------------------------------------------------------------------------------

        void
        Element::unflag_dofs()
        {
            for( uint d=0; d<mNumberOfDofs; ++d )
            {
                mDOFs[ d ]->unflag();
            }
        }

//------------------------------------------------------------------------------

        void
        Element::compute_edge_directions()
        {

            // number of edges of this element
            uint tNumEdges = mElement->number_of_edges() ;

            if( tNumEdges == 0 )
            {
                return;
            }

            BELFEM_ASSERT( mElement->has_edges(), "edges not allocated for this element %lu",
                           ( long unsigned int ) mElement->id() );

            Cell< mesh::Node * > tNodes ;

            // if this is a facet on a tape,
            // make sure that we take the original facet
            //mesh::Element * tElement = ( mFacet == nullptr ) ? mElement :
             //       mParent->parent()->mesh()->original_facet( mElement->id() )->element() ;

            for( uint e=0; e<tNumEdges; ++e )
            {
                // grab node IDs
                id_t tA = mElement->edge( e )->node( 0 )->id();
                id_t tB = mElement->edge( e )->node( 1 )->id();

                mElement->get_nodes_of_edge( e, tNodes );

                if ( tA == tNodes(0)->id() && tB == tNodes(1)->id() )
                {
                    mEdgeDirections.set( e );
                }
                else if ( tA == tNodes(1)->id() && tB == tNodes(0)->id() )
                {
                    mEdgeDirections.reset( e );
                }
                else
                {
                    // check if this is a cut
                    tA = mParent->parent()->mesh()->duplicate_node( tA )->id() ;
                    tB = mParent->parent()->mesh()->duplicate_node( tB )->id() ;

                    bool tSuccess = false ;

                    if ( tA == tNodes(0)->id() && tB == tNodes(1)->id() )
                    {
                        mEdgeDirections.set( e );
                        tSuccess = true ;
                    }
                    else if ( tA == tNodes(1)->id() && tB == tNodes(0)->id() )
                    {
                        mEdgeDirections.reset( e );
                        tSuccess = true ;
                    }

                    BELFEM_ERROR( tSuccess, "Error when trying to compute edge direction of element %lu of type %s",
                                 ( long unsigned int ) this->id(), to_string( mElement->type() ).c_str() );
                }

            }
        }

//------------------------------------------------------------------------------

        void
        Element::link_second_derivative_functions( const ElementType aElementType )
        {
            // link L-Function for second derivatives
            switch( mesh::dimension( aElementType ) )
            {
                case( 1 ):
                {
                    mL  = & Element::L1D;
                    break;
                }
                case( 2 ):
                {
                    mL  = & Element::L2D;
                    break;
                }
                case( 3 ):
                {
                    mL  = & Element::L3D;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid dimension for element type");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::set_rotation_data(
                const Matrix< real > & aRotationAxes,
                const Vector< real > & aRotationAngles )
        {
            // free memory if it has been set
            if( mRotationData != nullptr )
            {
                free( mRotationData );
            }

            // get the number of nodes from the element
            uint tNumNodes = mElement->number_of_nodes() ;

            // make sure that container sizes are correct
            BELFEM_ASSERT( aRotationAxes.n_rows() >= 3,
                          "Invalid number of rows for aRotationAxes. Is %u but expect at least 3",
                          ( unsigned int ) aRotationAngles.length() );

            BELFEM_ASSERT( aRotationAxes.n_cols() >= tNumNodes,
                          "Invalid number of columns for aRotationAngles. Is %u but expect at least %u.",
                          ( unsigned int ) aRotationAngles.length(),
                          ( unsigned int ) tNumNodes );


            BELFEM_ASSERT( aRotationAngles.length() >= tNumNodes,
                          "Invalid length of  aRotationAngles. Is %u but expect at least %u.",
                          ( unsigned int ) aRotationAngles.length(),
                          ( unsigned int ) tNumNodes );

            // allocate the memory
            mRotationData   = ( real * ) malloc( 4 * tNumNodes * sizeof( real ) );

            uint tCount = 0 ;

            // populate data
            for( uint k=0; k<tNumNodes; ++k )
            {
                mRotationData[ tCount++  ] = aRotationAngles( k );
            }

            // initialize counter
            for( uint k=0; k<tNumNodes; ++k )
            {
                mRotationData[ tCount++ ] = aRotationAxes( 0, k );
                mRotationData[ tCount++ ] = aRotationAxes( 1, k );
                mRotationData[ tCount++ ] = aRotationAxes( 2, k );
            }

        }

//------------------------------------------------------------------------------

        void
        Element::get_rotation_data( Matrix< real > & aRotationAxes, Vector< real > & aRotationAngles )
        {
            // get the number of nodes from the element
            uint tNumNodes = mElement->number_of_nodes() ;

            // make sure that container sizes are correct
            BELFEM_ASSERT( aRotationAxes.n_rows() >= 3,
                          "Invalid number of rows for aRotationAxes. Is %u but expect at least 3",
                          ( unsigned int ) aRotationAngles.length() );

            BELFEM_ASSERT( aRotationAxes.n_cols() >= tNumNodes,
                          "Invalid number of columns for aRotationAngles. Is %u but expect at least %u.",
                          ( unsigned int ) aRotationAngles.length(),
                          ( unsigned int ) tNumNodes );


            BELFEM_ASSERT( aRotationAngles.length() >= tNumNodes,
                          "Invalid length of  aRotationAngles. Is %u but expect at least %u.",
                          ( unsigned int ) aRotationAngles.length(),
                          ( unsigned int ) tNumNodes );

            // make sure that fields have been allocated
            BELFEM_ASSERT( mRotationData != nullptr,
                          "Rotation data has not been set for element %lu",
                          ( long unsigned int ) mElement->id() );

            uint tCount = 0 ;

            // populate data
            for( uint k=0; k<tNumNodes; ++k )
            {
                aRotationAngles( k ) = mRotationData[ tCount++  ] ;
            }

            // initialize counter
            for( uint k=0; k<tNumNodes; ++k )
            {
                aRotationAxes( 0, k ) = mRotationData[ tCount++ ] ;
                aRotationAxes( 1, k ) = mRotationData[ tCount++ ] ;
                aRotationAxes( 2, k ) = mRotationData[ tCount++ ] ;
            }

        }

//------------------------------------------------------------------------------
    }
}