//
// Created by Christian Messe on 24.10.19.
//

#ifndef BELFEM_CL_FEM_ELEMENT_HPP
#define BELFEM_CL_FEM_ELEMENT_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Element.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_Group.hpp"
#include "en_IWG_SideSetDofLinkMode.hpp"
#include "cl_Bitset.hpp"

namespace belfem
{
    namespace fem
    {
        class Field;
        class DofManager ;
        class Block;
        class SideSet ;

//------------------------------------------------------------------------------

        class Element
        {
            // parent block of this element
            Group * mParent;

            // pointer to master, if this is a side element
            Element * mMaster = nullptr ;

            // pointer to master, if this is a side element
            Element * mSlave = nullptr ;

            // Facet in mesh, if this is a side element
            mesh::Facet  * mFacet = nullptr ;

            // index on
            // actual element on mesh
            mesh::Element * mElement;

            // contains the dofs that are connected to this element
            // Cell< Dof * >  mDOFs;
            Dof ** mDOFs = nullptr ;
            uint mNumberOfDofs = 0 ;

            // pointer to L-Matrix function
            void
            ( Element::*mL )( Matrix< real > & aJ, Matrix< real > & aL );

            real * mRotationData = nullptr ;

            // edge directions, if this is a Nedelec element
            Bitset< 12 > mEdgeDirections ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // deprecated constructor
            Element( Group * aParent,
                     Field * aField,
                     mesh::Element * aElement );

//------------------------------------------------------------------------------

          Element(
                  Block         * aParent,
                  DofManager    * aDofManager,
                  mesh::Element * aElement );

//------------------------------------------------------------------------------

            /**
             * constructor for facet
             */
            Element(
                    SideSet       * aParent,
                    DofManager    * aDofManager,
                    mesh::Facet   * aFacet,
                    const SideSetDofLinkMode aSideSetDofLinkMode,
                    const id_t      aMasterBlockID,
                    const id_t      aSlaveBlockID );

//------------------------------------------------------------------------------

            /**
             * constructor for thin shells
             */
            Element(
                    SideSet       * aParent,
                    DofManager    * aDofManager,
                    mesh::Facet   * aFacet,
                    Cell< mesh::Facet * > & aLayers,
                    const id_t      aMasterBlockID,
                    const id_t      aSlaveBlockID ) ;

//------------------------------------------------------------------------------

            virtual ~Element();

//------------------------------------------------------------------------------

            /**
             * expose the underlying element on the mesh
             */
             mesh::Element *
             element();

//------------------------------------------------------------------------------

            /**
             * expose a specific dof
             */
            Dof * &
            dof( const index_t aIndex );

//------------------------------------------------------------------------------

            /**
             * expose the number of dofs
             */
            uint
            number_of_dofs() const ;

//------------------------------------------------------------------------------

            /**
             * flag all dofs that are connected to this element
             */
             void
             flag_dofs();

//------------------------------------------------------------------------------

            /**
             * unflag all dofs that are connected to this element
             */
             void
             unflag_dofs();

//------------------------------------------------------------------------------

            /**
             * get the precomputed shape function from parent block
             */
            const Matrix< real > &
            N( const uint aPointIndex ) const;

//------------------------------------------------------------------------------

            /**
             * get the precomputed shape function from parent block ( Geometry )
             */
            const Matrix< real > &
            G( const uint aPointIndex ) const;

//------------------------------------------------------------------------------

            // compute the transposed geometry Jacobian
            // rule: J * dNdX = dNdXi
            Matrix< real > &
            J( const uint aPointIndex );

//------------------------------------------------------------------------------

            // compute the K-Matrix for the geometry transformation
            // rule: L * d2NdX2 = d2NdXi2 - K * dNdXi
            Matrix< real > &
            K( const uint aPointIndex );

//------------------------------------------------------------------------------

            // compute the L-Matrix for the geometry transformation
            Matrix< real > &
            // rule: L * d2NdX2 = d2NdXi2 - K * dNdXi
            L( const uint aPointIndex );

//------------------------------------------------------------------------------

            /**
             * tell if this element has rotation properties
             */
             bool
             has_rotation() const;

//------------------------------------------------------------------------------

            void
            get_node_coors( Matrix< real > & aNodeCoords );

            /**
             * set the rotation axes
             */
            void
            set_rotation_data( const Matrix< real > & aRotationAxes, const Vector< real > & aRotationAngles );

            void
            get_rotation_data( Matrix< real > & aRotationAxes, Vector< real > & aRotationAngles );

//------------------------------------------------------------------------------

            /**
             * expose container for edge directions (obsolete)
             */
            // void
            // edge_directions( Vector< real > & aEdgeDirections ) const ;

//------------------------------------------------------------------------------

            /**
             * expose container for edge directions (obsolete)
             */
            void
            edge_directions( real * aEdgeDirections ) const ;

//------------------------------------------------------------------------------

            /**
             * expose container for edge directions
             */
            void
            edge_directions( Bitset< 12 > & aEdgeDirections ) const ;

//------------------------------------------------------------------------------

            /**
             * expose container for edge directions
             */
            bool
            edge_direction( const index_t aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the id of the underlying mesh element
             */
            id_t
            id() const ;

//------------------------------------------------------------------------------

            /**
             * return the master element, only for sidesets
             */
            Element *
            master();

//------------------------------------------------------------------------------

            /**
             * return the slave element, only for sidesets, returns null if
             * there is none
             */
            Element *
            slave();

//------------------------------------------------------------------------------

            mesh::Facet *
            facet() ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            link_dofs(
                    DofManager  * aDofManager,
                    const SideSetDofLinkMode aMode,
                    const id_t aSideSetID,
                    const id_t aMasterBlockID,
                    const id_t aSlaveBlockID );

//------------------------------------------------------------------------------

            /**
             * connect element with dofs of field (deprecated)
             */
            void
            link_dofs( Field * aField );

            /**
             * connect element with dofs of field
             */
            void
            link_dofs( DofManager * aDofManager, const id_t aBlockID );

//------------------------------------------------------------------------------

            void
            link_dofs_facet_only(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aMasterNodeDofTypes,
                    const Vector< index_t > & aSlaveNodeDofTypes,
                    const Vector< index_t > & aMasterEdgeDofTypes,
                    const Vector< index_t > & aSlaveEdgeDofTypes,
                    const Vector< index_t > & aMasterFaceDofTypes,
                    const Vector< index_t > & aSlaveFaceDofTypes,
                    const Vector< index_t > & aLambdaDofTypes );

//------------------------------------------------------------------------------

            void
            link_dofs_facet_and_master(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aMasterNodeDofTypes,
                    const Vector< index_t > & aSlaveNodeDofTypes,
                    const Vector< index_t > & aMasterEdgeDofTypes,
                    const Vector< index_t > & aSlaveEdgeDofTypes,
                    const Vector< index_t > & aMasterFaceDofTypes,
                    const Vector< index_t > & aSlaveFaceDofTypes,
                    const Vector< index_t > & aLambdaDofTypes  );

//------------------------------------------------------------------------------

            void
            link_dofs_facet_and_slave(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aMasterNodeDofTypes,
                    const Vector< index_t > & aSlaveNodeDofTypes,
                    const Vector< index_t > & aMasterEdgeDofTypes,
                    const Vector< index_t > & aSlaveEdgeDofTypes,
                    const Vector< index_t > & aMasterFaceDofTypes,
                    const Vector< index_t > & aSlaveFaceDofTypes,
                    const Vector< index_t > & aLambdaDofTypes );

//------------------------------------------------------------------------------

            void
            link_dofs_master_and_slave(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aMasterNodeDofTypes,
                    const Vector< index_t > & aSlaveNodeDofTypes,
                    const Vector< index_t > & aMasterEdgeDofTypes,
                    const Vector< index_t > & aSlaveEdgeDofTypes,
                    const Vector< index_t > & aMasterFaceDofTypes,
                    const Vector< index_t > & aSlaveFaceDofTypes,
                    const Vector< index_t > & aLambdaDofTypes );

//------------------------------------------------------------------------------

            void
            link_dofs_thin_shell(
                    DofManager  * aDofManager,
                    Cell< mesh::Facet * > & aLayers ,
                    const SideSetDofLinkMode aMode,
                    const id_t aMasterBlockID,
                    const id_t aSlaveBlockID );

//------------------------------------------------------------------------------

            void
            link_dofs_thin_shell_master_and_slave(
                    DofManager  * aDofManager,
                    Cell< mesh::Facet * > & aLayers ,
                    const Vector< index_t > & aMasterNodeDofTypes,
                    const Vector< index_t > & aSlaveNodeDofTypes,
                    const Vector< index_t > & aThinShellEdgeDofTypes,
                    const Vector< index_t > & aThinShellFaceDofTypes,
                    const Vector< index_t > & aLambdaDofTypes );

//------------------------------------------------------------------------------

            void
            link_dofs_thin_shell_facet_only(
                    DofManager  * aDofManager,
                    Cell< mesh::Facet * > & aLayers,
                    const Vector< index_t > & aThinShellNodeDofTypes,
                    const Vector< index_t > & aThinShellEdgeDofTypes,
                    const Vector< index_t > & aThinShellFaceDofTypes );

//------------------------------------------------------------------------------

            void
            grab_edge_directions_for_facet( Element * aReferenceElement );

//------------------------------------------------------------------------------

            void
            link_node_dofs(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aDofTypes,
                    mesh::Element * aReferenceElement,
                    const uint aNumberOfNodesPerElement,
                    index_t & aCount );

//------------------------------------------------------------------------------

            void
            link_edge_dofs(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aDofTypes,
                    mesh::Element * aReferenceElement,
                    const uint aNumberOfEdgesPerElement,
                    index_t & aCount );

//------------------------------------------------------------------------------

            void
            link_face_dofs(
                    DofManager  * aDofManager,
                    const Vector< index_t > & aDofTypes,
                    mesh::Element * aReferenceElement,
                    const uint aNumberOfFacesPerElement,
                    index_t & aCount );

//------------------------------------------------------------------------------

            void
            link_lambda_dofs(   DofManager  * aDofManager,
                                const Vector< index_t > & aDofTypes,
                                index_t & aCount );

//------------------------------------------------------------------------------

            void
            L1D( Matrix< real > & aJ, Matrix< real > & aL );

//------------------------------------------------------------------------------

            void
            L2D( Matrix< real > & aJ, Matrix< real > & aL );

//------------------------------------------------------------------------------

            void
            L3D( Matrix< real > & aJ, Matrix< real > & aL);

//------------------------------------------------------------------------------

            void
            compute_edge_directions();

//------------------------------------------------------------------------------

            void
            link_second_derivative_functions( const ElementType aElementType );

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        inline  Element *
        Element::master()
        {
            return mMaster ;
        }

//------------------------------------------------------------------------------

        inline  Element *
        Element::slave()
        {
            return mSlave ;
        }

//------------------------------------------------------------------------------

        inline  mesh::Facet *
        Element::facet()
        {
            return mFacet ;
        }

//------------------------------------------------------------------------------

        inline  mesh::Element *
        Element::element()
        {
            return mElement;
        }

//------------------------------------------------------------------------------

        inline Dof * &
        Element::dof( const index_t aIndex )
        {
            return mDOFs[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline uint
        Element::number_of_dofs() const
        {
            return mNumberOfDofs ;
        }

//------------------------------------------------------------------------------

        inline bool
        Element::has_rotation() const
        {
            return mRotationData != nullptr ;
        }


//------------------------------------------------------------------------------

        /*inline void
        Element::edge_directions( Vector< real > & aEdgeDirections ) const
        {
            for( uint e=0; e<mElement->number_of_edges(); ++e )
            {
                aEdgeDirections( e ) = mEdgeDirections.test( e ) ? 1.0 : -1.0 ;
            }
        }*/

//------------------------------------------------------------------------------

        inline void
        Element::edge_directions( real * aEdgeDirections ) const
        {
            for( uint e=0; e<mElement->number_of_edges(); ++e )
            {
                aEdgeDirections[ e ] = mEdgeDirections.test( e ) ? 1.0 : -1.0 ;
            }
        }

//------------------------------------------------------------------------------


        inline void
        Element::edge_directions( Bitset< 12 > & aEdgeDirections ) const
        {
            for( uint e=0; e<mElement->number_of_edges(); ++e )
            {
                if( mEdgeDirections.test( e ) )
                {
                    aEdgeDirections.set( e );
                }
                else
                {
                    aEdgeDirections.reset( e );
                }
            }
        }

//------------------------------------------------------------------------------

        inline bool
        Element::edge_direction( const index_t aIndex ) const
        {
            return mEdgeDirections.test( aIndex ) ;
        }

//------------------------------------------------------------------------------

        inline id_t
        Element::id() const
        {
            return mElement->id() ;
        }

//------------------------------------------------------------------------------

        inline void
        Element::get_node_coors( Matrix< real > & aNodeCoords )
        {
            uint tN = aNodeCoords.n_rows() ;
            uint tD = aNodeCoords.n_cols() ;

            for( uint i=0; i<tD; ++i )
            {
                for( uint k=0; k<tN; ++k )
                {
                    aNodeCoords( k, i ) = mElement->node( k )->x( i );
                }
            }
        }
    }

//------------------------------------------------------------------------------

}
#endif //BELFEM_CL_FEM_ELEMENT_HPP
