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

#ifndef CL_MESH_CONNECTIVITYCALCULATOR_HPP
#define CL_MESH_CONNECTIVITYCALCULATOR_HPP
#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "op_Node_Index.hpp"
namespace belfem
{
    namespace mesh
    {
        /**
         * @class ConnectivityCalculator
         * @brief Computes element-to-element and node-to-node connectivity for finite element meshes
         *
         * This class computes entity-to-entity connectivities (node/edge/face/facet/element/
         * control point, plus thin-shell element neighbours) on the mesh it is given, in serial
         * on the calling rank.
         *
         * Key features:
         * - Uses efficient vector+sort+unique approach instead of bitsets for sparse connectivity
         * - Handles thin-shell ghost facets when linking edges, faces and elements across layers
         */
        class ConnectivityCalculator
        {
            //! Current MPI process rank
            const proc_t mCommRank ;

            //! Total number of MPI processes
            const proc_t mCommSize ;

            //! Pointer to the mesh being processed
            Mesh * mMesh ;


            Cell< Node * >    & mNodes ;
            Cell< Edge * >    & mEdges  ;
            Cell< Face * >    & mFaces;
            Cell< Element * > & mElements ;
            Cell< Facet * >   & mFacets ;
            Cell< ControlPoint * > & mControlPoints ;

            //! IDs of mesh blocks selected for connectivity computation
            //! If empty, all blocks are processed
            Vector< id_t >    mSelectedBlockIDs;

            //! Internal flag checking if we only calculate a partial mesh
            bool mPartialOnly = false ;

            Vector< proc_t > mElementOwners ;
            Vector< proc_t > mNodeOwners ;

            key128_t ( ConnectivityCalculator:: * mKey )( Cell< Node * > & aNodes ) const ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

            /**
             * @brief Constructor
             * @param aMesh Pointer to the mesh to compute connectivity for
             */
            ConnectivityCalculator( Mesh * aMesh );

            /**
             * @brief Destructor (nothing to release; all containers are owned by the mesh)
             */
            ~ConnectivityCalculator();

            /**
            * @brief Compute node-to-element connectivity
            */
            void
            connect_nodes_to_elements() ;

            void
            connect_facets_to_elements();

            void
            connect_nodes_to_facets();

            void
            connect_facets_to_facets();

            void
            connect_faces_to_edges_and_edges_to_faces();

            void
            connect_faces_to_faces();

            void
            connect_nodes_to_edges();

            void
            connect_edges_to_elements( Cell< Element * > & aElements );

            void
            connect_edges_to_ghost_facets();

            void
            connect_faces_to_ghost_facets();

            void
            connect_edges_to_edges();

            /**
             * @brief Compute node-to-node connectivity
             *
             * For each node, finds all other nodes that share at least one element
             * with it (all elements attached to the node are considered).
             */
            void
            connect_nodes_to_nodes();

            /**
             * @brief Compute element-to-element connectivity
             *
             * For each element, finds all other elements that share a facet with it
             * (an edge in 2D, a face in 3D). When thin shells exist, the pass is repeated
             * on the thin-shell blocks alone so that layer elements are linked to their
             * in-layer neighbours.
             */
            void
            connect_elements_to_elements();

            void
            connect_thin_shells_to_thin_shells();

            void
            connect_control_points_to_elements() ;

            void
            connect_control_points_to_control_points() ;

            /**
             * @brief stores the computed connectivities for debugging
             *
             * @param aPath Path of the HDF5 file where the data is stored
             */
            void
            save_connectivity_data( const string & aPath );

        private:

            void
            connect_elements_to_elements_sub( Cell< Element * > & aElements );

            key128_t
            facet_key2( Cell< Node * > & aNodes ) const ;

            key128_t
            facet_key3( Cell< Node * > & aNodes ) const ;

        };

        inline key128_t ConnectivityCalculator::facet_key2( Cell< Node * > & aNodes ) const
        {
            sort( aNodes, opNodeIndex );
            key128_t tNumNodes = mMesh->number_of_nodes();
            return aNodes( 1 )->index() * tNumNodes + aNodes( 0 )->index();
        }

        inline key128_t ConnectivityCalculator::facet_key3( Cell< Node * > & aNodes ) const
        {
            sort( aNodes, opNodeIndex );
            key128_t tNumNodes = mMesh->number_of_nodes();
            return ( aNodes( 2 )->index() * tNumNodes + aNodes( 1 )->index() )
                * tNumNodes+ aNodes( 0 )->index();
        }

    }
}
#endif //CL_MESH_CONNECTIVITYCALCULATOR_HPP
