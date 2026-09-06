/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MESH_PERIODICITYFACTORY_HPP
#define BELFEM_CL_MESH_PERIODICITYFACTORY_HPP

#include "cl_Mesh.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "cl_ProtoMesh.hpp"

namespace belfem
{
    namespace mesh
    {
        /**
         * Creates a Periodicity object by matching nodes on two periodic
         * boundary planes. The planes are defined by three .geo Point IDs
         * each. Nodes on sidesets lying on each plane are detected
         * automatically. Master-to-slave correspondence is established
         * by projecting both sides into a shared in-plane coordinate
         * system: facets are paired through their centroids with a 2D
         * k-d tree nearest-neighbor search, and nodes are then matched
         * within each facet pair.
         */
        class PeriodicityFactory
        {
            Mesh * mMesh;
            ProtoMesh * mProtoMesh ;

            // plane equation: dot(X, N) = D, where N is row 2 of T
            Cell< Node * > mSourcePlane ;
            Cell< Node * > mTargetPlane ;

                    real   mSourceDistance = BELFEM_QUIET_NAN ;
            Matrix< real > mSourceTransform ;
                    real   mTargetDistance = BELFEM_QUIET_NAN ;
            Matrix< real > mTargetTransform ;

            Vector< real > mSourceHesse ;
            Vector< real > mTargetHesse ;

            DynamicBitset * mNodeBitset = nullptr ;
            DynamicBitset * mFaceBitset = nullptr ;
            DynamicBitset * mEdgeBitset = nullptr ;

        public:

            PeriodicityFactory( Mesh * aMesh, ProtoMesh * aProtoMesh=nullptr );

            ~PeriodicityFactory();

            void
            set_master_plane( const Vector< id_t > & aPointIDs );

            void
            set_slave_plane( const Vector< id_t > & aPointIDs );

            void
            set_master_plane( Cell< Node * > & aNodes );

            void
            set_slave_plane( Cell< Node * > & aNodes );

            /**
             * Define the master plane from three .geo Point IDs.
             * Points A, B, C must not be collinear.
             */
            void
            set_master_plane( const id_t A, const id_t B, const id_t C );

            /**
             * Define the slave plane from three .geo Point IDs.
             * The point mapping (e.g. 11->11, 14->15, 16->16) must
             * produce the same in-plane coordinate system as the master.
             */
            void
            set_slave_plane( const id_t A, const id_t B, const id_t C );

            Periodicity *
            create_periodicity() ;

            proto::PeriodicityData *
            to_proto( Periodicity * aPeriodicity );

            Periodicity *
            from_proto( proto::PeriodicityData * aPeriodicity, const bool aCrosslinkEntities );

            /**
             * Match nodes, then derive edge and face pairs.
             * The resulting Periodicity object is owned by the mesh.
             */
            void
            update_periodicity( Periodicity * aPeriodicity );

            void
            tag_periodic_sidesets();

        private:

            void
            map_facets(
                Cell< Facet * >   & aSourceFacets,
                Cell< Facet * >   & aTargetFacets );

            void
            collect_nodes(
                Cell< Facet * >   & aFacets,
                Cell< Node * >    & aNodes );

            void
            collect_edges(
                Cell< Facet * >   & aFacets,
                Cell< Edge * >    & aEdges );

            bool
            match_nodes(
                    Cell< Facet * >   & aSourceFacets,
                    Cell< Facet * >   & aTargetFacets,
                    Cell< Node * >    & aSourceNodes,
                    Cell< Node * >    & aTargetNodes );
            bool
            match_edges(
                Cell< Facet * > & aSourceFacets,
                Cell< Facet * > & aTargetFacets,
                Cell< Node * >  & aSourceNodes,
                Cell< Edge * >  & aSourceEdges,
                Cell< Node * >  & aTargetNodes,
                Cell< Edge * >  & aTargetEdges  );

            bool
            match_facets_and_faces(
                Cell< Node * >  & aSourceNodes,
                Cell< Facet * > & aSourceFacets,
                Cell< Face * >  & aSourceFaces,
                Cell< Node * >  & aTargetNodes,
                Cell< Facet * > & aTargetFacets,
                Cell< Face * >  & aTargetFaces );

            void
            create_edge_map( Cell< Node * > & aNodes, Cell< Edge * > & aEdges, Map< key_t, Edge * > & aEdgeMap );

            void
            create_facet_map( Cell< Node * > & aNodes, Cell< Facet * > & aFaces, Map< key128_t, Facet * > & aFacetMap );

            // builds orthonormal basis [P; Q; N] and plane offset D from 3 nodes
            void
            compute_transformation_matrix(
                const Node * A,
                const Node * B,
                const Node * C,
                      real & D,
                Matrix< real > & T,
                Vector< real > & H );

            // recursive k-d tree nearest-neighbor search (Euclidean distances)
            Node *
            find_closest_node( Node   * aNode,
                               real     aX,
                               real     aY,
                               bool     aFlip,
                               Node   * aBest,
                               real   & aBestDist );

            // collects the sidesets whose nodes all sit on a given plane
            void
            select_sidesets(
                const Matrix< real > & aTransform,
                const real aDistance,
                Cell< SideSet * > & aSideSets );

            void
            create_temporary_nodes(
                const Matrix< real > & aTransform,
                const real aDistance,
                const Cell< SideSet * > & aSideSets,
                Cell< Node * > & aNodes );

            // builds 2D k-d tree from transformed node positions
            Node *
            create_kdtree( Cell< Node * > & aNodes );

            // recursive median-based k-d tree construction, alternating x/y splits
            Node *
            kdsort( Cell< Node * > & aNodes, index_t aStart, index_t aEnd, const bool aFlip );

            // reassigns slave faces from master to slave orientation
            void
            fix_face_slaves( Periodicity * aPeriodicity );

            void
            create_hesse( Cell< Node * > & aNodes , Vector< real > & aHesse );

            void
            reset_bitsets( const index_t aNumNodes, const index_t aNumEdges, const index_t aNumFaces );

            void
            crosslink( Periodicity * aPeriodicity );

            template < typename T >
            void
            crosslink( Cell< T * > & aMaster, Cell< T * > & aSlave )
            {

                BELFEM_ASSERT( aMaster.size() == aSlave.size() , "Periodic entity sizes do not match");

                index_t n = aMaster.size();

                for ( index_t k=0; k<n; ++k )
                {
                    T * A = aMaster( k );
                    T * B = aSlave( k );

                    if ( A != B )
                    {
                        A->set_periodic( B );
                        B->set_periodic( A );
                    }
                }
            }
        };
    }
}
#endif //BELFEM_CL_MESH_PeriodicityFactory_HPP