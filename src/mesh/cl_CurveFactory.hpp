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

#ifndef CL_CURVEFACTORY_HPP
#define CL_CURVEFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Curve.hpp"
#include "cl_Mesh.hpp"
#include "cl_Segment.hpp"
#include "cl_SideSet.hpp"

namespace belfem
{
    namespace mesh
    {
        class CurveFactory
        {
            const proc_t mCommRank ;
            Mesh * mMesh ;
            id_t mMaxElementID ;
            id_t mMaxCurveID ;

            struct ProtoCurve
            {
                Cell< Node * > mNodes ;
                Cell< Edge * > mEdges ;
            };

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            CurveFactory( Mesh * aMesh ) ;

            ~CurveFactory() = default ;

            /*
             * This function just checks if two sidesets intersect
             */
            bool
            intersection_exists( const id_t aSideSetA, const id_t aSideSetB ) ;

            /*
             * This function creates a curve from an intersection between
             * two surfaces. The intersection must be continuous.
             */
            Curve *
            intersect( const id_t aSideSetA, const id_t aSideSetB, id_t aID ) ;

            /*
            * This function creates a curve from a 1D sideset
            */
            Curve *
            from_2d_sidesets( const Vector< id_t > & aSideSets, id_t aID ) ;

            Cell< Curve * >
            thin_shell_side_curves( const Vector< id_t > & aThinShellSideSets, const Vector< id_t > & aDomainBoundaries );

            id_t &
            max_segment_id() ;

            void
            collect_nodes( Curve * aCurve );

            void
            compute_coordinates( Curve * aCurve );

            /*
             * The three sorters below are one pipeline and must be called in
             * this order: sort_end_nodes, orient_segments, sort_segments.
             * On entry to sort_end_nodes, every node of aEndNodes carries its
             * own slot of that cell as index() and is unflagged, and every
             * segment of the curve carries its own slot of the segment cell
             * as index() and is unflagged. sort_end_nodes leaves the visited
             * nodes flagged, and its "not one connected chain" refusal also
             * fires on a chain whose node flags were dirty on entry.
             * sort_end_nodes gives every end node its walk index 0 ... n-1
             * and sets the closed flag; orient_segments directs each segment
             * along the walk, so that on a closed curve the segment joining
             * the last node to the first runs n-1 -> 0; sort_segments then
             * orders the segments by the walk index of their node 0.
             */
            void
            sort_end_nodes( Curve * aCurve, Matrix< index_t > & aAdjacency, Cell< Node * > & aEndNodes );

            void
            orient_segments( const bool aClosed, Cell< Node * > & aEndNodes, Cell< Segment * > & aSegments );

            void
            sort_segments( Cell< Segment * > & aSegments );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_end_nodes_from_intersection(
                SideSet * aSidesetA,
                SideSet * aSideSetB,
                Cell< Node * > & aEndNodes );

            /**
             * one pair per facet edge whose two corners are end nodes, ordered
             * by the key max( slot ) * n + min( slot ) over the end-node slots;
             * if aQuadratic, midpoint k is node 2 of the edge of pair k. If
             * aSharedWith is given, only the edges that also belong to one of
             * its facets count, and their midpoints must agree. Expects every
             * end node to carry its slot in aEndNodes as index() and every
             * other facet corner gNoIndex: a line facet has no other corner,
             * so from_2d_sidesets() meets that as is, while intersect() resets
             * all node indices first. On a refusal, the mesh node indices are
             * restored first
             */
            void
            create_pairs_from_facets(
                const Cell< Facet * > & aFacets,
                Cell< Node * >        & aEndNodes,
                const bool              aQuadratic,
                Cell< std::pair< Node *, Node * > > & aPairs,
                Cell< Node * > & aMidNodes,
                const Cell< Facet * > * aSharedWith = nullptr );

            /**
             * node 2 of the edge of aFacet's element whose corners are aA and
             * aB, or nullptr if there is no such edge or it has no midpoint.
             * aWork is scratch
             */
            Node *
            midnode_of_facet_edge(
                Facet * aFacet,
                const Node * aA,
                const Node * aB,
                Cell< Node * > & aWork );

            void
            create_segments( Curve * aCurve, Cell< std::pair< Node *, Node * > > & aPairs, Cell< Node * > & aMidNodes );

            Node *
            next( const Matrix< index_t > & aAdjacency, Cell< Segment * > & aSegments, const Node * aNode );


            void
            collect_nodes( const Vector< id_t > & aSideSets, Cell< Node * > & aNodes );

            void
            create_edges_on_sidesets(
                const Vector< id_t > & aSideSets ,
                      Cell< Node * > & aNodes,
                      Cell< Edge * >       & aEdges,
                      Map< key_t, Edge * > & aMap );

            void
            connect_facets_to_edges(
                const Vector< id_t > & aSideSets,
                const key_t aNumNodes,
                Cell< Edge * >       & aEdges,
                Map< key_t, Edge * > & aMap  );

            void
            select_edges(
                const Vector< id_t > & aThinShellSideSets,
                const Vector< id_t > & aDomainBoundaries,
                      Cell< Edge * > & aAllEdges,
                      Cell< Edge * > & aSelectedEdges );

            void
            select_node_subset( Cell< Edge * > & aEdges, Cell< Node * > & aNodes );

            void
            connect_edges_to_nodes( Cell< Edge * > & aEdges, Cell< Node * > & aNodes );

            index_t
            identify_subchains( Cell< Node * > & aNodes );

            void
            create_protocurves(  Cell< Node * > & aNodes, Cell< Edge * > & aEdges, Cell< ProtoCurve * > & aProtoCurves );



        };

        inline id_t &
        CurveFactory::max_segment_id()
        {
            return mMaxElementID ;
        }

    }
}

#endif //CL_CURVEFACTORY_HPP
