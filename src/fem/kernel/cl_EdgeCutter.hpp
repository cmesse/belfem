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

#ifndef BELFEM_CL_EDGECUTTER_HPP
#define BELFEM_CL_EDGECUTTER_HPP
#include "cl_Mesh.hpp"
#include "cl_IF_ElementMapper.hpp"
#include "cl_Curve.hpp"
namespace belfem
{
    namespace mesh
    {
        class EdgeCutter
        {
            Mesh * mMesh;
            fem::ElementMapper * mElementMapper ;

            Cell< SideSet * > mSideSets ;

            DynamicBitset * mFacetBitset = nullptr ;
            Cell< Facet * > mFacets ;
            Vector< index_t > mFacetsPerLayerBegin ;
            Vector< index_t > mFacetsPerLayerEnd ;

            real mConnectorWitdh = 7.5e-4 ;

            // work vectors
            Vector< real > mP ;
            Vector< real > mQ ;
            Vector< real > mU ;
            Vector< real > mV ;
            Vector< real > mW ;
        public:

            EdgeCutter( Mesh * aMesh );

            ~EdgeCutter();

            void
            select_sidesets( const Vector< id_t > & aSideSets );

            void
            process_curve( Curve * aCurve );

        private:

            void
            facet_bfs( Curve * aCurve );

            void
            compute_node_vectors( Curve * aCurve,
                Matrix< real > & aNormals,
                Matrix< real > & aTangents,
                Matrix< real > & aBinomials,
                Vector< real > & aDistances );

            void
            compute_temporary_nodes( Curve * aCurve, Matrix< real > & aBinomials, Cell< Node * > & aNodes );

            void
            project_temporary_nodes(
                Cell< Facet * > & aFacets,
                Cell< Node * > & aNodes );

            real
            determine_sign( Curve * aCurve, Matrix< real > & aBinomials );

            void
            save_testmesh( Cell< Node * > & aNodes );

            bool
            inside_bounding_box( const Vector< real > & aPoint, const Vector< real > & aXmin, const Vector< real > & aXmax ) const;

            void
            project_to_plane( Facet * aFacet ) const;
        };
    }
}
#endif //BELFEM_CL_EDGECUTTER_HPP
