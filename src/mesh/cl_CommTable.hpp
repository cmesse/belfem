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

#ifndef BELFEM_CL_COMMTABLE_HPP
#define BELFEM_CL_COMMTABLE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace mesh
    {
        class Distributor;

        class CommTable
        {
            // contains the indices
            Cell< index_t > mNodes ;
            Cell< index_t > mEdges ;
            Cell< index_t > mFaces ;
            Cell< index_t > mElements ;
            Cell< index_t > mFacets ;
            Cell< index_t > mVertices ;
            Cell< index_t > mControlPoints ;

            // give distributor access to writable functions
            friend Distributor ;

        public:

            CommTable() = default ;

            ~CommTable() = default;

            const Cell< index_t > & nodes() const;
            const Cell< index_t > & edges() const;
            const Cell< index_t > & faces() const;
            const Cell< index_t > & elements() const;
            const Cell< index_t > & facets() const;
            const Cell< index_t > & vertices() const;
            const Cell< index_t > & control_points() const;

            index_t number_of_nodes() const;
            index_t number_of_edges() const;
            index_t number_of_faces() const;
            index_t number_of_elements() const;
            index_t number_of_facets() const;
            index_t number_of_vertices() const;
            index_t number_of_control_points() const;

        protected:

            Cell< index_t > & nodes();
            Cell< index_t > & edges();
            Cell< index_t > & faces();
            Cell< index_t > & elements();
            Cell< index_t > & facets();
            Cell< index_t > & vertices();
            Cell< index_t > & control_points();
        };

//-----------------------------------------------------------------------

        inline const Cell< index_t > &
        CommTable::nodes() const
        {
            return mNodes ;
        }

        inline const Cell< index_t > &
        CommTable::edges() const
        {
            return mEdges ;
        }

        inline const Cell< index_t > &
        CommTable::faces() const
        {
            return mFaces ;
        }

        inline const Cell< index_t > &
        CommTable::elements() const
        {
            return mElements ;
        }
        inline const Cell< index_t > &
        CommTable::facets() const
        {
            return mFacets ;
        }
        inline const Cell< index_t > &
        CommTable::vertices() const
        {
            return mVertices ;
        }

        inline Cell< index_t > &
        CommTable::nodes()
        {
            return mNodes ;
        }

        inline Cell< index_t > &
        CommTable::edges()
        {
            return mEdges ;
        }

        inline Cell< index_t > &
        CommTable::faces()
        {
            return mFaces ;
        }

        inline Cell< index_t > &
        CommTable::elements()
        {
            return mElements ;
        }
        inline Cell< index_t > &
        CommTable::facets()
        {
            return mFacets ;
        }
        inline Cell< index_t > &
        CommTable::vertices()
        {
            return mVertices ;
        }

        inline index_t
        CommTable::number_of_nodes() const
        {
            return mNodes.size();
        }

        inline index_t
        CommTable::number_of_edges() const
        {
            return mEdges.size();
        }

        inline index_t
        CommTable::number_of_faces() const
        {
            return mFaces.size();
        }

        inline index_t
        CommTable::number_of_elements() const
        {
            return mElements.size();
        }

        inline index_t
        CommTable::number_of_facets() const
        {
            return mFacets.size();
        }

        inline index_t
        CommTable::number_of_vertices() const
        {
            return mVertices.size();
        }

        inline index_t
        CommTable::number_of_control_points() const
        {
            return mControlPoints.size();
        }

        inline const
        Cell< index_t > & CommTable::control_points() const
        {
            return mControlPoints ;
        }

        inline
        Cell< index_t > & CommTable::control_points()
        {
            return mControlPoints ;
        }
    }
}
#endif //BELFEM_CL_COMMTABLE_HPP