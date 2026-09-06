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

#ifndef CL_CUT_DATA_HPP
#define CL_CUT_DATA_HPP


#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Mesh.hpp"
#include "cl_Cohomology.hpp"

namespace belfem
{
    namespace mesh
    {
        /**
         * @brief Per-cut topology and metadata.
         *
         * @ingroup grp_homology
         * @see @ref homology_homology_usage_guide
         */
        class CutData
        {
            Mesh             * mMesh ;
            Cochain          * mCochain ;

            const index_t      mIndex ;
            const id_t         mID ;

            const uint         mNumberOfCuts ;

            Cell< Edge * >     mCohomologyEdges ;

            // needed for weight computation, for each edge
            DynamicBitset *    mCohomologyPlus = nullptr ;
            DynamicBitset *    mCohomologyMinus = nullptr ;

            Cell< Element * >  mCutElements ; // elements that contribute to this cut
            Vector< int >      mCutCases ;
            Map< id_t, int >   mCutCaseMap ;
            Cell< Edge * >     mThinCutEdges ; // edges of the thin cut (2D): define the nodes to duplicate and the emitted cut sideset
            Cell< Face * >     mThinCutFaces ; // faces of the thin cut (3D): define the nodes to duplicate and the emitted cut sideset

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            CutData(
                      Mesh           * aMesh,
                      Cohomology     * aCohomology,
                const index_t          aCohomologyIndex );

            ~CutData();

            void
            add_thin_cut_sidesets_to_mesh( id_t & aMaxSideSetID, id_t & aMaxElementID );

            Cell< Element * > &
            elements();

            Cell< Edge * > &
            cohomology_edges();

            int
            weight( const Edge * aEdge );

            int
            cut_case( const index_t aLocalElementIndex );

            void
            collect_coefficients( const index_t aNumEdges );

            void
            collect_elements( const Vector< id_t > & aNonPhiDomains );

            void
            collect_edges();

            Mesh *
            create_debug_mesh();

//-----------------------------------------------------------------------------

            void
            flag_nodes();

//-----------------------------------------------------------------------------

            Cell< Edge * > &
            thin_cut_edges();

//-----------------------------------------------------------------------------

            Cell< Face * > &
            thin_cut_faces();

//-----------------------------------------------------------------------------

            index_t
            index() const ;

//-----------------------------------------------------------------------------

            void
            flag_all_elements();

//-----------------------------------------------------------------------------

            void
            unflag_all_elements();

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------


            //void
            //flag_elements_on_cuts( const Vector< id_t > & aNonPhiDomains );

            int
            determine_cut_case_2d( Element * aElement );

            int
            determine_cut_case_3d( Element * aElement );

        };

//-----------------------------------------------------------------------------

        inline int
        CutData::weight( const Edge * aEdge )
        {
            index_t tIndex = aEdge->index();
            if ( mCohomologyPlus->test( tIndex ) )
            {
                return 1;
            }
            if ( mCohomologyMinus->test( tIndex ) )
            {
                return -1;
            }
            return 0;
        }

//-----------------------------------------------------------------------------

        inline int
        CutData::cut_case( const index_t aLocalElementIndex )
        {
            return mCutCases( aLocalElementIndex );
        }

//-----------------------------------------------------------------------------

        inline
        Cell< Element * > &
        CutData::elements()
        {
            return mCutElements;
        }

//-----------------------------------------------------------------------------

        inline
        Cell< Edge * > &
        CutData::cohomology_edges()
        {
            return mCohomologyEdges;
        }

//-----------------------------------------------------------------------------

        inline
        Cell< Edge * > &
        CutData::thin_cut_edges()
        {
            return mThinCutEdges ;
        }

//-----------------------------------------------------------------------------

        inline
        Cell< Face * > &
        CutData::thin_cut_faces()
        {
            return mThinCutFaces ;
        }

//-----------------------------------------------------------------------------

        inline index_t
        CutData::index() const
        {
            return mIndex;
        }

//-----------------------------------------------------------------------------
    }
}


#endif // CL_CUT_DATA_HPP
