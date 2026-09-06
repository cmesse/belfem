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

#ifndef CL_CUTSET_HPP
#define CL_CUTSET_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        class CutSet
        {
            Mesh * mMesh ;

            //! list of node candidates relevant for the duplication
            Cell< Node * > & mNodeOriginals ;

            //! the bitset containing the cut pattern
            DynamicBitset * mBitset ;

            //! the bitset containing the node switches
            DynamicBitset * mNodeBitset ;

            //! the map containing the node duplicates
            Map< id_t, Node * > mNodeDuplicates ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            CutSet(
                Mesh * aMesh,
                Cell< Node * >          & aNodeOriginals,
                const string & aHexString,
                const index_t  aNumberOfCuts );

            ~CutSet();

            DynamicBitset *
            node_bitset();

            bool
            test( const Node * aNode ) const ;

//-----------------------------------------------------------------------------

            void
            create_duplicates( id_t & aMaxNodeID, Cell< Node * > & aAbstractNodes );

//-----------------------------------------------------------------------------

            Node *
            duplicate( Node * aNode );

//-----------------------------------------------------------------------------

            Map< id_t, Node * > &
            duplicate_map() ;

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------

        inline bool
        CutSet::test( const Node * aNode ) const
        {
            return mBitset->test( aNode->index() );
        }

//-----------------------------------------------------------------------------

        inline DynamicBitset *
        CutSet::node_bitset()
        {
            return mNodeBitset;
        }

//-----------------------------------------------------------------------------

        inline Node *
        CutSet::duplicate( Node * aNode )
        {
            return mNodeDuplicates[ aNode->id() ];
        }

//-----------------------------------------------------------------------------

        inline Map< id_t, Node * > &
        CutSet::duplicate_map()
        {
            return mNodeDuplicates;
        }

//-----------------------------------------------------------------------------
    }
}
#endif //CL_CUTSET_HPP
