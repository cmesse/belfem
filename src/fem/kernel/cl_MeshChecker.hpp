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
#ifndef CL_MESHCHECKER_HPP
#define CL_MESHCHECKER_HPP
#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Pipette.hpp"
namespace belfem
{
    class MeshChecker
    {
        const proc_t mCommRank ;

        Mesh * mMesh = nullptr;

        mesh::Pipette * mPipette = nullptr;

        index_t mElementCount = 0 ;

        //! signature of the element-swapping routine selected by the constructor
        using SwapFunction = void ( MeshChecker::* )( mesh::Element * aElement );

        SwapFunction mFunSwap;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        MeshChecker( Mesh * aMesh );

//------------------------------------------------------------------------------

        ~MeshChecker();

//------------------------------------------------------------------------------

        index_t
        element_count() const ;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        link_to_block( mesh::Block * aBlock );

//------------------------------------------------------------------------------

        void
        process_block( mesh::Block * aBlock );

        void
        swap( mesh::Element * aElement );

        void
        swap( mesh::Element * aElement, uint aNode1, uint aNode2 );

        void
        swap_tri3( mesh::Element * aElement );

        void
        swap_tri6( mesh::Element * aElement );

        void
        swap_tri10( mesh::Element * aElement );

        void
        swap_tri15( mesh::Element * aElement );

        void
        swap_quad4( mesh::Element * aElement );

        void
        swap_quad9( mesh::Element * aElement );

        void
        swap_quad16( mesh::Element * aElement );

        void
        swap_tet4( mesh::Element * aElement );

        void
        swap_tet10( mesh::Element * aElement );

        void
        swap_tet20( mesh::Element * aElement );

        void
        swap_tet35( mesh::Element * aElement );

        void
        swap_penta6( mesh::Element * aElement );

        void
        swap_penta15( mesh::Element * aElement );

        void
        swap_penta18( mesh::Element * aElement );

        void
        swap_pyra5( mesh::Element * aElement );

        void
        swap_pyra14( mesh::Element * aElement );

        void
        swap_hex8( mesh::Element * aElement );

        void
        swap_hex8tb( mesh::Element * aElement );

        void
        swap_hex20( mesh::Element * aElement );

        void
        swap_hex27( mesh::Element * aElement );

        void
        swap_hex64( mesh::Element * aElement );

//------------------------------------------------------------------------------
    };
    
//------------------------------------------------------------------------------

    inline
    index_t
    MeshChecker::element_count() const
    {
        return mElementCount;
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap( mesh::Element * aElement )
    {
        (this->*mFunSwap)( aElement );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap( mesh::Element * aElement, uint aNode1, uint aNode2 )
    {
        mesh::Node * tSwap = aElement->node( aNode1 );
        aElement->insert_node( aElement->node( aNode2 ), aNode1 );
        aElement->insert_node( tSwap, aNode2 );
    }
    
//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tri3( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement, 1, 2 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tri6( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement, 1, 2 );
        this->swap( aElement, 3, 5 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tri10( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement, 1, 2 );
        this->swap( aElement, 3, 8 );
        this->swap( aElement, 4, 7 );
        this->swap( aElement, 5, 6 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tri15( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement,  1,  2 );
        this->swap( aElement,  3, 11 );
        this->swap( aElement,  4, 10 );
        this->swap( aElement,  5,  9 );
        this->swap( aElement,  6,  8 );
        this->swap( aElement, 13, 14 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_quad4( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement,  1,  3 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_quad9( mesh::Element * aElement )
    {
        // flip nodes
        this->swap( aElement,  1,  3 );
        this->swap( aElement,  4,  7 );
        this->swap( aElement,  5,  6 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_quad16( mesh::Element * aElement )
    {
        this->swap( aElement,  1,  3 );
        this->swap( aElement,  4, 11 );
        this->swap( aElement,  5, 10 );
        this->swap( aElement,  6,  9 );
        this->swap( aElement,  7,  8 );
        this->swap( aElement, 13, 15 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tet4( mesh::Element * aElement )
    {
        this->swap( aElement, 1, 2 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tet10( mesh::Element * aElement )
    {
        this->swap( aElement, 1, 2 );
        this->swap( aElement, 4, 6 );
        this->swap( aElement, 8, 9 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tet20( mesh::Element * aElement )
    {
        this->swap( aElement,  1,  2 );
        this->swap( aElement,  4,  9 );
        this->swap( aElement,  5,  8 );
        this->swap( aElement,  6,  7 );
        this->swap( aElement, 12, 14 );
        this->swap( aElement, 13, 15 );
        this->swap( aElement, 17, 18 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_tet35( mesh::Element * aElement )
    {
        this->swap( aElement,  1,  2 );
        this->swap( aElement,  4, 12 );
        this->swap( aElement,  5, 11 );
        this->swap( aElement,  6, 10 );
        this->swap( aElement,  7,  9 );
        this->swap( aElement, 16, 19 );
        this->swap( aElement, 17, 20 );
        this->swap( aElement, 18, 21 );
        this->swap( aElement, 23, 24 );
        this->swap( aElement, 25, 28 );
        this->swap( aElement, 26, 30 );
        this->swap( aElement, 27, 29 );
        this->swap( aElement, 32, 33 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_penta6( mesh::Element * aElement )
    {
        this->swap( aElement,  1,  2 );
        this->swap( aElement,  4,  5 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_penta15( mesh::Element * aElement )
    {
        this->swap( aElement,  1, 2 );
        this->swap( aElement,  4, 5 );
        this->swap( aElement,  6, 9 );
        this->swap( aElement, 10,11 );
        this->swap( aElement, 12,14 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_penta18( mesh::Element * aElement )
    {
        this->swap( aElement,  1, 2 );
        this->swap( aElement,  4, 5 );
        this->swap( aElement,  6, 9 );
        this->swap( aElement, 10,11 );
        this->swap( aElement, 12,14 );
        this->swap( aElement, 15,17 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_pyra5( mesh::Element * aElement )
    {
        this->swap( aElement,  0,  2 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_pyra14( mesh::Element * aElement )
    {
        this->swap( aElement,  0,  2 );
        this->swap( aElement,  5,  6 );
        this->swap( aElement,  7,  8 );
        this->swap( aElement,  9,  11 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_hex8( mesh::Element * aElement )
    {
        this->swap( aElement,  1, 3 );
        this->swap( aElement,  5, 7 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_hex8tb( mesh::Element * aElement )
    {
        // side connector elements carry edges that are tied to fixed node
        // slots; swapping the nodes would silently corrupt that pairing.
        // A negative volume means the factory built a left-handed element
        // and must be fixed there.
        BELFEM_ERROR( false,
            "side connector element %lu has negative volume",
            ( long unsigned int ) aElement->id() );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_hex20( mesh::Element * aElement )
    {
        this->swap( aElement,  1, 3 );
        this->swap( aElement,  5, 7 );
        this->swap( aElement,  8,11 );
        this->swap( aElement,  9,12 );
        this->swap( aElement, 10,15 );
        this->swap( aElement, 16,18 );
        this->swap( aElement, 19,17 );
    }

//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_hex27( mesh::Element * aElement )
    {
        this->swap( aElement,  1, 3 );
        this->swap( aElement,  5, 7 );
        this->swap( aElement,  8,11 );
        this->swap( aElement,  9,12 );
        this->swap( aElement, 10,15 );
        this->swap( aElement, 16,18 );
        this->swap( aElement, 19,17 );
        this->swap( aElement, 20,25 );
        this->swap( aElement, 22,23 );
    }
//------------------------------------------------------------------------------

    inline void
    MeshChecker::swap_hex64( mesh::Element * aElement )
    {
        this->swap( aElement, 1, 3 );
        this->swap( aElement, 5, 7 );
        this->swap( aElement, 8, 10 );
        this->swap( aElement, 9, 11 );
        this->swap( aElement, 14, 19 );
        this->swap( aElement, 15, 18 );
        this->swap( aElement, 16, 22 );
        this->swap( aElement, 17, 23 );
        this->swap( aElement, 24, 26 );
        this->swap( aElement, 25, 27 );
        this->swap( aElement, 28, 31 );
        this->swap( aElement, 29, 30 );
        this->swap( aElement, 33, 35 );
        this->swap( aElement, 36, 40 );
        this->swap( aElement, 37, 43 );
        this->swap( aElement, 38, 42 );
        this->swap( aElement, 39, 41 );
        this->swap( aElement, 44, 49 );
        this->swap( aElement, 45, 48 );
        this->swap( aElement, 46, 51 );
        this->swap( aElement, 47, 50 );
        this->swap( aElement, 53, 55 );
        this->swap( aElement, 57, 59 );
        this->swap( aElement, 61, 63 );
    }


}
#endif //CL_MESHCHECKER_HPP
