//
// Created by christian on 9/2/22.
//

#ifndef BELFEM_CL_MESH_ORIENTATION_CHECKER_HPP
#define BELFEM_CL_MESH_ORIENTATION_CHECKER_HPP

#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        class OrientationChecker
        {
            ElementType mElementType = ElementType::UNDEFINED ;

            //Cell< Node * > mNodes ;

            // node coordinates
            Vector< real > mX ;
            Vector< real > mY ;
            Vector< real > mZ ;

            // rotation matrix, needed for 3D
            Matrix< real > mT ;

            // normal vector
            Vector< real > mN ;

            // vector with node coordinates
            Vector< real > mP ;
            Vector< real > mQ ;

            // rotation axis
            Vector< real > mR ;

            uint mNumNodesPerElement       = BELFEM_UINT_MAX ;
            uint mNumCornerNodesPerElement = BELFEM_UINT_MAX ;

            // function that collects node coordinates
            void
            ( OrientationChecker::*mFunCollect )
                    (       Element       * aElement );


            // function that checks if the element orientation must be changed
            void
            ( OrientationChecker::*mFunProcess )
            (       Element       * aElement );

            // function that checks if the element orientation must be changed
            void
            ( OrientationChecker::*mFunFlip )
                    (       Element       * aElement );

            // map for flipping sidesets
            Map< uint, uint > mSideSetMap ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            OrientationChecker();

//------------------------------------------------------------------------------

            ~OrientationChecker() = default ;

//------------------------------------------------------------------------------

            void
            set_element_type( const ElementType aType );

//------------------------------------------------------------------------------

            void
            process_element( Element * aElement );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_node_coords_tri( Element * aElement );

//------------------------------------------------------------------------------

            void
            collect_node_coords_quad( Element * aElement );

//------------------------------------------------------------------------------

            void
            collect_node_coords_tet( Element * aElement );

//------------------------------------------------------------------------------

            void
            collect_node_coords_penta( Element * aElement );

//------------------------------------------------------------------------------

            void
            collect_node_coords_hex( Element * aElement );

//------------------------------------------------------------------------------

            void
            process_2d_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            process_tet_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            process_penta_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            process_hex_element( Element * aElement );

//-----------------------------------------------------------------------------

            void
            swap( Element * aElement, uint aNode1, uint aNode2 );

//------------------------------------------------------------------------------

            void
            flip_tri3( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tri6( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tri10( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tri15( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_quad4( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_quad9( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_quad16( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tet4( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tet10( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tet20( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_tet35( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_penta6( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_penta15( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_penta18( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_hex8( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_hex20( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_hex27( Element * aElement );

//------------------------------------------------------------------------------

            void
            flip_hex64( Element * aElement );

//------------------------------------------------------------------------------

            void
            populate_sideset_map( ElementType aType );

//------------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------

        inline void
        OrientationChecker::process_element( Element * aElement )
        {
            BELFEM_ASSERT( aElement->type() == mElementType,
                           "Element type does not match" );

            ( this->*mFunProcess )( aElement );
        }

//-----------------------------------------------------------------------------

        inline void
        OrientationChecker::swap( Element * aElement, uint aNode1, uint aNode2 )
        {
            Node * tSwap = aElement->node( aNode1 );
            aElement->insert_node( aElement->node( aNode2 ), aNode1 );
            aElement->insert_node( tSwap, aNode2 );
        }

//-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_ORIENTATION_CHECKER_HPP
