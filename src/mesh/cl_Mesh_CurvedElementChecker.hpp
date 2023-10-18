//
// Created by Christian Messe on 05.01.22.
//

#ifndef BELFEM_CL_MESH_CURVEDELEMENTCHECKER_HPP
#define BELFEM_CL_MESH_CURVEDELEMENTCHECKER_HPP

#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        /**
         * this temporary class is created by mesh->flag_curved_elements().
         * it scans all elements and checks if the element is curved or not
         */
        class CurvedElementChecker
        {
            const proc_t mMyRank ;

            const uint mMeshDimension ;

            // blocks on mesh
            Cell< Block * > & mBlocks ;

            // sidesets on mesh
            Cell< SideSet * > & mSideSets ;

            // epsilon environment for mesh
            const real mTwoMeshEpsilon   = 2e-9 ;
            const real mFourMeshEpsilon  = mTwoMeshEpsilon  + mTwoMeshEpsilon ;
            const real mEightMeshEpsilon = mFourMeshEpsilon + mFourMeshEpsilon ;

            // for quad4
            Cell< Node * > mNodes ;
            real mX[ 4 ];
            real mY[ 4 ];
            real mZ[ 4 ];


            // function pointer to checker
            bool
            ( CurvedElementChecker::*mFunCheck )( Element * aElement );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            CurvedElementChecker(
                    const uint aMeshDimension,
                    Cell< Block * >   & aBlocks,
                    Cell< SideSet * > & aSideSets );

            ~CurvedElementChecker() = default ;

            index_t
            flag_curved_elements();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            link_check_function( const ElementType aElementType );

//------------------------------------------------------------------------------

            bool
            check_linear( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_tri6_2d( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_quad8_2d( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_quad9_2d( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_tri6( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_quad4( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_quad4_face( Element * aElement, const uint aFaceIndex );

//------------------------------------------------------------------------------

            bool
            check_quad8( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_quad9( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_tet10( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_penta6( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_penta15( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_penta18( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_hex8( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_hex20( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_hex27( Element * aElement );

//------------------------------------------------------------------------------

            bool
            check_midpoint( const real aA,
                            const real aB,
                            const real aC ) const ;

//------------------------------------------------------------------------------

            bool
            check_midpoint( const real aA,
                            const real aB,
                            const real aC,
                            const real aD,
                            const real aE ) const ;

//------------------------------------------------------------------------------

            bool
            check_midpoint( const real aA,
                            const real aB,
                            const real aC,
                            const real aD,
                            const real aE,
                            const real aF,
                            const real aG,
                            const real aH,
                            const real aI ) const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline bool
        CurvedElementChecker::check_linear( Element * aElement )
        {
            // do nothing, linear elements are never curved
            return false ;
        }

//------------------------------------------------------------------------------
        inline bool
        CurvedElementChecker::check_midpoint( const real aA,
                                              const real aB,
                                              const real aC ) const
        {
            return std::abs( aA + aB - aC - aC ) > mTwoMeshEpsilon ;
        }

//------------------------------------------------------------------------------

        inline bool
        CurvedElementChecker::check_midpoint( const real aA,
                                              const real aB,
                                              const real aC,
                                              const real aD,
                                              const real aE ) const
        {
            return std::abs( aA + aB + aC + aD - aE - aE - aE - aE ) > mFourMeshEpsilon ;
        }

//------------------------------------------------------------------------------

        inline bool
        CurvedElementChecker::check_midpoint( const real aA,
                                              const real aB,
                                              const real aC,
                                              const real aD,
                                              const real aE,
                                              const real aF,
                                              const real aG,
                                              const real aH,
                                              const real aI ) const
        {
            return std::abs( aA + aB + aC + aD + aE + aF + aG + aH
                -aI - aI - aI - aI -aI - aI - aI - aI ) > mEightMeshEpsilon ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_CURVEDELEMENTCHECKER_HPP
