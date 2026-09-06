//
// Created by christian on 8/20/25.
//

#ifndef BELFEM_CL_EF_PENTA6TS_HPP
#define BELFEM_CL_EF_PENTA6TS_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        // the friend below takes a Calculator; today the declaration arrives
        // through cl_FEM_Element.hpp's include chain, but this header must
        // not depend on an include order it does not control
        class Calculator ;

        namespace maxwell
        {
            const Matrix< real > &
            side_connector_edge_function( Calculator * aCalc, const real xi, const real eta, const real zeta );
        }

        class EF_PENTA6TS : public EdgeFunction
        {
            // edge directions
            real mS[ 6 ];

            // derivatives
            real mNablaXi[ 3 ];
            real mNablaEta[ 3 ];
            real mNablaZeta[ 3 ];

            // constant in-plane gradient tensor for G(), layout i + 3*j
            // ( W = nXi (x) nEta - nEta (x) nXi; identical for all three
            //   edge pairs because nZeta = -nXi - nEta ); filled in link()
            real mGradW[ 9 ];

            // normal
            real mN[ 3 ];
            real mDetGram ;
            real mThickness ;
            real mSurface ;

            // coefficients for shape function
            Matrix< real > mF ;
            Matrix< real > mG ;
            Matrix< real > mH ;

            // pseudo inverse
            Matrix< real > mPseudoInvJ ;

            // gram matrix
            Matrix< real > mGram ;
            Matrix< real > mInvGram ;

            // coefficient matrix
            Matrix< real > mCoeffs ;

            // parameter matrix for curl
            Cell< Matrix< real > > mCurlPars ;

            friend
            const Matrix< real > & maxwell::side_connector_edge_function(
                Calculator * aCalc, const real xi, const real eta, const real zeta );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_PENTA6TS();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_PENTA6TS() override = default ;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement ) override ;

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 ) override;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( point-dependent, computed per integration point;
              *   layout contract in EdgeFunction::mGrad )
              */
            const Matrix< real > &
            G( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_PENTA6TS_HPP