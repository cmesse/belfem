//
// Created by christian on 12/3/21.
//

#ifndef BELFEM_CL_EF_TET4_HPP
#define BELFEM_CL_EF_TET4_HPP


#include "cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_TET4 : public EdgeFunction
        {
            real mNablaXi[ 3 ];
            real mNablaEta[ 3 ];
            real mNablaZeta[ 3 ];
            real mNablaTau[ 3 ];
            real mS[ 6 ];

            Matrix< real > mG ;
            Matrix< real > mH ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_TET4();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_TET4() = default;

//------------------------------------------------------------------------------

            /**
             * links the shape function with the element and precomputes data
             * @param aElement
             * @param aLinkB    if the B-Operator is used (phi and az)
             * @param aLinkBA   if the B-A-Operator is used (ax, ay, az)
             * @param aLinkC   if the curl operator is used
             */
            /**
             * links the shape function with the element and precomputes data
             * @param aElement
             * @param aOperatorCurlH    if the B-Operator is used (phi and az)
             * @param aOperatorGradPhi  if the B-A-Operator is used (ax, ay, az)
             * @param aOperatorCurlA    if the curl operator is used
             */
            virtual void
            link( Element * aElement,
                  const bool aOperatorCurlH=true,
                  const bool aOperatorGradPhi=false,
                  const bool aOperatorCurlA=false );

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi );

//------------------------------------------------------------------------------

            // compute the edge function
            const Matrix< real > &
            E( const uint aIndex );

//------------------------------------------------------------------------------

            // compute the curl function
            const Matrix< real > &
            C( const uint aIndex = 0 );

//------------------------------------------------------------------------------

            // compute the B-Operator
            const Matrix< real > &
            B( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------

            // compute the Curl operator for the a-field
            const Matrix< real > &
            CA( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------
        } ;

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET4::B( const uint aIndex )
        {
            return mB ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET4::CA( const uint aIndex )
        {
            return mCurlA ;
        }
//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET4::C( const uint aIndex )
        {
            return mC ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_EF_TET4_HPP
