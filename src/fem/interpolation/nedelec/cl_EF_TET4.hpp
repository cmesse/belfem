//
// Created by christian on 12/3/21.
//

#ifndef BELFEM_CL_EF_TET4_HPP
#define BELFEM_CL_EF_TET4_HPP


#include "nedelec/cl_EF_EdgeFunction.hpp"
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
            ~EF_TET4() override = default;

//------------------------------------------------------------------------------

            /**
             * links the shape function with the element and precomputes data
             * @param aElement
             */
            void
            link( Element * aElement ) override ;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            // compute the edge function
            const Matrix< real > &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            // compute the curl function
            const Matrix< real > &
            C( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( constant for this element, filled in link();
              *   layout contract in EdgeFunction::mGrad )
              */
            const Matrix< real > &
            G( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------
        } ;

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET4::C( const uint aIndex )
        {
            return mC ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_TET4::G( const uint aIndex )
        {
            return mGrad ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_EF_TET4_HPP
