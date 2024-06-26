//
// Created by Christian Messe on 18.01.22.
//

#ifndef BELFEM_CL_IF_INTEGRATIONDATA_HPP
#define BELFEM_CL_IF_INTEGRATIONDATA_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "Mesh_Enums.hpp"
#include "en_IntegrationScheme.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_Facet.hpp"

namespace belfem
{
    namespace fem
    {
        class IntegrationData
        {
            const ElementType mElementType ;

            bool mOwnShapeFunction = false ;
            InterpolationFunction * mShapeFunction ;

            Vector< real > mWeights ;
            Matrix< real > mPoints ;

            // nodal shape function
            Cell< Vector< real > > mPhi ;
            Cell< Vector< real > > mdPhidxi ; // line elements only
            Cell< Matrix< real > > mN ;
            Cell< Matrix< real > > mNvector ;
            Cell< Matrix< real > > mdNdXi;
            Cell< Matrix< real > > md2NdXi2;

            uint mNumberOfIntegrationPoints = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IntegrationData( const ElementType aElementType,
                             const InterpolationType aType = InterpolationType::LAGRANGE,
                             InterpolationFunction * aShapeFunction = nullptr );

//------------------------------------------------------------------------------

            ~IntegrationData();

//------------------------------------------------------------------------------

            /**
             * default function for popularization
             */
            void
            populate( const uint aIntegrationOrder,
                      const IntegrationScheme aScheme = IntegrationScheme::GAUSS  );

//------------------------------------------------------------------------------

            /**
             * special function for popularization if this is a sideset
             */
            void
            populate_for_master( const uint aMasterIndex,
                                 const uint aIntegrationOrder=0,
                                 const IntegrationScheme aScheme = IntegrationScheme::GAUSS );

//------------------------------------------------------------------------------

            /**
             * special function for popularization if this is a sideset
             */
            void
            populate_for_slave( const uint aSlaveIndex,
                                const uint aOrientation=0,
                                const uint aIntegrationOrder=0,
                                const IntegrationScheme aScheme = IntegrationScheme::GAUSS  );


//------------------------------------------------------------------------------

            /**
             * access the shape function
             */
             const InterpolationFunction *
             function() const ;

//------------------------------------------------------------------------------

            /**
             * tells how many integration points are used
             */
            uint
            number_of_integration_points() const ;

//------------------------------------------------------------------------------

             /**
              * return the integration weights
              */
             const Vector< real > &
             weights() const ;

//------------------------------------------------------------------------------

             /**
              * return the integration points
              */
             const Matrix< real > &
             points() const ;

//------------------------------------------------------------------------------

            /**
             * return the node shape function as vector
             */
             const Vector< real > &
             phi( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * return the node shape function deriivative as vector (line element only)
             */
            const Vector< real > &
            dphidxi( const uint aIndex ) const ;

//------------------------------------------------------------------------------

             /**
              * return the node shape function as Matrix (scalar field)
              */
             const Matrix< real > &
             N( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * return the node shape function as Matrix, for vector field
             */
            const Matrix< real > &
            Nvector( const uint aIndex ) const ;

//------------------------------------------------------------------------------

             /**
              * return the first derivative of N
              */
             const Matrix< real > &
             dNdXi( const uint aIndex ) const ;

//------------------------------------------------------------------------------

             /**
              * return the second derivative of N
              */
             const Matrix< real > &
             d2NdXi2( const uint aIndex ) const ;

//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

            /**
             * special function for popularization if this is a sideset
             */
            void
            populate_for_slave_tri( const uint aFaceIndex,
                                    const uint aIntegrationOrder,
                                    const IntegrationScheme aScheme  );

//------------------------------------------------------------------------------

            /**
             * special function for popularization if this is a sideset
             */
            void
            populate_for_slave_tet(
                    const uint aSlaveIndex,
                    const uint aOrientation,
                    const uint aIntegrationOrder,
                    const IntegrationScheme aScheme );

//------------------------------------------------------------------------------

            /**
             * special function for popularization if this is a sideset
             */
            void
            populate_for_slave_hex(
                    const uint aSlaveIndex,
                    const uint aOrientation,
                    const uint aIntegrationOrder,
                    const IntegrationScheme aScheme );

//------------------------------------------------------------------------------

            void
            evaluate_function();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

            inline const InterpolationFunction *
            IntegrationData::function() const
            {
                return mShapeFunction ;
            }

//------------------------------------------------------------------------------

             inline uint
             IntegrationData::number_of_integration_points() const
             {
                return mNumberOfIntegrationPoints ;
             }

//------------------------------------------------------------------------------

            inline const Vector< real > &
            IntegrationData::weights() const
            {
                return mWeights ;
            }

//------------------------------------------------------------------------------

            inline const Matrix< real > &
            IntegrationData::points() const
            {
                return mPoints ;
            }

//------------------------------------------------------------------------------

            inline const Vector< real > &
            IntegrationData::phi( const uint aIndex ) const
            {
                return mPhi( aIndex ) ;
            }
//------------------------------------------------------------------------------

            inline const Vector< real > &
            IntegrationData::dphidxi( const uint aIndex ) const
            {
                BELFEM_ASSERT( mdPhidxi.size() > 0,
                "dphidxi can only be called for line elements" );

                return mdPhidxi( aIndex ) ;
            }

//------------------------------------------------------------------------------

            inline const Matrix< real > &
            IntegrationData::N( const uint aIndex ) const
            {
                return mN( aIndex ) ;
            }

//------------------------------------------------------------------------------

            inline const Matrix< real > &
            IntegrationData::Nvector( const uint aIndex ) const
            {
                return mNvector( aIndex ) ;
            }

//------------------------------------------------------------------------------

            inline const Matrix< real > &
            IntegrationData::dNdXi( const uint aIndex ) const
            {
                return mdNdXi( aIndex );
            }

//------------------------------------------------------------------------------

            inline const Matrix< real > &
            IntegrationData::d2NdXi2( const uint aIndex ) const
            {
                return md2NdXi2( aIndex );
            }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_INTEGRATIONDATA_HPP
