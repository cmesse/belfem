//
// Created by gregorygiard on 1/12/26.
//

#include "cl_TimestepMatrices.hpp"

namespace belfem
{
    namespace fem
    {

//----------------------------------------------------------------------------------------

        void
        TimestepMatrices::initialize( const index_t aNumDofs )
        {
            if (mFlags.test( static_cast< index_t >( MatrixFlag::M )))
            {
                mM.set_size( aNumDofs, aNumDofs, 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::D )))
            {
                mD.set_size( aNumDofs, aNumDofs, 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::K )))
            {
                mK.set_size( aNumDofs, aNumDofs, 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::F )))
            {
                mF.set_size( aNumDofs, 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )))
            {
                mdMdX_times_x.set_size( aNumDofs, aNumDofs, 0.0  ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_h )))
            {
                mdMdX_times_h.set_size( aNumDofs, aNumDofs, 0.0  ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dKdX_times_x )))
            {
                mdKdX_times_x.set_size( aNumDofs, aNumDofs, 0.0  ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dFdX )))
            {
                mdFdX.set_size( aNumDofs, aNumDofs, 0.0  ) ;
            }

            mJ.set_size( aNumDofs, aNumDofs, 0.0  ) ;
            mdJdx.set_size( aNumDofs, aNumDofs, 0.0  ) ;
        }

//----------------------------------------------------------------------------------------

        void TimestepMatrices::reset()
        {
            if (mFlags.test( static_cast< index_t >( MatrixFlag::M )))
            {
                mM.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::D )))
            {
                mD.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::K )))
            {
                mK.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::F )))
            {
                mF.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )))
            {
                mdMdX_times_x.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_h )))
            {
                mdMdX_times_h.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dKdX_times_x )))
            {
                mdKdX_times_x.fill( 0.0 ) ;
            }

            if (mFlags.test( static_cast< index_t >( MatrixFlag::dFdX )))
            {
                mdFdX.fill( 0.0 ) ;
            }

            mJ.fill( 0.0 ) ;
            mdJdx.fill( 0.0 ) ;
        }

//----------------------------------------------------------------------------------------

        void
        TimestepMatrices::assemble_J( const real adt )
        {
            mJ.fill( 0.0 ) ;

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::M )) )
            {
                mJ += mM ;
            }

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::K )) )
            {
                mJ += mK*adt ;
            }

        }

//----------------------------------------------------------------------------------------

        void
        TimestepMatrices::assemble_dJdx( const real adt, const real aAlpha )
        {
            mdJdx.fill( 0.0 ) ;

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_x )) )
            {
                mdJdx += mdMdX_times_x*aAlpha ;
            }

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::dMdX_times_h )) )
            {
                mdJdx -= mdMdX_times_h ;
            }

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::dKdX_times_x )) )
            {
                mdJdx += mdKdX_times_x*adt ;
            }

            if ( mFlags.test( static_cast< index_t >( MatrixFlag::dFdX )) )
            {
                mdJdx -= mdFdX*adt ;
            }

        }

//----------------------------------------------------------------------------------------


    }
}