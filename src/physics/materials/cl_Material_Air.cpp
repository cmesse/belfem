//
// Created by Christian Messe on 03.08.20.
//

#include "commtools.hpp"

#include "cl_Gas.hpp"
#include "cl_Material_Air.hpp"
namespace belfem
{
    namespace material
    {
        Air::Air() :
                mHeatSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
                mConductivitySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax )
        {

            // set maximum temperature
            mTmax = 2200.0;

            this->create_splines() ;

            mHasThermal = true;
        }

//----------------------------------------------------------------------------

        void
        Air::create_splines()
        {
            if( comm_rank() == 0 )
            {
                // if not defined otherwise, gas is automatically air
                Gas tAir;

                // copy splines
                mHeatSpline.matrix_data()
                    = tAir.heat_spline().matrix_data();

                mConductivitySpline.matrix_data()
                    = tAir.conductivity_spline().matrix_data();

                mR = tAir.R( BELFEM_TREF, BELFEM_PREF );
            }

            // synchronize data with other procs
            mHeatSpline.broadcast( 0 );
            mConductivitySpline.broadcast( 0 );
            broadcast( 0, mR );

        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */