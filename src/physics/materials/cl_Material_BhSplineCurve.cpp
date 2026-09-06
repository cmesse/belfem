/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "constants.hpp"
#include "commtools.hpp"
#include "cl_Material_BhSplineCurve.hpp"
#include "cl_HDF5.hpp"


namespace belfem
{
    namespace material
    {

        BhSplineCurve::BhSplineCurve( const string & aPath, const string & aLabel ) :
            BhCurve( aPath, aLabel )
        {
            this->load_data( aPath, aLabel );
        }

        BhSplineCurve::~BhSplineCurve()
        {
            if( mNuSpline != nullptr ) delete mNuSpline ;
            if( mMuSpline != nullptr ) delete mMuSpline ;
        }

        void
        BhSplineCurve::load_data( const string & aPath, const string & aLabel )
        {
            Vector< real > tRdata( 2 );

            if( mCommRank == 0 )
            {
                HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

                tFile.select_group( aLabel );
                tFile.select_group( "bnur");;
                mNuSpline = new Spline( tFile.active_group(), 0 );
                tFile.close_active_group() ;


                tFile.select_group( "hnur" );
                mMuSpline = new Spline( tFile.active_group(), 0 );
                tFile.close_active_group() ;



                tFile.load_data( "bsat", tRdata( 0 ) );
                tFile.load_data( "hsat", tRdata( 1 ) );

                tFile.close_active_group() ;
                tFile.close() ;
                comm_barrier() ;

                broadcast( tRdata );
            }
            else
            {
                mNuSpline = new Spline( 0 );
                mMuSpline = new Spline( 0 );

                comm_barrier() ;
                broadcast( tRdata );
            }

            mBsat = tRdata( 0 );
            mHsat = tRdata( 1 );
            mMsat = mBsat * constant::nu0 - mHsat ;

        }
    }
}