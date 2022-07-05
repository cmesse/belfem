//
// Created by Christian Messe on 27.12.19.
//

#include "fn_ODE_RK45.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace ode
    {
//------------------------------------------------------------------------------

        void
        RK45_init( ODE & aODE, Cell< Vector< real > > & aWork )
        {
            // number of entries in Y-vector
            index_t tN = aODE.dimension();

            aWork.set_size( 8, {} );

            // k0 - k5, y=ya, yb
            for( index_t k=0; k<8; ++k )
            {
                aWork( k ).set_size( tN, 0.0 );
            }
        }

//------------------------------------------------------------------------------

        Status
        RK45(    ODE                    & aODE,
                 real                   & aT,
                 Vector< real >         & aY,
                 real                   & aStep,
                 Cell< Vector< real > > & aWork,
                 const real              aEpsilon,
                 const uint              aMaxIterations,
                 const real              aTmax,
                 const bool              aAutoTimestep )
        {
            // test length of work array
            BELFEM_ASSERT( aWork.size() >= 8, "Size of working cell does not match" );

            Vector< real > & k0 = aWork(  0 );
            Vector< real > & k1 = aWork(  1 );
            Vector< real > & k2 = aWork(  2 );
            Vector< real > & k3 = aWork(  3 );
            Vector< real > & k4 = aWork(  4 );
            Vector< real > & k5 = aWork(  5 );

            Vector< real > & y0 = aY;
            Vector< real > & y1 = aWork(  6 );
            Vector< real > & y2 = aWork(  6 );
            Vector< real > & y3 = aWork(  6 );
            Vector< real > & y4 = aWork(  6 );
            Vector< real > & y5 = aWork(  6 );

            Vector< real > & ya = aWork( 6 );
            Vector< real > & yb = aWork( 7 );

            real h = aStep;

            real t0 = aT;
            real t1 ;
            real t2 ;
            real t3 ;
            real t4 ;
            real t5 ;

            // initial timestep
            aODE.compute( t0, y0, k0 );

            real tolb;
            Status aStatus = Status::OK ;

            //real tEvent0 = aODE.check_events( t0, y0 );

            for( uint i=0; i<aMaxIterations; ++i )
            {
                // trap h
                if( t0 + h >= aTmax )
                {
                    h = aTmax - t0;
                    aStatus = Status::TRAPPED ;
                }

                t1 = t0 + 0.25*h;
                y1 = y0 + 0.25*h*k0;

                aODE.compute( t1, y1, k1 );

                t2 = t0 + 0.375*h;
                y2 = y0 + (h/32)*(3*k0 + 9*k1);

                aODE.compute( t2, y2, k2 );

                t3 = t0 + (12./13.)*h;
                y3 = y0 + (h/2197.)*(1932.*k0 - 7200.*k1 + 7296.*k2);

                aODE.compute( t3, y3, k3 );

                t4 = t0 + h;
                y4 = y0 + (h/4104.)*(8341.*k0 - 32832.*k1 + 29440.*k2 - 845.*k3);

                aODE.compute( t4, y4, k4 );

                t5 = t0 + 0.5*h;

                y5 = y0 + (h/20520.)*(-6080.*k0 + 41040.*k1 - 28352.*k2 + 9295.*k3 - 5643.*k4);

                aODE.compute( t5, y5, k5 );

                // fourth order solution
                ya = y0 + (h/20520.)*(2375.*k0 + 11264.*k2 + 10985.*k3 - 4104.*k4);

                // fifth order solution
                yb = y0 + (h/282150.)*(33440.*k0 + 146432.*k2 + 142805.*k3 - 50787.*k4 + 10260.*k5);

                tolb = 0.0;

                // compute the norm
                for( uint k=0; k<ya.length(); ++k )
                {
                    tolb += std::pow( ya( k ) - yb( k ), 2 );
                }

                tolb = std::sqrt( tolb ) / h;

                /*
                 *  real tEvent1 = aODE.check_events( t1, yb );
                 *
                 *  if ( ( tEvent0 * tEvent1 < 0 ) && ( std::abs( tEvent0 - tEvent1 ) > tEpsilon )
                 *  {
                 *      --> adapt timestep using Regula Falsi
                 *      aStatus = Status::TRAPPED_EVENT;
                 *  }
                 */

                if( aAutoTimestep )
                {
                    // adapt timestep
                    if ( tolb > 0 ) // todo: why is this zero ?
                    {
                        h *= 0.9 * std::pow( aEpsilon / tolb, 0.25 );
                    }
                    else
                    {
                        h *= 1.05;
                    }
                }

                if( tolb < aEpsilon )
                {
                    aT = t4;
                    aY.vector_data() = yb.vector_data();

                    if( aStatus == Status::OK && aAutoTimestep )
                    {
                        aStep = h;
                    }
                    return aStatus;
                }

            }

            return Status::MAXIT ;
        }

//------------------------------------------------------------------------------
    }
}