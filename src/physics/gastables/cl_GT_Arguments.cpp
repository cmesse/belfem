//
// Created by Christian Messe on 01.09.19.
//

#include "cl_GT_Arguments.hpp"
#include "commtools.hpp"
namespace belfem
{
    namespace gastables
    {

//------------------------------------------------------------------------------

        Arguments::Arguments( int & argc, char * argv[] ) :
                belfem::Arguments( argc, argv )
        {
            if ( comm_rank() == 0 )
            {
                this->check_arguments();
            }
        }

//------------------------------------------------------------------------------

        void
        Arguments::check_arguments()
        {
            if ( mArguments.size() <= 1 )
            {
                mState = State::PrintUsage;
            }
            else
            {
                mState = State::PrintTable ;

                uint tCount = 0;

                for ( string tArg: mArguments )
                {
                    if ( tArg == "-h" || tArg == "--help" )
                    {
                        mState = State::PrintHelp;
                        break;
                    }
                    else if ( tArg == "-g" || tArg == "--gas" )
                    {
                        mGasName = mArguments( tCount + 1 );
                    }
                    else if ( tArg == "-d" || tArg == "--deltaT" )
                    {
                        mDeltaT = std::stod( mArguments( tCount + 1 ));
                    }
                    else if ( tArg == "-a" || tArg == "--Tmin" )
                    {
                        mTmin = std::stod( mArguments( tCount + 1 ));
                    }
                    else if ( tArg == "-b" || tArg == "--Tmax" )
                    {
                        mTmax = std::stod( mArguments( tCount + 1 ));
                    }
                    else if ( tArg == "-m" || tArg == "--molar" )
                    {
                        mMolarFlag = true;
                    }
                    tCount++;
                }
            }
        }

//------------------------------------------------------------------------------

        const State &
        Arguments::state() const
        {
            return mState;
        }

//------------------------------------------------------------------------------

        const real &
        Arguments::T_min() const
        {
            return mTmin;
        }

//------------------------------------------------------------------------------

        const real &
        Arguments::T_max() const
        {
            return mTmax;
        }

//------------------------------------------------------------------------------

        const real &
        Arguments::delta_T() const
        {
            return mDeltaT;
        }

//------------------------------------------------------------------------------

        const string &
        Arguments::gasname() const
        {
            return mGasName;
        }

//------------------------------------------------------------------------------

        bool
        Arguments::molar() const
        {
            return mMolarFlag;
        }

//------------------------------------------------------------------------------
    }
}