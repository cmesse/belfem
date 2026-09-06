/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_GT_Arguments.hpp"
#include "commtools.hpp"
#include "fn_GT_parse.hpp"
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

                const uint tNumArgs = mArguments.size();

                for ( const string & tArg: mArguments )
                {
                    if ( tArg == "-h" || tArg == "--help" )
                    {
                        mState = State::PrintHelp;
                        break;
                    }
                    else if ( tArg == "-V" || tArg == "--version" )
                    {
                        mState = State::PrintBanner;
                        break;
                    }
                    else if ( tArg == "-g" || tArg == "--gas" )
                    {
                        BELFEM_ERROR( tCount + 1 < tNumArgs,
                                "flag %s needs a value", tArg.c_str() );
                        mGasName = mArguments( tCount + 1 );
                    }
                    else if ( tArg == "-d" || tArg == "--deltaT" )
                    {
                        BELFEM_ERROR( tCount + 1 < tNumArgs,
                                "flag %s needs a value", tArg.c_str() );
                        mDeltaT = parse_real( mArguments( tCount + 1 ), "--deltaT" );
                    }
                    else if ( tArg == "-a" || tArg == "--Tmin" )
                    {
                        BELFEM_ERROR( tCount + 1 < tNumArgs,
                                "flag %s needs a value", tArg.c_str() );
                        mTmin = parse_real( mArguments( tCount + 1 ), "--Tmin" );
                    }
                    else if ( tArg == "-b" || tArg == "--Tmax" )
                    {
                        BELFEM_ERROR( tCount + 1 < tNumArgs,
                                "flag %s needs a value", tArg.c_str() );
                        mTmax = parse_real( mArguments( tCount + 1 ), "--Tmax" );
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