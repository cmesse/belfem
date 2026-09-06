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

#include <dlfcn.h>

#include "constants.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "stringtools.hpp"
#include "filetools.hpp"
#include "fn_sprint.hpp"

#include "cl_SourceFunction.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    SourceFunctionType
    boundary_condition_function_type( const string & aString )
    {
        string tString = string_to_lower( aString );

        if ( tString == "constant" )
        {
            return SourceFunctionType::Constant;
        }
        else if ( tString == "ramp" )
        {
            return SourceFunctionType::Ramp;
        }
        else if ( tString == "sigmoid" )
        {
            return SourceFunctionType::Sigmoid;
        }
        else if ( tString == "sine" )
        {
            return SourceFunctionType::Sine;
        }
        else if ( tString == "square" )
        {
            return SourceFunctionType::Square;
        }
        else if ( tString == "triangle" )
        {
            return SourceFunctionType::Triangle;
        }
        else if ( tString == "sawtooth" )
        {
            return SourceFunctionType::Sawtooth;
        }
        else if ( tString == "userdefined" )
        {
            return SourceFunctionType::UserDefined;
        }
        else
        {
            BELFEM_ERROR( false, "Unknown Boundary Condition Function Type: %s", aString.c_str());
            return SourceFunctionType::UNDEFINED;
        }
    }

//-----------------------------------------------------------------------------

    SourceFunction::SourceFunction() :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() )
    {
        this->set_defaults() ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_defaults()
    {
        mValues.set_size( 7, 0.0 );

        // amplitude
        mValues( BELFEM_BCVAL_AMPLITUDE ) = BELFEM_QUIET_NAN ;

        // period
        mValues( BELFEM_BCVAL_PERIOD ) = BELFEM_INFINITY ;

        // frequency
        mValues( BELFEM_BCVAL_FREQUENCY ) = 0.0 ;

        // omega
        mValues( BELFEM_BCVAL_OMEGA ) = 0.0 ;

        // phase
        mValues( BELFEM_BCVAL_PHASE ) = 0.0 ;

        // time offset
        mValues( BELFEM_BCVAL_TIMEOFFSET  ) = 0.0 ;

        // fuzzyness
        mValues( BELFEM_BCVAL_FUZZYNESS ) = BELFEM_QUIET_NAN ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_constant( const real aAmplitude )
    {
        this->set_defaults() ;
        this->set_amplitude( aAmplitude );
        this->synch() ;

        mType = SourceFunctionType::Constant ;
        mFun = & SourceFunction::function_constant ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_ramp(
            const real aAmplitude,
            const real aPeriod,
            const real aTimeOffset )
    {
        this->set_defaults() ;
        this->set_amplitude( aAmplitude );
        this->set_period( aPeriod );
        this->set_time_offset( aTimeOffset );
        this->synch() ;

        mType = SourceFunctionType::Ramp ;
        mFun = & SourceFunction::function_ramp ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_sigmoid(
            const real aAmplitude,
            const real aPeriod,
            const real aTimeOffset,
            const real aFuzzyness )
    {
        this->set_defaults() ;
        this->set_amplitude( aAmplitude );
        this->set_period( aPeriod );
        this->set_time_offset( aTimeOffset );
        this->set_fuzziness( aFuzzyness );
        this->synch() ;
        mType = SourceFunctionType::Sigmoid ;
        mFun = & SourceFunction::function_sigmoid ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_periodic(
            const SourceFunctionType aType,
            const real aAmplitude,
            const real aPeriod,
            const real aPhase )
    {
        this->set_defaults() ;
        this->set_amplitude( aAmplitude );
        this->set_period( aPeriod );
        // set_phase writes BOTH slots: the phase itself, which the sine
        // reads, and the time offset phase/omega, which triangle, square
        // and sawtooth read. No second conversion belongs here
        this->set_phase( aPhase );

        switch ( aType )
        {
            case( SourceFunctionType::Sine ) :
            {
                mFun = & SourceFunction::function_sine ;
                break ;
            }
            case( SourceFunctionType::Square ) :
            {
                mFun = & SourceFunction::function_square ;
                break ;
            }
            case( SourceFunctionType::Triangle ) :
            {
                mFun = & SourceFunction::function_triangle ;
                break ;
            }
            case( SourceFunctionType::Sawtooth ) :
            {
                mFun = & SourceFunction::function_sawtooth ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid BC type");
            }
        }

        mType = aType ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_user_defined(
              UserFunc * Function )
    {
        mUserFunction = Function ;
        mFun = & SourceFunction::function_userdefined ;
        mType = SourceFunctionType::UserDefined;
    }
//-----------------------------------------------------------------------------

    void
    SourceFunction::read_user_defined( const string & aLibraryPath,
                                        const string & aLabel )
    {
        // one plugin per object. A second call would strand the first mapping,
        // and mUserFunction would silently start pointing into a library this
        // object no longer knows it holds
        BELFEM_ERROR( mHandle == nullptr,
            "SourceFunction already loaded a plugin; cannot load %s on top of it",
            aLibraryPath.c_str() );

        // the same search the material and defect plugins get: run directory
        // first, then the material subdirectory of the shared data path. An
        // unresolved name is handed on unchanged, so dlopen still gets its
        // $LD_LIBRARY_PATH search for a name carrying no slash
        const string tPath = search_data_file( aLibraryPath, "material" );

        mHandle = dlopen( tPath.c_str(), RTLD_LAZY );

        BELFEM_ERROR( mHandle != nullptr, "Could not load library %s ( error: %s )",
            tPath.c_str(), dlerror() );

        string tFunctionName = sprint( "%s_init", aLabel.c_str() );

        void (* tInitFunc )( SourceFunction * ) =
            reinterpret_cast< void(*)( SourceFunction * ) >( dlsym( mHandle, tFunctionName.c_str() ) );

        BELFEM_ERROR( tInitFunc != nullptr, "Could not find function %s in library %s ( error: %s )",
             tFunctionName.c_str(), tPath.c_str(), dlerror() );

        tInitFunc( this );
    }

//-----------------------------------------------------------------------------

    SourceFunction::~SourceFunction()
    {
        // mUserFunction points into the plugin's text, so the mapping must
        // outlive every compute() call. It does: the owners delete this object
        // at the end of the boundary condition's or component's life, after
        // which nothing can reach mUserFunction any more.
        //
        // dlclose( nullptr ) is not a POSIX no-op, and mHandle is null for the
        // seven builtin waveforms, so the guard is required rather than tidy.
        // The Maxwell factory builds one SourceFunction per group member from
        // the same ( file, label ), which is N dlopen calls on one library;
        // dlopen refcounts, so N guarded closes are exactly right
        if( mHandle != nullptr )
        {
            dlclose( mHandle );
            mHandle = nullptr ;
        }
    }

    void
    SourceFunction::synch()
    {
        broadcast( mValues );
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_amplitude( const real aValue )
    {
        mValues( BELFEM_BCVAL_AMPLITUDE ) = aValue ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_period(  const real aValue )
    {
        if( std::abs( aValue ) < BELFEM_EPSILON )
        {
            mValues( BELFEM_BCVAL_PERIOD ) = 0.0 ;

            // frequency
            mValues( BELFEM_BCVAL_FREQUENCY ) = BELFEM_INFINITY ;

            // omega
            mValues( BELFEM_BCVAL_OMEGA ) = BELFEM_INFINITY ;
        }
        else
        {
            // period
            mValues( BELFEM_BCVAL_PERIOD ) = aValue ;

            // frequency
            mValues( BELFEM_BCVAL_FREQUENCY ) = 1.0 / aValue ;

            // omega
            mValues( BELFEM_BCVAL_OMEGA ) = 2.0 * constant::pi / aValue ;
        }
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_frequency(  const real aValue )
    {
        if( std::abs( aValue ) < BELFEM_EPSILON )
        {
            // period
            mValues( BELFEM_BCVAL_PERIOD ) = BELFEM_INFINITY ;

            // frequency
            mValues( BELFEM_BCVAL_FREQUENCY ) = 0.0 ;

            // omega
            mValues( BELFEM_BCVAL_OMEGA ) = 0.0 ;
        }
        else
        {
            mValues( BELFEM_BCVAL_PERIOD ) = 1.0 / aValue ;
            mValues( BELFEM_BCVAL_FREQUENCY ) = aValue ;
            mValues( BELFEM_BCVAL_OMEGA ) = 2.0 * constant::pi * aValue ;
        }
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_omega(  const real aValue )
    {
        if( std::abs( aValue ) < BELFEM_EPSILON )
        {
            // period
            mValues( BELFEM_BCVAL_PERIOD ) = BELFEM_INFINITY ;

            // frequency
            mValues( BELFEM_BCVAL_FREQUENCY ) = 0.0 ;

            // omega
            mValues( BELFEM_BCVAL_OMEGA ) = 0.0 ;
        }
        else
        {
            // period
            mValues( BELFEM_BCVAL_PERIOD ) = 2.0 * constant::pi / aValue ;

            // frequency
            mValues( BELFEM_BCVAL_FREQUENCY ) =  0.5 * aValue / constant::pi ;

            // omega
            mValues( BELFEM_BCVAL_OMEGA ) = aValue ;
        }
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_phase( const real aValue )
    {
        // phase
        mValues( BELFEM_BCVAL_PHASE ) = aValue ;

        // time offset
        if( this->omega() == BELFEM_INFINITY )
        {
            mValues( BELFEM_BCVAL_TIMEOFFSET  ) = 0.0 ;
        }
        else if( std::abs( this->omega() ) < BELFEM_EPSILON )
        {
            mValues( BELFEM_BCVAL_TIMEOFFSET  ) = BELFEM_INFINITY ;
        }
        else
        {
            mValues( BELFEM_BCVAL_TIMEOFFSET  ) = aValue / this->omega() ;
        }
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_time_offset( const real aValue )
    {
        // phase
        mValues( BELFEM_BCVAL_PHASE ) = aValue * this->omega() ;

        // time offset
        mValues( BELFEM_BCVAL_TIMEOFFSET  ) = aValue ;
    }

//-----------------------------------------------------------------------------

    void
    SourceFunction::set_fuzziness( const real aValue )
    {
        mValues( BELFEM_BCVAL_FUZZYNESS ) = aValue ;
    }

//-----------------------------------------------------------------------------


}
