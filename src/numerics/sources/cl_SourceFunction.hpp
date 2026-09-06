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

#ifndef BELFEM_CL_SOURCEFUNCTION_HPP
#define BELFEM_CL_SOURCEFUNCTION_HPP

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#pragma clang diagnostic ignored "-Wformat"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"

#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#elif defined(BELFEM_INTEL)
    // Intel compiler diagnostics
#pragma warning(push)
    // disable format warnings
#pragma warning(disable: 1011)
    // disable unused variable warnings
#pragma warning(disable: 177)
    // keep your original one
#pragma warning(disable: 1595)

#endif

#include "typedefs.hpp"
#include "cl_Cell.hpp"

// note: NEVER include cl_Vector.hpp or cl_Matrix.hpp here, or any class that uses it.
//       Doing so would break the API for the user defined materials.

#define BELFEM_BCVAL_AMPLITUDE  0
#define BELFEM_BCVAL_PERIOD     1
#define BELFEM_BCVAL_FREQUENCY  2
#define BELFEM_BCVAL_OMEGA      3
#define BELFEM_BCVAL_PHASE      4
#define BELFEM_BCVAL_TIMEOFFSET 5
#define BELFEM_BCVAL_FUZZYNESS  6
namespace belfem
{
//-----------------------------------------------------------------------------

    enum class SourceFunctionType
    {
        Constant    = 0,
        Ramp        = 1,
        Sigmoid     = 2,
        Sine        = 3,
        Square      = 4,
        Triangle    = 5,
        Sawtooth    = 6,
        UserDefined = 7,
        UNDEFINED   = 8
    };

//-----------------------------------------------------------------------------

    SourceFunctionType
    boundary_condition_function_type( const string & aString );

//------------------------------------------------------------------------------

    typedef real ( UserFunc )( const real ) ;

    class SourceFunction
    {
        const proc_t mCommRank ;
        const proc_t mCommSize ;

        SourceFunctionType mType = SourceFunctionType::UNDEFINED ;

        Cell< real > mValues ;

        real
        ( SourceFunction::*mFun )( const real aTime ) const ;

        UserFunc *mUserFunction = nullptr ;

        //! dlopen handle of the user plugin, null for the builtin waveforms.
        //! Owned: closed by the destructor, which is why the copy and move
        //! operations below are deleted
        void *mHandle = nullptr ;

//-----------------------------------------------------------------------------
    public:
//-----------------------------------------------------------------------------

        SourceFunction();

//-----------------------------------------------------------------------------

        ~SourceFunction() ;

//-----------------------------------------------------------------------------

        // This object can own a dlopen handle ( mHandle ) and a function
        // pointer into that library ( mUserFunction ). Copying either would
        // give two objects one mapping and close it twice. Nothing in the tree
        // copies or moves one -- production always allocates and passes a
        // pointer, and the tests default-construct on the stack, which needs
        // neither -- so deleting these costs nothing today and keeps it that way
        SourceFunction( const SourceFunction & ) = delete ;
        SourceFunction & operator=( const SourceFunction & ) = delete ;
        SourceFunction( SourceFunction && ) = delete ;
        SourceFunction & operator=( SourceFunction && ) = delete ;

//-----------------------------------------------------------------------------

        SourceFunctionType
        type() const ;

//-----------------------------------------------------------------------------

        void
        set_constant( const real aAmplitude );

//-----------------------------------------------------------------------------

        void
        set_ramp(
                const real aAmplitude,
                const real aPeriod,
                const real aTimeOffset=0.0 );

//-----------------------------------------------------------------------------

        void
        set_sigmoid(
                const real aAmplitude,
                const real aPeriod,
                const real aTimeOffset,
                const real aFuzzyness=0.01 );

//-----------------------------------------------------------------------------

        void
        set_periodic(
                const SourceFunctionType aType,
                const real aAmplitude,
                const real aPeriod,
                const real aPhase );

//-----------------------------------------------------------------------------

        void
        set_user_defined( UserFunc * Function ) ;

//-----------------------------------------------------------------------------

        void
        read_user_defined( const string & aLibraryPath,
                            const string & aLabel ) ;

//-----------------------------------------------------------------------------

        real
        compute( const real aTime ) const ;

//-----------------------------------------------------------------------------

        real
        amplitude() const ;

        real
        period() const ;

        real
        frequency() const ;

        real
        omega() const ;

        real
        phase() const ;

        real
        time_offset() const ;

        real
        fuzziness() const ;

//-----------------------------------------------------------------------------
    private:
//-----------------------------------------------------------------------------

        void
        set_defaults();

//-----------------------------------------------------------------------------

        void
        synch();

//-----------------------------------------------------------------------------

        void
        set_amplitude( const real aValue ) ;

        void
        set_period(  const real aValue ) ;

        void
        set_frequency(  const real aValue ) ;

        void
        set_omega(  const real aValue ) ;

        void
        set_phase( const real aValue ) ;

        void
        set_time_offset(  const real aValue ) ;

        void
        set_fuzziness(  const real aValue ) ;

//-----------------------------------------------------------------------------

        real
        function_constant( const real aTime ) const ;

        real
        function_ramp( const real aTime ) const ;

        real
        function_sigmoid( const real aTime ) const ;

        real
        function_sine( const real aTime ) const ;

        /**
         * normalized time within one period, on [0,1)
         *
         * std::fmod keeps the sign of its dividend, so a shifted time that
         * lands before the start of the period comes back negative and the
         * waveform below it reads off the wrong branch. Shifting a negative remainder up by one
         * period gives the periodic extension f(t) == f(t+T).
         * Precondition: a positive, finite period -- a zero, infinite or
         * negative period keeps whatever fmod returns, exactly as before
         */
        real
        wrap_tau( const real aTime ) const ;

        real
        function_square( const real aTime ) const ;

        real
        function_triangle( const real aTime ) const ;

        real
        function_sawtooth( const real aTime ) const ;

        real
        function_userdefined( const real aTime ) const ;

    };

//-----------------------------------------------------------------------------

    inline SourceFunctionType
    SourceFunction::type() const
    {
        return mType ;
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::compute( const real aTime ) const
    {
        return ( this->*mFun )( aTime ) ;
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::amplitude() const
    {
        return mValues( BELFEM_BCVAL_AMPLITUDE ) ;
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::period() const
    {
        return mValues( BELFEM_BCVAL_PERIOD );
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::frequency() const
    {
        return mValues( BELFEM_BCVAL_FREQUENCY );
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::omega() const
    {
        return mValues( BELFEM_BCVAL_OMEGA );
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::phase() const
    {
        return mValues( BELFEM_BCVAL_PHASE );
    }


//-----------------------------------------------------------------------------

    inline real
    SourceFunction::time_offset() const
    {
        return mValues( BELFEM_BCVAL_TIMEOFFSET );
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::fuzziness() const
    {
        return mValues( BELFEM_BCVAL_FUZZYNESS );
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::function_constant( const real aTime ) const
    {
        return this->amplitude() ;
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::function_ramp( const real aTime ) const
    {
        if( aTime < this->time_offset() )
        {
            return 0.0 ;
        }
        else if ( aTime > this->time_offset() + this->period() )
        {
            return this->amplitude() ;
        }
        else
        {
            return this->amplitude() * ( aTime - this->time_offset() ) * this->frequency() ;
        }
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::function_sine( const real aTime ) const
    {
        return this->amplitude() * std::sin( this->omega() * aTime  + this->phase() );
    }

//------------------------------------------------------------------------------

    inline real
    SourceFunction::function_sigmoid( const real aTime ) const
    {
        // see doi:10.1016/j.fss.2005.02.016

        real d = this->time_offset() + this->period()/2.0;
        real beta = -2.0*log(this->fuzziness()/(1.0-this->fuzziness()))/this->period() ;
        real f = 1./(1+std::exp(beta * (d - aTime)));

        return this->amplitude() * f ;
    }

//------------------------------------------------------------------------------

    inline real
    SourceFunction::wrap_tau( const real aTime ) const
    {
        real tTau = std::fmod( aTime, this->period() ) / this->period() ;

        // strictly less: an exact zero must stay zero, or the square would
        // flip at its own rising edge
        return tTau < 0.0 ? tTau + 1.0 : tTau ;
    }

//------------------------------------------------------------------------------

    inline real
    SourceFunction::function_square( const real aTime ) const
    {
        real tTau = this->wrap_tau( aTime + this->time_offset() ) ;


        if( tTau < 0.5 )
        {
            return this->amplitude() ;
        }
        else
        {
            return -this->amplitude() ;
        }
    }

//------------------------------------------------------------------------------

    inline real
    SourceFunction::function_triangle( const real aTime ) const
    {
        real tTau = this->wrap_tau( aTime + this->time_offset() + 3.0*this->period()/4.0 );

        // Shift tTau to range from -0.5 to 0.5
        tTau -= 0.5;

        // Calculate triangle wave value
        return 4.0 * this->amplitude() * std::fabs(tTau)  - this->amplitude();
    }

//------------------------------------------------------------------------------

    inline real
    SourceFunction::function_sawtooth( const real aTime ) const
    {
        return this->amplitude()
             * this->wrap_tau( aTime + this->time_offset() ) ;
    }

//-----------------------------------------------------------------------------

    inline real
    SourceFunction::function_userdefined( const real aTime ) const
    {
        return mUserFunction(aTime) ;
    }

    //-----------------------------------------------------------------------------
}

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_INTEL
#pragma warning pop
#endif

#endif //BELFEM_CL_SOURCEFUNCTION_HPP
