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

#include <cstdio>
#include <iostream>

#include "cl_Progressbar.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    Progressbar::Progressbar( const uint       aNumSteps, FILE * aFile ) :
        mNumSteps( aNumSteps ),
        mFile( aFile )
    {
        this->reset();
    }

//------------------------------------------------------------------------------

    Progressbar::~Progressbar()
    {
        // reset() hid the cursor. Give it back even if the caller never
        // reached finish(), or the terminal stays blind after the run.
        this->show_cursor();
    }

//------------------------------------------------------------------------------
    void
    Progressbar::reset()
    {
        mStep     = 0;
        mProgress = 0;
        mCursorRestored = false;

        // hide cursor
        std::fprintf( mFile, "\033[?25l" );
        if( gLog.info_level() > 1 )
        {
            std::fprintf( mFile, "\n" );
        }

        this->flush();
    }

//------------------------------------------------------------------------------

    void
    Progressbar::step( const uint & aProgress )
    {
        // remember progress
        mProgress = aProgress ;

        if( gLog.info_level() > 1 && mNumSteps > 0 )
        {
            // width of the bar in characters. computed in index_t and
            // clamped, so that an overshooting caller can neither overflow
            // the product nor run past the end of the field
            index_t tStep = ( ( index_t ) aProgress * mWidth ) / mNumSteps ;

            if( tStep > mWidth ) tStep = mWidth ;

            if( tStep > mStep )
            {
                this->draw( ( uint ) tStep, aProgress );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Progressbar::step()
    {
        this->step( mProgress + 1 );
    }

//------------------------------------------------------------------------------

    void
    Progressbar::finish()
    {
        if( gLog.info_level() > 1 && mNumSteps > 0 )
        {
            // force the closing frame, unless the last step already drew it.
            // a full-width bar can only come from aProgress >= mNumSteps, so
            // mStep == mWidth means 100.0 % is already on screen
            if( mStep < mWidth )
            {
                this->draw( ( uint ) mWidth, ( uint ) mNumSteps );
            }

            std::fprintf( mFile, "\n" );
        }

        this->show_cursor();
    }

//------------------------------------------------------------------------------

    void
    Progressbar::show_cursor()
    {
        if( ! mCursorRestored )
        {
            mCursorRestored = true;
            std::fprintf( mFile, "\033[?25h" );
            this->flush();
        }
    }

//------------------------------------------------------------------------------

    void
    Progressbar::draw( const uint aStep, const uint aProgress )
    {
        mStep = aStep ;

        std::fprintf( mFile, "\r" );

        std::fprintf( mFile, "    " );

        for ( uint k = 1; k < aStep; ++k )
        {
            std::fprintf( mFile, "=" );
        }
        std::fprintf( mFile, ">" );

        for( uint k=aStep+1; k<mWidth; ++k )
        {
            std::fprintf( mFile, " " );
        }
        std::fprintf( mFile, "%-#4.1f %% ", aProgress * 100. / ( ( float ) mNumSteps ) );

        this->flush();
    }

//------------------------------------------------------------------------------

    void
    Progressbar::flush()
    {
        // The bar redraws with '\r' and never emits a newline, so a
        // line-buffered stream would hold the frame back. Under mpirun stdout
        // is a pipe rather than a terminal and is therefore FULLY buffered:
        // without this flush the whole run's output sits in the 4 KiB block
        // buffer and reaches the user in one burst at exit -- which is the
        // "nothing, then a finished bar" symptom. Flushing once per redraw
        // pushes the complete frame out in a single write(), so the I/O
        // forwarder cannot tear it either.
        std::fflush( mFile );
    }

//------------------------------------------------------------------------------
}
