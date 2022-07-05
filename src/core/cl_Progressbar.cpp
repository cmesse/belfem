//
// Created by Christian Messe on 03.11.19.
//

#include <cstdio>
#include <iostream>

#include "cl_Progressbar.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    Progressbar::Progressbar( const index_t aNumSteps, FILE * aFile ) :
        mNumSteps( aNumSteps ),
        mFile( aFile )
    {
        this->reset();
    }

//------------------------------------------------------------------------------
    void
    Progressbar::reset()
    {
        mStep = 0;

        // hide cursor
        std::fprintf( mFile, "\033[?25l" );
        if( gLog.info_level() > 1 )
        {
            std::fprintf( mFile, "\n" );
        }
    }

//------------------------------------------------------------------------------

    void
    Progressbar::step( const uint & aProgress )
    {
        // remember progress
        mProgress = aProgress ;

        if( gLog.info_level() > 1 )
        {
            uint tStep = ( aProgress * 71 ) / mNumSteps;

            if( tStep > mStep )
            {
                for ( uint k = 0; k < 80; ++k )
                {
                    std::fprintf( mFile, "\b" );
                }

                std::fprintf( mFile, " " );

                for ( uint k = 1; k < tStep; ++k )
                {
                    std::fprintf( mFile, "=" );
                }
                std::fprintf( mFile, ">" );

                for( uint k=tStep+1; k<72; ++k )
                {
                    std::fprintf( mFile, " " );
                }
                std::fprintf( mFile, "%-#4.1f %% ", aProgress * 100 / ( ( float) mNumSteps ) );
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
        if( gLog.info_level() > 1 )
        {
            this->step( mNumSteps );
            std::fprintf( mFile, "\n" );
        }

        // show cursor
        std::fprintf( mFile, "\033[?25h" );
    }
//------------------------------------------------------------------------------
}