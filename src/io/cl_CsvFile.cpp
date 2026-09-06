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

#include <algorithm>
#include "cl_CsvFile.hpp"

namespace belfem
{
    CsvFile::CsvFile( const string & aPath, const char aDelimiter ) :
        Ascii( aPath, FileMode::OPEN_RDONLY )
    {
        if ( mBuffer.size() == 0 ) return;

        index_t m = mBuffer.size();
        index_t n = std::count( mBuffer(0).begin(), mBuffer(0).end(), aDelimiter ) + 1 ;

        mData.set_size( m, n );

        for ( index_t i=0; i<m; ++i )
        {
            string tLine = clean_string( mBuffer( i ) );
            Cell< string > tWords = string_to_words( tLine, aDelimiter );

            for ( index_t j=0; j<n; ++j )
            {
                mData( i, j ) = std::stod( tWords( j ) );
            }
        }
    }
    CsvFile::~CsvFile()
    {

    }

}