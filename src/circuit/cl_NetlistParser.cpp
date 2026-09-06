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

#include <cctype>

#include "cl_NetlistParser.hpp"
#include "assert.hpp"
#include "stringtools.hpp"
#include "cl_Ascii.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
//  local helpers
//-----------------------------------------------------------------------------

        namespace
        {
            bool
            is_space( const char aChar )
            {
                return std::isspace( static_cast< unsigned char >( aChar ) ) != 0;
            }

            string
            trim( const string & aLine )
            {
                size_t tBegin = 0;
                size_t tEnd   = aLine.length();
                while ( tBegin < tEnd && is_space( aLine[ tBegin ] ) )
                {
                    ++tBegin;
                }
                while ( tEnd > tBegin && is_space( aLine[ tEnd - 1 ] ) )
                {
                    --tEnd;
                }
                return aLine.substr( tBegin, tEnd - tBegin );
            }

            //! true if the line is a "* belfem:" directive; on success
            //! aPayload receives the text after the marker
            bool
            is_directive( const string & aLine, string & aPayload )
            {
                size_t tPos = 0;
                const size_t tLength = aLine.length();
                while ( tPos < tLength && is_space( aLine[ tPos ] ) )
                {
                    ++tPos;
                }
                if ( tPos >= tLength || aLine[ tPos ] != '*' )
                {
                    return false;
                }
                ++tPos;
                while ( tPos < tLength && is_space( aLine[ tPos ] ) )
                {
                    ++tPos;
                }
                const char * tMarker = "belfem";
                for ( size_t k = 0; tMarker[ k ] != '\0'; ++k )
                {
                    if ( tPos + k >= tLength ||
                         std::tolower( static_cast< unsigned char >( aLine[ tPos + k ] ) )
                         != tMarker[ k ] )
                    {
                        return false;
                    }
                }
                tPos += 6;
                // tolerate whitespace before the colon ( "* belfem : order" )
                while ( tPos < tLength && is_space( aLine[ tPos ] ) )
                {
                    ++tPos;
                }
                if ( tPos >= tLength || aLine[ tPos ] != ':' )
                {
                    return false;
                }
                aPayload = aLine.substr( tPos + 1 );
                return true;
            }
        }

//-----------------------------------------------------------------------------
//  construction
//-----------------------------------------------------------------------------

        NetlistParser::NetlistParser( const string & aPath ) :
            mSource( aPath )
        {
            Ascii tFile( aPath, FileMode::OPEN_RDONLY );
            const index_t tNumLines = tFile.length();
            Cell< string > tLines( tNumLines, "" );
            for ( index_t k = 0; k < tNumLines; ++k )
            {
                tLines( k ) = tFile.line( k );
            }
            this->parse( tLines );
        }

//-----------------------------------------------------------------------------

        NetlistParser::NetlistParser( const Cell< string > & aLines,
                                      const string & aSource ) :
            mSource( aSource )
        {
            this->parse( aLines );
        }

//-----------------------------------------------------------------------------
//  the line loop
//-----------------------------------------------------------------------------

        void
        NetlistParser::parse( const Cell< string > & aLines )
        {
            BELFEM_ERROR( aLines.size() > 0,
                          "netlist %s is empty ( not even a title line )",
                          mSource.c_str() );

            // the first line is always the title -- consumed, never parsed
            mTitle = trim( aLines( 0 ) );

            string  tPending;          // logical card under assembly
            index_t tPendingLine = 0;  // 1-based line where it started

            // classify-and-clear, shared by the loop and the epilogue
            auto tFlush = [ this, &tPending, &tPendingLine ]()
            {
                if ( !tPending.empty() )
                {
                    this->classify_card( tPending, tPendingLine );
                    tPending.clear();
                }
            };

            for ( index_t k = 1; k < aLines.size(); ++k )
            {
                const index_t tLineNumber = k + 1;  // 1-based, incl. title

                // drop a trailing carriage return ( CRLF files )
                string tLine = aLines( k );
                if ( !tLine.empty() && tLine.back() == '\r' )
                {
                    tLine.pop_back();
                }

                // extension directives hide behind the comment marker and
                // must be recognized before the comment skip. A directive is
                // a STATEMENT, not a comment: it completes the pending card
                // first, so node collection follows line order ( O9 ), and
                // it ends any continuation chain
                string tPayload;
                if ( is_directive( tLine, tPayload ) )
                {
                    tFlush();
                    this->parse_directive(
                            trim( this->strip_eol_comment( tPayload ) ),
                            tLineNumber );
                    continue;
                }

                // full-line comment; does not break a continuation chain
                const string tTrimmedRaw = trim( tLine );
                if ( !tTrimmedRaw.empty() && tTrimmedRaw[ 0 ] == '*' )
                {
                    continue;
                }

                const string tTrimmed =
                        trim( this->strip_eol_comment( tLine ) );

                // blank line; does not break a continuation chain
                if ( tTrimmed.empty() )
                {
                    continue;
                }

                // leading '+' continues the pending card
                if ( tTrimmed[ 0 ] == '+' )
                {
                    BELFEM_ERROR( !tPending.empty(),
                                  "%s:%lu: continuation line with no card to continue: '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) tLineNumber,
                                  tTrimmed.c_str() );
                    tPending += " ";
                    tPending += tTrimmed.substr( 1 );
                    continue;
                }

                // a new card begins: the previous one is complete
                tFlush();

                // ".end" terminates the deck; anything after it is ignored.
                // The whitespace test at [ 4 ] keeps .endc/.ends/.endif on
                // the unsupported-control path -- do not "simplify" it away
                const string tFolded = string_to_lower( tTrimmed );
                if ( tFolded == ".end" ||
                     ( tFolded.length() > 4 &&
                       tFolded.compare( 0, 4, ".end" ) == 0 &&
                       is_space( tFolded[ 4 ] ) ) )
                {
                    // the manual's form is exactly ".end" -- trailing junk
                    // silently ending a deck would hide a malformed card
                    BELFEM_ERROR( tFolded == ".end",
                                  "%s:%lu: unexpected tokens after .end: '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) tLineNumber,
                                  tTrimmed.c_str() );
                    mHasEnd = true;
                    break;
                }

                tPending     = tTrimmed;
                tPendingLine = tLineNumber;
            }

            // EOF also terminates the deck ( .end is optional, documented )
            tFlush();

            BELFEM_ERROR( mElements.size() > 0,
                          "netlist %s contains no element cards",
                          mSource.c_str() );
        }

//-----------------------------------------------------------------------------
//  comment stripping
//-----------------------------------------------------------------------------

        string
        NetlistParser::strip_eol_comment( const string & aLine ) const
        {
            const size_t tLength = aLine.length();
            for ( size_t k = 0; k < tLength; ++k )
            {
                const char tChar = aLine[ k ];
                if ( tChar == ';' || tChar == '$' )
                {
                    return aLine.substr( 0, k );
                }
                // '//' only at the start of a token, so "1e-//2" cannot
                // lose its tail to a false comment
                if ( tChar == '/' && k + 1 < tLength && aLine[ k + 1 ] == '/' &&
                     ( k == 0 || is_space( aLine[ k - 1 ] ) ) )
                {
                    return aLine.substr( 0, k );
                }
            }
            return aLine;
        }

//-----------------------------------------------------------------------------
//  tokenizing
//-----------------------------------------------------------------------------

        void
        NetlistParser::tokenize( const string & aCard,
                                 Cell< string > & aTokens ) const
        {
            // pass 1: fold case, turn '(' ')' ',' into spaces
            string tNorm;
            tNorm.reserve( aCard.length() );
            for ( const char tChar : aCard )
            {
                if ( tChar == '(' || tChar == ')' || tChar == ',' )
                {
                    tNorm += ' ';
                }
                else
                {
                    tNorm += static_cast< char >(
                            std::tolower( static_cast< unsigned char >( tChar ) ) );
                }
            }

            // pass 2: drop whitespace around '=' so "r = 5" and "r=5"
            // tokenize identically
            string tClean;
            tClean.reserve( tNorm.length() );
            for ( size_t k = 0; k < tNorm.length(); ++k )
            {
                if ( tNorm[ k ] == '=' )
                {
                    while ( !tClean.empty() && is_space( tClean.back() ) )
                    {
                        tClean.pop_back();
                    }
                    tClean += '=';
                    while ( k + 1 < tNorm.length() && is_space( tNorm[ k + 1 ] ) )
                    {
                        ++k;
                    }
                }
                else
                {
                    tClean += tNorm[ k ];
                }
            }

            // pass 3: split on whitespace
            aTokens.clear();
            size_t tPos = 0;
            while ( tPos < tClean.length() )
            {
                while ( tPos < tClean.length() && is_space( tClean[ tPos ] ) )
                {
                    ++tPos;
                }
                size_t tBegin = tPos;
                while ( tPos < tClean.length() && !is_space( tClean[ tPos ] ) )
                {
                    ++tPos;
                }
                if ( tPos > tBegin )
                {
                    aTokens.push( tClean.substr( tBegin, tPos - tBegin ) );
                }
            }
        }

//-----------------------------------------------------------------------------
//  classification
//-----------------------------------------------------------------------------

        void
        NetlistParser::classify_card( const string & aCard,
                                      const index_t aLine )
        {
            Cell< string > tTokens;
            this->tokenize( aCard, tTokens );

            // a card of only parentheses/commas normalizes to nothing --
            // reject it here so the token accesses below stay in bounds
            BELFEM_ERROR( tTokens.size() > 0,
                          "%s:%lu: card contains no tokens: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aCard.c_str() );

            const char tFirst = tTokens( 0 )[ 0 ];

            if ( tFirst == '.' )
            {
                this->classify_control_card( tTokens, aCard, aLine );
            }
            else if ( std::isalpha( static_cast< unsigned char >( tFirst ) ) )
            {
                this->classify_element_card( tTokens, aCard, aLine );
            }
            else
            {
                BELFEM_ERROR( false,
                              "%s:%lu: card must start with an element letter or '.': '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aLine,
                              aCard.c_str() );
            }
        }

//-----------------------------------------------------------------------------

        void
        NetlistParser::classify_control_card( const Cell< string > & aTokens,
                                              const string & aCard,
                                              const index_t aLine )
        {
            const string tKeyword = aTokens( 0 ).substr( 1 );

            if ( tKeyword == "model" )
            {
                BELFEM_ERROR( aTokens.size() >= 3,
                              "%s:%lu: .model needs a name and a type: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aLine,
                              aCard.c_str() );

                NetlistModel tModel;
                tModel.mName = aTokens( 1 );
                tModel.mType = aTokens( 2 );
                tModel.mLine = aLine;
                tModel.mCard = aCard;
                for ( index_t k = 3; k < aTokens.size(); ++k )
                {
                    BELFEM_ERROR( aTokens( k ).find( '=' ) != string::npos,
                                  "%s:%lu: .model parameters must be key=value, got '%s' in '%s'",
                                  mSource.c_str(),
                                  ( long unsigned int ) aLine,
                                  aTokens( k ).c_str(),
                                  aCard.c_str() );
                    this->store_kwarg( aTokens( k ), aCard, aLine, tModel.mKwargs );
                }
                mModels.push( tModel );
            }
            else if ( tKeyword == "tran" )
            {
                NetlistControl tControl;
                tControl.mKeyword = tKeyword;
                tControl.mLine    = aLine;
                tControl.mCard    = aCard;
                for ( index_t k = 1; k < aTokens.size(); ++k )
                {
                    tControl.mArgs.push( aTokens( k ) );
                }
                mControls.push( tControl );
            }
            else
            {
                // .subckt/.ends/.control/.ac/.dc/.op/.param/.ic/.include/...
                // v1 refuses rather than silently running an altered deck
                BELFEM_ERROR( false,
                              "%s:%lu: control card '.%s' is not supported in v1: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aLine,
                              tKeyword.c_str(),
                              aCard.c_str() );
            }
        }

//-----------------------------------------------------------------------------

        void
        NetlistParser::classify_element_card( const Cell< string > & aTokens,
                                              const string & aCard,
                                              const index_t aLine )
        {
            const char tType = aTokens( 0 )[ 0 ];  // already folded

            // the mismatches get their own messages
            BELFEM_ERROR( tType != 's' && tType != 'w',
                          "%s:%lu: controlled switches ( S/W cards ) are not supported -- "
                          "BELFEM's switch is timed; use the '* belfem: switch' directive: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aCard.c_str() );

            BELFEM_ERROR( tType != 'x',
                          "%s:%lu: subcircuit instances ( X cards / .subckt ) are not "
                          "supported in v1: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aCard.c_str() );

            BELFEM_ERROR( tType == 'r' || tType == 'c' || tType == 'l' ||
                          tType == 'v' || tType == 'i' || tType == 'd',
                          "%s:%lu: element type '%c' is not supported ( v1 handles "
                          "R, C, L, V, I, D ): '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          tType,
                          aCard.c_str() );

            BELFEM_ERROR( aTokens.size() >= 3,
                          "%s:%lu: element card needs a name and two nodes: '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aCard.c_str() );

            NetlistElement tElement;
            tElement.mType = tType;
            tElement.mName = aTokens( 0 );
            tElement.mLine = aLine;
            tElement.mCard = aCard;

            for ( index_t k = 1; k < 3; ++k )
            {
                const string & tNode = aTokens( k );
                BELFEM_ERROR( tNode.find( '=' ) == string::npos,
                              "%s:%lu: expected a node name, got '%s' -- element "
                              "cards need two nodes before any parameter: '%s'",
                              mSource.c_str(),
                              ( long unsigned int ) aLine,
                              tNode.c_str(),
                              aCard.c_str() );
                tElement.mNodes.push( tNode );
                this->collect_node( tNode );
            }

            for ( index_t k = 3; k < aTokens.size(); ++k )
            {
                const string & tToken = aTokens( k );
                if ( tToken.find( '=' ) == string::npos )
                {
                    tElement.mValues.push( tToken );
                }
                else
                {
                    this->store_kwarg( tToken, aCard, aLine, tElement.mKwargs );
                }
            }

            mElements.push( tElement );
        }

//-----------------------------------------------------------------------------
//  directives
//-----------------------------------------------------------------------------

        void
        NetlistParser::parse_directive( const string & aPayload,
                                        const index_t aLine )
        {
            Cell< string > tTokens;
            this->tokenize( aPayload, tTokens );

            BELFEM_ERROR( tTokens.size() > 0,
                          "%s:%lu: empty '* belfem:' directive",
                          mSource.c_str(),
                          ( long unsigned int ) aLine );

            NetlistDirective tDirective;
            tDirective.mKind = tTokens( 0 );
            tDirective.mLine = aLine;
            tDirective.mCard = aPayload;

            for ( index_t k = 1; k < tTokens.size(); ++k )
            {
                const string & tToken = tTokens( k );
                const size_t tEq = tToken.find( '=' );
                if ( tEq == string::npos )
                {
                    tDirective.mArgs.push( tToken );
                }
                else
                {
                    this->store_kwarg( tToken, aPayload, aLine, tDirective.mKwargs );

                    // directives can introduce nodes ( superconductor,
                    // switch ) -- they join the O9 first-appearance list
                    // in written order, like element-card nodes
                    const string tKey = tToken.substr( 0, tEq );
                    if ( tKey == "n+" || tKey == "n-" )
                    {
                        this->collect_node( tToken.substr( tEq + 1 ) );
                    }
                }
            }

            mDirectives.push( tDirective );
        }

//-----------------------------------------------------------------------------
//  kwarg storage
//-----------------------------------------------------------------------------

        void
        NetlistParser::store_kwarg( const string & aToken,
                                    const string & aCard,
                                    const index_t aLine,
                                    Map< string, string > & aKwargs ) const
        {
            const size_t tEq = aToken.find( '=' );

            BELFEM_ERROR( tEq != string::npos && tEq > 0 &&
                          tEq + 1 < aToken.length(),
                          "%s:%lu: malformed parameter '%s' in '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aToken.c_str(),
                          aCard.c_str() );

            const string tKey   = aToken.substr( 0, tEq );
            const string tValue = aToken.substr( tEq + 1 );

            // "r==5" would otherwise become r -> "=5"
            BELFEM_ERROR( tValue.find( '=' ) == string::npos,
                          "%s:%lu: malformed parameter '%s' in '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          aToken.c_str(),
                          aCard.c_str() );

            // last-wins would silently change the circuit
            BELFEM_ERROR( !aKwargs.key_exists( tKey ),
                          "%s:%lu: duplicate parameter '%s' in '%s'",
                          mSource.c_str(),
                          ( long unsigned int ) aLine,
                          tKey.c_str(),
                          aCard.c_str() );

            aKwargs[ tKey ] = tValue;
        }

//-----------------------------------------------------------------------------
//  node collection
//-----------------------------------------------------------------------------

        void
        NetlistParser::collect_node( const string & aName )
        {
            // linear scan keeps first-appearance order; node counts are tiny
            for ( index_t k = 0; k < mNodeNames.size(); ++k )
            {
                if ( mNodeNames( k ) == aName )
                {
                    return;
                }
            }
            mNodeNames.push( aName );
        }

//-----------------------------------------------------------------------------
    }
}
