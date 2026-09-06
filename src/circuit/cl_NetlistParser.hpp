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

#ifndef BELFEM_CL_NETLISTPARSER_HPP
#define BELFEM_CL_NETLISTPARSER_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        /**
         * one element card ( R/C/L/V/I/D ), lexed but not interpreted:
         * all payloads stay strings; numeric conversion and semantic
         * validation are the netlist factory's job
         */
        struct NetlistElement
        {
            //! element letter, case-folded: 'r' 'c' 'l' 'v' 'i' 'd'
            char                    mType = '?';

            //! full instance name, case-folded ( "is", "l1" ) -- becomes
            //! the component label
            string                  mName;

            //! the two node names, case-folded, as written ( "0", "gnd",
            //! "in" ); the ground remap happens in the factory
            Cell< string >          mNodes;

            //! positional tokens after the nodes, case-folded, parentheses
            //! and commas normalized to spaces ( "sin", "0", "500", "10" )
            Cell< string >          mValues;

            //! key=value tokens ( "r" -> "5k", "ic" -> "0" ), case-folded
            Map< string, string >   mKwargs;

            //! 1-based physical line of the card's first line
            index_t                 mLine = 0;

            //! the joined logical card as written, for error reporting
            string                  mCard;
        };

//-----------------------------------------------------------------------------

        /**
         * one .model card ( ".model dmod d( is=1e-14 )" )
         */
        struct NetlistModel
        {
            string                  mName;    // case-folded model name
            string                  mType;    // case-folded type token ( "d" )
            Map< string, string >   mKwargs;
            index_t                 mLine = 0;
            string                  mCard;
        };

//-----------------------------------------------------------------------------

        /**
         * one supported control card ( only ".tran" today; ".end" sets a
         * flag instead of a record )
         */
        struct NetlistControl
        {
            string                  mKeyword; // without the dot ( "tran" )
            Cell< string >          mArgs;
            index_t                 mLine = 0;
            string                  mCard;
        };

//-----------------------------------------------------------------------------

        /**
         * one "* belfem:" extension directive ( plan §5, Option A ):
         * superconductor, timed switch, per-instance BDF order
         */
        struct NetlistDirective
        {
            string                  mKind;    // case-folded ( "order" )
            Cell< string >          mArgs;    // positional tokens
            Map< string, string >   mKwargs;
            index_t                 mLine = 0;
            string                  mCard;
        };

//-----------------------------------------------------------------------------

        /**
         * Lexes an ngspice netlist into the intermediate representation the
         * NgspiceCircuitFactory consumes ( todo/ngspice_parser_plan.md §9.1 ).
         *
         * Handled here: the mandatory title line ( always consumed, never
         * parsed ); full-line '*' comments; "* belfem:" directives ( single
         * line, recognized before the comment strip, whitespace tolerated
         * before the colon ); end-of-line comments ';', '$' and
         * token-leading '//' ( PSPICE compatibility decks, where '$' is an
         * ordinary character, are not supported ); leading '+'
         * continuations; '(' ')' ',' normalized to spaces and whitespace
         * around '=' removed before tokenizing; case-folding of every
         * identifier; ".end" ( exactly that token -- trailing tokens
         * hard-error ) after which remaining lines are ignored; the ordered
         * first-appearance node-name list ( restart contract O9 ),
         * including nodes introduced only by a directive's n+/n- keys.
         *
         * Two deliberate BELFEM extensions beyond the manual's strict
         * wording ( the manual both demands that continuations "immediately
         * follow" and declares empty lines ignored and comment lines legal
         * anywhere ): comment and blank lines may sit between a card and
         * its '+' continuation, and EOF terminates a deck that has no
         * ".end". Both accept strictly more than stock ngspice and never
         * reinterpret a deck ngspice accepts. A "* belfem:" DIRECTIVE is a
         * statement, not a comment: it completes the pending card ( so
         * node collection follows line order, O9 ) and ends the
         * continuation chain.
         *
         * Hard errors ( BELFEM_ERROR, each naming source, line and card ):
         * unsupported element letters -- controlled switches S/W ( use the
         * "* belfem: switch" directive ), subcircuit instances X, and every
         * other letter beyond R/C/L/V/I/D; every control card except
         * ".model", ".tran" and ".end" ( .subckt/.control/.ac/.dc/.op/
         * .param/.ic/... -- v1 refuses rather than silently altering the
         * deck, plan §12.2 ); a continuation with no card to continue; an
         * element card with fewer than two nodes; a deck with no cards.
         *
         * Parameter-level policy ( e.g. rejecting "ic=" on a supported
         * card ) is deliberately NOT enforced here -- the factory owns
         * card semantics.
         */
        class NetlistParser
        {
            string                    mSource;  // path or "<memory>", for errors
            string                    mTitle;
            Cell< NetlistElement >    mElements;
            Cell< NetlistModel >      mModels;
            Cell< NetlistControl >    mControls;
            Cell< NetlistDirective >  mDirectives;
            Cell< string >            mNodeNames; // first appearance, folded
            bool                      mHasEnd = false;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            /**
             * lex a netlist file from disk
             */
            NetlistParser( const string & aPath );

            /**
             * lex an in-memory netlist ( unit tests ); aLines are the
             * physical lines including the title line
             */
            NetlistParser( const Cell< string > & aLines,
                           const string & aSource = "<memory>" );

            ~NetlistParser() = default;

//-----------------------------------------------------------------------------

            const string &
            title() const
            {
                return mTitle;
            }

            /**
             * the path or "<memory>" label errors are reported under
             */
            const string &
            source() const
            {
                return mSource;
            }

            const Cell< NetlistElement > &
            elements() const
            {
                return mElements;
            }

            const Cell< NetlistModel > &
            models() const
            {
                return mModels;
            }

            const Cell< NetlistControl > &
            controls() const
            {
                return mControls;
            }

            const Cell< NetlistDirective > &
            directives() const
            {
                return mDirectives;
            }

            /**
             * distinct node names in order of first appearance -- the
             * factory's packing order ( ground moves to the last index
             * there, not here )
             */
            const Cell< string > &
            node_names() const
            {
                return mNodeNames;
            }

            /**
             * true if the deck carried an explicit .end
             */
            bool
            has_end() const
            {
                return mHasEnd;
            }

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            parse( const Cell< string > & aLines );

            //! strip ';', '$' and token-leading '//' comments from one line
            string
            strip_eol_comment( const string & aLine ) const;

            //! classify and store one joined logical card
            void
            classify_card( const string & aCard,
                           const index_t aLine );

            void
            classify_control_card( const Cell< string > & aTokens,
                                   const string & aCard,
                                   const index_t aLine );

            void
            classify_element_card( const Cell< string > & aTokens,
                                   const string & aCard,
                                   const index_t aLine );

            void
            parse_directive( const string & aPayload,
                             const index_t aLine );

            //! '('/')'/',' to spaces, drop whitespace around '=', fold case,
            //! then whitespace-tokenize
            void
            tokenize( const string & aCard,
                      Cell< string > & aTokens ) const;

            //! split one "key=value" token into aKwargs; hard-errors on a
            //! missing key or value, a second '=' in the value, and a
            //! duplicate key ( last-wins would be silent )
            void
            store_kwarg( const string & aToken,
                         const string & aCard,
                         const index_t aLine,
                         Map< string, string > & aKwargs ) const;

            //! register a node name if it has not appeared yet
            void
            collect_node( const string & aName );

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_NETLISTPARSER_HPP
