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

#ifndef BELFEM_FN_MESH_CONFIG_TAG_HPP
#define BELFEM_FN_MESH_CONFIG_TAG_HPP

#include <cstdio>
#include <algorithm>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
// cl_Input_Section.hpp declares Vector< id_t > without including it, so this
// header must come first for this file to be includable on its own
#include "cl_Vector.hpp"
#include "stringtools.hpp"
#include "cl_InputFile.hpp"
#include "fn_FEM_ghost_switch.hpp"
#include "cl_Input_Section.hpp"

namespace belfem
{
    namespace fem
    {
    namespace maxwell
    {
//------------------------------------------------------------------------------

        /**
         * The mesh CONFIGURATION tag: a fingerprint of the settings that decide
         * how a .msh is turned into an enriched mesh -- cuts, thin-shell layers,
         * edge-coating walls, periodicity.
         *
         * It exists because the .bfm reuse test compares only Mesh::checksum(),
         * which is the identity of the BASE mesh ( node coordinates and element
         * connectivity ). The two answer different questions and both are needed:
         * the checksum catches a changed .msh under an unchanged deck, this tag
         * catches a changed deck under an unchanged .msh.
         *
         * Deliberately NOT covered, so the user can retune a run against a cached
         * mesh: solver, timestepping, tolerances, output, and boundary-condition
         * amplitudes. Terminal ids ARE covered, because they feed cut construction,
         * and so is the thin-shell ghost switch ( nonlinear magnetic {
         * nitsche ghost penalty { eta } } ), because it decides whether the
         * layer interfaces carry duplicate dofs -- a discretization, not a
         * tuning.
         *
         * Known gap, accepted: a change INSIDE a material definition ( one that
         * adds or removes rho ) can alter thin-shell node duplication while every
         * layer label stays the same. Material labels are covered, definitions are
         * not.
         *
         * The tag must be reproducible on another machine, because a .bfm travels
         * with its deck and mesh. Hence FNV-1a over fixed-precision TEXT rather
         * than std::hash over raw values: std::hash is unspecified across standard
         * libraries, and raw IEEE bytes would carry endianness plus last-ulp noise
         * ( 100 um and 0.1 mm must give the same tag ).
         */

//------------------------------------------------------------------------------

        /**
         * FNV-1a, 64 bit. A defined byte algorithm: identical text yields an
         * identical value on every compiler, library and platform.
         */
        inline uint64_t
        mesh_config_hash( const string & aText )
        {
            uint64_t tHash = 0xcbf29ce484222325ULL ;

            for ( const char tChar : aText )
            {
                tHash ^= static_cast< uint64_t >(
                        static_cast< unsigned char >( tChar ) );
                tHash *= 0x100000001b3ULL ;
            }

            return tHash ;
        }

//------------------------------------------------------------------------------

        namespace config_tag
        {
            /**
             * True if the WHOLE word is a single decimal number.
             *
             * The whole word matters: id lists are written "1,3,5,7,9,11" and
             * ranges "1:3", both of which begin with a digit. Testing only the
             * first character would send them through to_real(), which stops at
             * the first separator -- so "1,3,5,7,9,11" and "1,3,5,7,9,13" would
             * both collapse to 1 and produce the SAME tag. Anything carrying a
             * separator must stay a string and be compared verbatim.
             */
            inline bool
            looks_numeric( const string & aWord )
            {
                if ( aWord.size() == 0 ) return false ;

                bool tHasDigit = false ;

                for ( const char tChar : aWord )
                {
                    if ( tChar >= '0' && tChar <= '9' )
                    {
                        tHasDigit = true ;
                    }
                    else if ( tChar != '+' && tChar != '-' && tChar != '.'
                           && tChar != 'e' && tChar != 'E' )
                    {
                        return false ;
                    }
                }

                return tHasDigit ;
            }

//------------------------------------------------------------------------------

            /**
             * Canonical form of one input value. A quantity with a unit becomes
             * its SI magnitude, so 100 um and 0.1 mm compare equal; anything else
             * is lower-cased and whitespace-collapsed. Twelve significant digits
             * absorb unit-conversion round-off while still separating any value a
             * user would type.
             */
            inline string
            value_to_canonical( const string & aString )
            {
                Cell< string > tWords = string_to_words( aString );

                if ( tWords.size() == 0 ) return "" ;

                // a number, optionally followed by a unit
                if ( tWords.size() <= 2 && looks_numeric( tWords( 0 ) ) )
                {
                    real tValue = to_real( tWords( 0 ) );

                    if ( tWords.size() == 2 )
                    {
                        // scale into SI; unit_to_si returns 1.0 for anything it
                        // does not recognize, which is the identity we want here
                        tValue *= unit_to_si( tWords( 1 ) ).first ;
                    }

                    char tBuffer[ 32 ];
                    std::snprintf( tBuffer, sizeof( tBuffer ), "%.12g", tValue );
                    return string( tBuffer );
                }

                // not a quantity: join the words so spacing cannot matter
                string tResult = string_to_lower( tWords( 0 ) );
                for ( index_t k = 1 ; k < tWords.size() ; ++k )
                {
                    tResult += " " + string_to_lower( tWords( k ) );
                }

                // fold the boolean spellings onto one form, so "edge coating :
                // on" and "... : true" are one configuration rather than two
                // tags. The NUMERIC spellings are deliberately left alone: "1"
                // is also a perfectly good entity id, and mapping it to "true"
                // would let an id collide with a flag
                if ( tResult == "on" || tResult == "true" || tResult == "yes" )
                {
                    return "true" ;
                }
                if ( tResult == "off" || tResult == "false" || tResult == "no" )
                {
                    return "false" ;
                }

                return tResult ;
            }

//------------------------------------------------------------------------------

            /**
             * Emit every key of a section, then recurse into its subsections.
             * Enumerating keys rather than naming them keeps a newly added
             * topology key inside the tag automatically.
             */
            inline void
            collect_section(
                    const input::Section * aSection,
                    const string         & aPath,
                    Cell< string >       & aLines,
                    const index_t          aIndex = 0 )
            {
                // see collect_terminals: the index keeps same-type siblings apart
                const string tPath = aPath + "/" + string_to_lower( aSection->type() )
                                   + ":" + string_to_lower( aSection->label() )
                                   + "[" + std::to_string( aIndex ) + "]" ;

                for ( index_t k = 0 ; k < aSection->num_keys() ; ++k )
                {
                    const string & tKey = aSection->key( k );

                    aLines.push( tPath + "." + string_to_lower( tKey ) + " = "
                                 + value_to_canonical( aSection->get_string( tKey ) ) );
                }

                for ( index_t k = 0 ; k < aSection->num_sections() ; ++k )
                {
                    collect_section( aSection->section( k ), tPath, aLines, k );
                }
            }

//------------------------------------------------------------------------------

            /**
             * Collect the current-injection domains from a section tree.
             *
             * Recursive and applied to BOTH `boundary conditions` and `circuit`,
             * because the ids that feed cut construction hide in three places:
             *   - `boundary conditions/current { input terminals ... }`
             *   - the nested `boundary conditions/maxwell/...` form the factory
             *     also accepts
             *   - `circuit/.../{ input curves ... }`, whose terminal pairs become
             *     CircuitVoltage / CircuitCurrent BC domains and end up in the
             *     same terminal list
             * A deck can carry cut-defining ids in the circuit section and none
             * under `boundary conditions` at all ( examples/tapestack_circuit ).
             *
             * Only terminal-ish keys are taken. Amplitude, period and waveform
             * stay out so the excitation can be retuned against a cached mesh.
             */
            inline void
            collect_terminals(
                    const input::Section * aSection,
                    const string         & aPath,
                    Cell< string >       & aLines,
                    const index_t          aIndex = 0 )
            {
                // the index disambiguates same-type siblings: two unlabelled
                // `terminal pair` sections would otherwise share a path, and
                // swapping ids BETWEEN them would leave the line set - and the
                // tag - unchanged. Order sensitivity is the safe trade here:
                // reordering costs a rebuild, a collision would cost a stale mesh
                const string tPath = aPath + "/" + string_to_lower( aSection->type() )
                                   + ":" + string_to_lower( aSection->label() )
                                   + "[" + std::to_string( aIndex ) + "]" ;

                for ( index_t k = 0 ; k < aSection->num_keys() ; ++k )
                {
                    const string tKey = string_to_lower( aSection->key( k ) );

                    // both spellings are in use: helix writes "input terminals",
                    // corc and the circuit decks write "input curves"
                    if ( tKey.find( "terminal" ) == string::npos
                      && tKey.find( "curve" )    == string::npos ) continue ;

                    aLines.push( tPath + "." + tKey + " = "
                                 + value_to_canonical(
                                         aSection->get_string( aSection->key( k ) ) ) );
                }

                for ( index_t k = 0 ; k < aSection->num_sections() ; ++k )
                {
                    collect_terminals( aSection->section( k ), tPath, aLines, k );
                }
            }

//------------------------------------------------------------------------------

            /**
             * The layer stack of one thin shell. Parsed from the raw buffer, not
             * from the key API: read_thin_shell_data() word-splits these lines,
             * so a key-based walk would miss the thicknesses entirely.
             */
            inline void
            collect_layers(
                    const input::Section * aSection,
                    Cell< string >       & aLines )
            {
                const string tPath = "layers:" + string_to_lower( aSection->label() );

                const Cell< string > & tBuffer = aSection->buffer();

                index_t tCount = 0 ;

                for ( index_t k = aSection->start() ; k < aSection->end() ; ++k )
                {
                    Cell< string > tWords = string_to_words( tBuffer( k ) );

                    // material <sep> thickness unit
                    if ( tWords.size() < 4 ) continue ;

                    real tThickness = to_real( tWords( 2 ) )
                                    * unit_to_si( tWords( 3 ) ).first ;

                    char tBuf[ 32 ];
                    std::snprintf( tBuf, sizeof( tBuf ), "%.12g", tThickness );

                    aLines.push( tPath + ".layer[" + std::to_string( tCount ) + "].material = "
                                 + string_to_lower( tWords( 0 ) ) );
                    aLines.push( tPath + ".layer[" + std::to_string( tCount ) + "].thickness = "
                                 + string( tBuf ) );
                    ++tCount ;
                }

                // the layer COUNT is itself mesh-defining: it sets how many layer
                // blocks are built
                aLines.push( tPath + ".count = " + std::to_string( tCount ) );
            }
        }

//------------------------------------------------------------------------------

        /**
         * The canonical text behind the tag. Stored in the .bfm next to the tag
         * so a mismatch can show WHICH setting changed instead of only reporting
         * that two numbers differ.
         *
         * Lines are sorted, so neither the order of sections in the deck nor the
         * iteration order of any container can change the result.
         */
        inline string
        mesh_config_text( const InputFile & aInputFile )
        {
            Cell< string > tLines ;

            // -- topology: domain types, thin shells and their sidesets, edge
            //    coating, periodic source/target, curves
            if ( aInputFile.section_exists( "topology" ) )
            {
                config_tag::collect_section(
                        aInputFile.section( "topology" ), "", tLines );
            }

            // -- homology: the cut algorithm changes the cuts baked into the file
            if ( aInputFile.section_exists( "homology" ) )
            {
                const input::Section * tSection = aInputFile.section( "homology" );

                if ( tSection->key_exists( "algorithm" ) )
                {
                    tLines.push( "homology.algorithm = " + config_tag::value_to_canonical(
                            tSection->get_string( "algorithm" ) ) );
                }
            }

            // -- terminal ids that feed the cuts. Both the boundary-condition
            //    tree and the circuit tree, recursively: a deck can put its
            //    cut-defining ids in either ( examples/tapestack_circuit has them ONLY
            //    under `circuit` ). Amplitudes and waveforms stay out
            if ( aInputFile.section_exists( "boundary conditions" ) )
            {
                config_tag::collect_terminals(
                        aInputFile.section( "boundary conditions" ), "", tLines );
            }

            if ( aInputFile.section_exists( "circuit" ) )
            {
                config_tag::collect_terminals(
                        aInputFile.section( "circuit" ), "", tLines );
            }

            // -- layer stacks, matched to a thin shell by label
            for ( index_t k = 0 ; k < aInputFile.num_sections() ; ++k )
            {
                const input::Section * tSection = aInputFile.section( k );

                if ( string_to_lower( tSection->type() ) == "layers" )
                {
                    config_tag::collect_layers( tSection, tLines );
                }
            }

            // -- the ONE solver-section setting that changes the discretization:
            //    whether layer interfaces get duplicate dofs and ghost facets.
            //    Added 2026-09-01; every cache built before then misses once
            //    and is rebuilt, which is the point ( a stale cache silently
            //    kept the other layout )
            tLines.push( string( "thinshell.ghost = " )
                         + ( fem::ghost_facets_requested( aInputFile ) ? "on" : "off" ) );

            // sort so deck order cannot change the tag
            std::sort( tLines.vector_data().begin(), tLines.vector_data().end() );

            string tText ;
            for ( const string & tLine : tLines )
            {
                tText += tLine + "\n" ;
            }

            return tText ;
        }

//------------------------------------------------------------------------------

        /**
         * Convenience: the tag of an input file.
         */
        inline uint64_t
        mesh_config_tag( const InputFile & aInputFile )
        {
            return mesh_config_hash( mesh_config_text( aInputFile ) );
        }

//------------------------------------------------------------------------------
    }
    }
}
#endif //BELFEM_FN_MESH_CONFIG_TAG_HPP
