//
// Created by christian on 8/28/26.
//

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
#include <cmath>
#include <iostream>

#include "assert.hpp"
#include "banner.hpp"
#include "cl_Communicator.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Arguments.hpp"
#include "cl_Mesh.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

namespace
{
//------------------------------------------------------------------------------

    void
    print_help()
    {
        std::cout << "Usage: db2exo <hdf5file>" << std::endl;
        std::cout << "  converts the tensor tables of a BELFEM database into an" << std::endl;
        std::cout << "  exodus file of the same name, one node field per table" << std::endl;
        std::cout << "  -h, --help                                print this help message" << std::endl;
        std::cout << "  -v, --verbose [level]                     set the info level" << std::endl;
    }

//------------------------------------------------------------------------------

    bool
    is_uint( const string & aString )
    {
        if ( aString.empty() )
        {
            return false ;
        }
        for ( const char tChar : aString )
        {
            if ( ! std::isdigit( static_cast< unsigned char >( tChar ) ) )
            {
                return false ;
            }
        }
        return true ;
    }

//------------------------------------------------------------------------------

    /**
     * the shared -v / --verbose flags are read by Arguments but not removed
     * from the list, so the file is the first argument that is neither a flag
     * nor the operand of one. Returns an empty string if there is none
     */
    string
    input_file_from_arguments( const Arguments & aArguments )
    {
        const Cell< string > & tArgs = aArguments.data() ;

        index_t tNumArgs = tArgs.size() ;

        for ( index_t k=1; k<tNumArgs; ++k )
        {
            const string & tArg = tArgs( k ) ;

            if ( tArg == "-v" || tArg == "--verbose" )
            {
                // an unsigned integer behind the bare flag is its operand
                if ( k+1 < tNumArgs && is_uint( tArgs( k+1 ) ) )
                {
                    ++k ;
                }
                continue ;
            }

            if ( ! tArg.empty() && tArg.front() == '-' )
            {
                continue ;
            }

            return tArg ;
        }

        return "" ;
    }

//------------------------------------------------------------------------------

    bool
    help_requested( const Arguments & aArguments )
    {
        for ( const string & tArg : aArguments.data() )
        {
            if ( tArg == "-h" || tArg == "--help" )
            {
                return true ;
            }
        }
        return false ;
    }

//------------------------------------------------------------------------------

    /**
     * replaces the extension of the file name, not of the path : a directory
     * in the path may carry a dot while the file name does not, and the
     * leading dot of a hidden file belongs to its name rather than marking an
     * extension
     */
    string
    exodus_path( const string & aPath )
    {
        std::size_t tSlash = aPath.find_last_of( "/" ) ;

        // first character of the file name
        std::size_t tBegin = ( tSlash == string::npos ) ? 0 : tSlash + 1 ;

        std::size_t tDot = aPath.find_last_of( "." ) ;

        // no extension : either no dot in the file name at all, or the only
        // dot is the one that makes the file hidden
        if ( tDot == string::npos || tDot <= tBegin )
        {
            return aPath + ".exo" ;
        }

        return aPath.substr( 0, tDot ) + ".exo" ;
    }

//------------------------------------------------------------------------------

    /**
     * a table is a group that carries the complete record written by
     * Database::save()
     */
    bool
    is_table( const hid_t aLoc )
    {
        return hdf5::dataset_exists( aLoc, "order" )
            && hdf5::dataset_exists( aLoc, "origin" )
            && hdf5::dataset_exists( aLoc, "step" )
            && hdf5::dataset_exists( aLoc, "points" )
            && hdf5::dataset_exists( aLoc, "values" ) ;
    }

//------------------------------------------------------------------------------

    /**
     * element-wise comparison that tolerates a length mismatch instead of
     * aborting on it. Armadillo's operator-= asserts equal sizes, and one file
     * may well hold a 2D table next to a 3D one
     */
    bool
    grids_match( const Vector< real > & aLeft,
                 const Vector< real > & aRight )
    {
        if ( aLeft.length() != aRight.length() )
        {
            return false ;
        }

        for ( index_t k=0; k<aLeft.length(); ++k )
        {
            // relative, so that a coordinate of order 1e3 is not held to an
            // absolute 2e-15
            real tScale = std::max( 1.0, std::abs( aRight( k ) ) ) ;

            if ( std::abs( aLeft( k ) - aRight( k ) ) > BELFEM_EPSILON * tScale )
            {
                return false ;
            }
        }
        return true ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    bool
    grids_match( const Vector< index_t > & aLeft,
                 const Vector< index_t > & aRight )
    {
        if ( aLeft.length() != aRight.length() )
        {
            return false ;
        }

        for ( index_t k=0; k<aLeft.length(); ++k )
        {
            if ( aLeft( k ) != aRight( k ) )
            {
                return false ;
            }
        }
        return true ;
    }

//------------------------------------------------------------------------------
}

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // create Arguments ( this also reads the shared -v / --verbose flags )
    Arguments tArguments( argc, argv );

    if ( gComm.rank() == 0 )
    {
        string tPath = input_file_from_arguments( tArguments );

        if ( tPath.empty() || help_requested( tArguments ) )
        {
            print_help();
            return gComm.finalize();
        }

        print_banner( "db2exo" );

        string tName = exodus_path( tPath );

        // the exodus writer opens its target with EX_CLOBBER
        // ( cl_Mesh_ExodusWriter.cpp ), so an input that already carries the
        // exodus extension would be destroyed by its own conversion. The test
        // is on the extension rather than on tName != tPath, because the
        // latter passes for file.EXO -- a different string that is the same
        // file on the case-insensitive volumes this project also ships on
        BELFEM_ERROR( string_to_lower( filetype( tPath ) ) != "exo",
            "Input %s already carries the exodus extension and would be overwritten by its own output. Rename it first.",
            tPath.c_str() );

        HDF5 tFile( tPath, FileMode::OPEN_RDONLY );

        Cell< string > tGroups = hdf5::get_groups( tFile.active_group() );

        Mesh * tMesh = nullptr ;

        uint              tOrder = 0 ;
        Vector< index_t > tPoints ;
        Vector< real >    tOrigin ;
        Vector< real >    tStep ;

        // number of tables that became a field on the mesh
        uint tNumTables = 0 ;

        for ( const string & tGroup : tGroups )
        {
            tFile.select_group( tGroup ) ;
            hid_t tLoc = tFile.active_group() ;

            if ( ! is_table( tLoc ) )
            {
                tFile.close_active_group() ;
                continue ;
            }

            tFile.load_data( "order", tOrder );
            tFile.load_data( "origin", tOrigin );
            tFile.load_data( "step", tStep );
            tFile.load_data( "points", tPoints );

            // a table exodus cannot hold is skipped, not fatal -- otherwise a
            // single odd table sitting first in name order would abort the run
            // and hide every convertible table behind it
            if ( tPoints.length() != 2 && tPoints.length() != 3 )
            {
                message( InfoLevel::Default,
                    "    skipping table %s : it is %u-dimensional, exodus takes 2D and 3D only",
                    tGroup.c_str(),
                    ( unsigned int ) tPoints.length() );

                tFile.close_active_group() ;
                continue ;
            }

            if ( tMesh == nullptr )
            {
                // the first convertible table defines the grid all others share
                BELFEM_ERROR( tOrigin.length() == tPoints.length()
                           && tStep.length() == tPoints.length(),
                    "Table %s of %s is inconsistent : points, origin and step are of length %u, %u and %u",
                    tGroup.c_str(),
                    tPath.c_str(),
                    ( unsigned int ) tPoints.length(),
                    ( unsigned int ) tOrigin.length(),
                    ( unsigned int ) tStep.length() );

                // TensorMeshConfig evaluates ( points - 1 ) % order while
                // initializing mNumElementsPerDim, which is declared ahead of
                // mElementType and therefore runs BEFORE the order is
                // validated : order 0 divides by zero and a point count of 0
                // underflows the unsigned subtraction. The file is external
                // input, so it is checked here, on the always-active tier
                BELFEM_ERROR( tOrder >= 1 && tOrder <= 3,
                    "Table %s of %s has interpolation order %u, only 1, 2 and 3 are supported",
                    tGroup.c_str(),
                    tPath.c_str(),
                    ( unsigned int ) tOrder );

                for ( index_t i=0; i<tPoints.length(); ++i )
                {
                    BELFEM_ERROR( tPoints( i ) > 1,
                        "Table %s of %s has %lu grid points in direction %u, at least 2 are needed",
                        tGroup.c_str(),
                        tPath.c_str(),
                        ( long unsigned int ) tPoints( i ),
                        ( unsigned int ) i );

                    // a zero step turns into an infinite inverse element step
                    // in TensorMeshConfig and a degenerate mesh here
                    BELFEM_ERROR( std::isfinite( tStep( i ) ) && tStep( i ) > 0.0,
                        "Table %s of %s has a step of %g in direction %u, it must be finite and positive",
                        tGroup.c_str(),
                        tPath.c_str(),
                        tStep( i ),
                        ( unsigned int ) i );

                    BELFEM_ERROR( std::isfinite( tOrigin( i ) ),
                        "Table %s of %s has a non-finite origin in direction %u",
                        tGroup.c_str(),
                        tPath.c_str(),
                        ( unsigned int ) i );
                }

                tMesh = new Mesh( tOrder, tPoints, tStep, tOrigin );

                // exodus carries the corner connectivity only for QUAD16 and
                // HEX64 ( ExodusWriter::fix_num_nodes ), so a cubic table
                // keeps all of its nodes and their values but connects 8 of
                // every 64. Accepted for now -- said out loud rather than left
                // for the reader to discover in ParaView
                if ( tOrder > 2 )
                {
                    message( InfoLevel::Default,
                        "    warning : order %u tables are written with corner connectivity only,\n"
                        "              the interior sample points will not be connected to an element",
                        ( unsigned int ) tOrder );
                }
            }
            else if ( tOrder != tMesh->tensorconf()->order()
                   || ! grids_match( tPoints, tMesh->tensorconf()->num_nodes_vector() )
                   || ! grids_match( tOrigin, tMesh->tensorconf()->origin() )
                   || ! grids_match( tStep,   tMesh->tensorconf()->step() ) )
            {
                // one exodus file carries one mesh, so a table on a different
                // grid cannot travel with the others. Say so rather than drop
                // it in silence
                message( InfoLevel::Default,
                    "    skipping table %s : its grid differs from the one of the first table",
                    tGroup.c_str() );

                tFile.close_active_group() ;
                continue ;
            }

            Vector< real > & tValues = tMesh->create_field( tGroup );
            tFile.load_data( "values", tValues );

            // load_data resizes the field, and the exodus writer only asserts
            // the length, so a release build would write a short buffer
            BELFEM_ERROR( tValues.length() == tMesh->number_of_nodes(),
                "Table %s of %s holds %lu values, but its grid has %lu nodes",
                tGroup.c_str(),
                tPath.c_str(),
                ( long unsigned int ) tValues.length(),
                ( long unsigned int ) tMesh->number_of_nodes() );

            ++tNumTables ;

            tFile.close_active_group() ;
        }
        tFile.close();

        BELFEM_ERROR( tMesh != nullptr,
            "No tensor table found in %s. A table is a group holding the datasets order, origin, step, points and values.",
            tPath.c_str() );

        tMesh->save( tName ) ;

        message( InfoLevel::Default,
            "    wrote %u table%s into %s",
            ( unsigned int ) tNumTables,
            tNumTables == 1 ? "" : "s",
            tName.c_str() );

        delete tMesh ;
    }

    return gComm.finalize();
}
