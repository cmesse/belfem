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

#include "fn_read_node_tmatrix_from_hdf5.hpp"

namespace belfem
{
    namespace mesh
    {
        void
        read_node_tmatrix_from_hdf5( Mesh * aMesh, const string & aFile )
        {
            // open the file
            HDF5 tFile( aFile, FileMode::OPEN_RDONLY );

            tFile.select_group("HangingNodes" );

            Vector< id_t > tNodeIDs ;
            tFile.load_data( "NodeIDs", tNodeIDs );

            Matrix< id_t > tBasisIDs ;
            tFile.load_data( "BasisIDs", tBasisIDs );

            Matrix< real > tCoeffients ;
            tFile.load_data( "Coefficients", tCoeffients );

            uint tNumMatrices = tNodeIDs.length() ;


            uint tNumSources = tBasisIDs.n_cols() ;

            Cell< Node * > tSources( tNumSources, nullptr );
            Vector< real > tCoeffs( tNumSources );

            for( uint k=0; k<tNumMatrices; ++k )
            {
                // get target node
                Node * tTarget = aMesh->node( tNodeIDs( k ) );

                // populate containers for sources and coefficients
                for( uint i=0; i<tNumSources; ++i )
                {
                    tSources( i ) = aMesh->node( tBasisIDs( k, i ) );
                    tCoeffs( i )  = tCoeffients( k, i );
                }

                // create the matrix
                tTarget->set_sources( tSources, tCoeffs );
            }

            // close the file
            tFile.close() ;
        }
    }
}