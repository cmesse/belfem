//
// Created by Christian Messe on 5/28/26.
//
#include "banner.hpp"
#include "cl_Communicator.hpp"
#include "cl_CsvFile.hpp"
#include "cl_HDF5.hpp"
#include "cl_Logger.hpp"
#include "cl_Matrix.hpp"
#include "cl_Timer.hpp"
#include "cl_Vector.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "constants.hpp"
#include "cl_DatabaseProjector.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

//------------------------------------------------------------------------------

/**
 * B-spline projection driver for the HTS lookup tables.
 *
 * usage: mattest <extended-table.hdf5>
 *
 * The file is written by the theta-extension script: it carries the groups
 * jc and n with the sampled nodal values, plus a root vector "grid" holding
 * [ nT, nB, nA, dT, dB, dA, oT, oB, oA ]. Reading the grid from the file
 * keeps this driver free of hardcoded table dimensions -- the tables differ
 * in the number of field nodes.
 *
 * The projector imposes no boundary condition. The gradient condition in
 * theta comes from the extension itself, and which one depends on the table:
 * three TRANSLATED copies of one period for a pi-periodic table, so the
 * solution one period away from the open ends matches in value AND slope at
 * theta = 0 and pi; EVEN reflection at both planes for a table folded about
 * pi/2, which gives a zero gradient there instead. The driver does not care
 * which -- it projects whatever extension it is handed.
 *
 * The projected nodal values are written back as "values_bspline" beside
 * the input; the fold script picks them up.
 */
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    BELFEM_ERROR( argc > 1, "usage: mattest <extended-table.hdf5>" );

    HDF5 tFile( argv[ 1 ], FileMode::OPEN_RDWR );

    Vector< real > tGrid ;
    tFile.load_data( "grid", tGrid );

    BELFEM_ERROR( tGrid.length() == 9,
        "the grid vector of %s must have 9 entries, but has %u",
        argv[ 1 ], ( unsigned int ) tGrid.length() );

    Mesh tMesh ( 2,
        { ( index_t ) tGrid( 0 ), ( index_t ) tGrid( 1 ), ( index_t ) tGrid( 2 ) },
        {             tGrid( 3 ),             tGrid( 4 ),             tGrid( 5 ) },
        {             tGrid( 6 ),             tGrid( 7 ),             tGrid( 8 ) } );

    database::Projector tProjector( &tMesh );

    Vector< real > & jc = tMesh.create_field( "jc" );
    Vector< real > & n  = tMesh.create_field( "n" );

    tFile.select_group( "jc" );
    tFile.load_data( "values", jc );
    tFile.close_active_group();

    tFile.select_group( "n" );
    tFile.load_data( "values", n );
    tFile.close_active_group();

    Vector< real > & jc2 = tMesh.create_field( "jc2" );
    Vector< real > & n2  = tMesh.create_field( "n2" );

    tProjector.project( "jc", jc2 );
    tProjector.project( "n", n2 );

    tFile.select_group( "jc" );
    tFile.save_data( "values_bspline", jc2 );
    tFile.close_active_group();

    tFile.select_group( "n" );
    tFile.save_data( "values_bspline", n2 );
    tFile.close_active_group();

    tFile.close();

    std::cout << "projected " << tMesh.number_of_nodes() << " nodes" << std::endl;

    return gComm.finalize();
}
