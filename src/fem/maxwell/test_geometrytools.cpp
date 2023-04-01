//
// Created by christian on 3/1/23.
//
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_HDF5.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Block.hpp"
#include "FEM_geometry.hpp"
#include "fn_norm.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "cl_Pipette.hpp"
#include "cl_EdgeFunctionFactory.hpp"

using namespace belfem ;
using namespace fem ;
using namespace geometry ;

Communicator gComm;
Logger       gLog( 5 );

bool
EXPECT_NEAR( const real aA, const real aB, const real aEpsilon )
{
    const bool aCheck = std::abs( aA - aB ) < aEpsilon ;
    if( aCheck )
    {
        std::cout << "    pass! " << std::endl ;
    }
    else
    {
        std::cout << "    FAIL! | " << aA << " - " << aB << " | = " << std::abs( aA - aB )  <<  " > " << aEpsilon << std::endl ;
        BELFEM_ERROR( false, "fail");
    }
    return aCheck ;
}

bool
EXPECT_TRUE( const bool aBool )
{
    if( aBool )
    {
        std::cout << "    pass! " << std::endl ;
    }
    else
    {
        std::cout << "    FAIL!" << std::endl ;
        BELFEM_ERROR( false, "fail");
    }
    return aBool ;
}

//------------------------------------------------------------------------------


void test_edge_function( HDF5 * aDatabase  , ElementType aType )
{
    Mesh * tMesh = geometry::create_test_mesh( aDatabase, "edgefunctions", aType );
    tMesh->create_edges( true );
    tMesh->finalize_edges( tMesh->elements() );


    // create a new group and one element
    Group * tGroup = new Group( nullptr, GroupType::BLOCK, aType, 1, 1, true );
    fem::Element * tElement = new fem::Element( tMesh->elements()( 0 ) );
    tGroup->elements().set_size( 1, tElement );

    uint tNumDimensions = mesh::dimension( aType );

    // open the tree in the database
    aDatabase->select_group( "edgefunctions" );
    aDatabase->select_group( to_string( aType ) );


    // get the integration points
    Matrix< real > tPoints ;
    aDatabase->load_data( "Points", tPoints );
    uint tNumPoints = tPoints.n_cols() ;

    // results from edge interpolation
    Matrix< real > tEdge ;
    aDatabase->load_data( "Edge", tEdge );

    // results from curl interpolation
    Matrix< real > tCurl ;
    aDatabase->load_data( "Curl", tCurl );
    aDatabase->close_tree() ;

    fem::EdgeFunctionFactory tFactory ;
    fem::EdgeFunction * tFunction = tFactory.create_edge_function( aType );

    // matrix with node coordinates
    Matrix< real > & tX = tGroup->work_X() ;
    tElement->get_node_coors( tX );

    // number of points

    // perform precomputation
    tFunction->precompute( tPoints );

    tFunction->link( tElement, true );

    for( uint e=0 ; e<tElement->element()->number_of_edges(); ++e )
    {
        std::cout << "edge " << e << " : " << tElement->edge_direction( e ) << std::endl ;

    }
    uint tCount = 0 ;
    uint tNumDofs = tFunction->E( 0 ).n_cols() ;
    Matrix< real > tExpect( tNumDimensions, tNumDofs );
    Matrix< real > tCompute( tNumDimensions, tNumDofs );

    Vector< real > tError( tNumDimensions );



    if( tNumDimensions == 2 )
    {
        // reset counter
        tCount = 0 ;

        // test edge function
        for( uint k=0; k<tNumPoints; ++k )
        {
            tExpect.set_row( 0, tEdge.row( tCount++ ) );
            tExpect.set_row( 1, tEdge.row(tCount++ ) );

            tCompute = tFunction->E( k ) ;
            std::cout << "======= " << k << " ======" << std::endl ;
            
            tExpect.print("Expect");
            tCompute.print("Compute");

            tError = tCompute - tExpect ;

            tError( 0 ) = norm( tCompute.row( 0 ) - tExpect.row( 0 ) ) ;
            tError( 1 ) = norm( tCompute.row( 1 ) - tExpect.row( 1 ) ) ;

            tError.print("Err");
            EXPECT_NEAR( norm( tError ), 0, BELFEM_EPSILON );
        }

        tExpect.set_size( 1, tNumDofs );
        tCompute.set_size( 1, tNumDofs );

        // test curl function
        tCount = 0 ;
        for( uint k=0; k<tNumPoints; ++k )
        {
            tExpect.set_row( 0, tCurl.row( tCount++ ) );
            tCompute = tFunction->C( k );

            tExpect.print("E");
            tCompute.print("C");
            real tR = norm( tExpect.row( 0 ) - tCompute.row( 0 ) );
            EXPECT_NEAR( tR, 0, BELFEM_EPSILON );
        }
    }

    delete tFunction ;

    // delete the group ( will also delete the element )
    delete tGroup ;

    // delete the mesh
    delete tMesh ;
}

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    // todo: only test on one core
    // todo: make path relative
    const string tPath ="/home/christian/codes/belfem/test/fem/test_database.hdf5" ;

    // load the file
    HDF5 tFile( tPath, FileMode::OPEN_RDONLY );

    //test_edge_function( & tFile, ElementType::TRI3 );
    test_edge_function( & tFile, ElementType::TRI6 );

    //EXPECT_TRUE( test_interface( & tFile, ElementType::TRI3 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TRI6 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::QUAD4 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::QUAD9 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TET4 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TET10 ) );
    //EXPECT_TRUE(  test_interface( & tFile, ElementType::HEX8 ) );
    //EXPECT_TRUE(  test_interface( & tFile, ElementType::HEX20 ) );
    //EXPECT_TRUE(  test_interface( & tFile, ElementType::HEX27 ) );

    tFile.close() ;


    // close communicator
    return gComm.finalize();
}
