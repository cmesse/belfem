//
// Created by christian on 3/1/23.
//
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_HDF5.hpp"
#include "cl_FEM_SideSet.hpp"
#include "FEM_geometry.hpp"
#include "fn_norm.hpp"
#include "fn_det.hpp"
#include "fn_dot.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "cl_Pipette.hpp"
#include "cl_Element_Factory.hpp"

using namespace belfem ;
using namespace fem ;

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
#
bool
test_interface( HDF5 * aDatabase , ElementType aType )
{
    uint tNumPermutations = 0 ;
    switch (  mesh::geometry_type( aType )  )
    {
        case( GeometryType::TRI ):
        {
            tNumPermutations = 3 ;
            break;
        }
        case( GeometryType::QUAD ):
        {
            tNumPermutations = 4 ;
            break;
        }
        case( GeometryType::TET ):
        {
            tNumPermutations = 12 ;
            break;
        }
        case( GeometryType::HEX ):
        {
            tNumPermutations = 24 ;
            break;
        }
        default:
        {
            BELFEM_ERROR( false, "unsupported geometry type");
        }
    }


    uint tNumFacets = mesh::number_of_facets( aType );

    for( uint p=0; p<tNumPermutations; ++p )
    {
        // load mesh
        Mesh * tMesh = geometry::test_create_mesh( aDatabase,
                                                   aType,
                                                   p );

        tMesh->save( "/tmp/mesh.exo");

        for( uint f=0; f<tNumFacets; ++f )
        {
            mesh::Facet * tFacet = tMesh->facets()( f );

            IntegrationData * tSurface = new IntegrationData( tFacet->element()->type() );
            tSurface->populate( 0 );

            IntegrationData * tMaster = new IntegrationData( tFacet->master()->type() );
            tMaster->populate_for_master( tFacet->master_index(), 0 );

            IntegrationData * tSlave = new IntegrationData( tFacet->slave()->type() );

            switch (  mesh::geometry_type( aType )  )
            {
                case( GeometryType::TRI ):
                case( GeometryType::QUAD ) :
                {
                    tSlave->populate_for_slave_tri( tFacet->slave_index(), 0 );
                    break;
                }
                case( GeometryType::TET ) :
                {
                    tSlave->populate_for_slave_tet( tFacet->slave_index(),
                                                    tFacet->orientation_on_slave(),
                                                    0 );
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    tSlave->populate_for_slave_hex( tFacet->slave_index(),
                                                    tFacet->orientation_on_slave(),
                                                    0 );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "unsupported geometry type");
                }
            }

            // surface
            uint tNumNodesFacet = tFacet->number_of_nodes() ;
            Matrix< real > tX( tNumNodesFacet, 3 );

            // collect surface nodes
            for( uint k=0; k<tFacet->number_of_nodes(); ++k )
            {
                tX( k, 0 ) = tFacet->node( k )->x() ;
                tX( k, 1 ) = tFacet->node( k )->y() ;
                tX( k, 2 ) = tFacet->node( k )->z() ;
            }


            // master coords
            uint tNumNodesMaster = tFacet->master()->number_of_nodes() ;
            Matrix< real > tXm( tNumNodesMaster, 3 );
            for( uint k=0; k<tFacet->master()->number_of_nodes(); ++k )
            {
                tXm( k, 0 ) = tFacet->master()->node( k )->x() ;
                tXm( k, 1 ) = tFacet->master()->node( k )->y() ;
                tXm( k, 2 ) = tFacet->master()->node( k )->z() ;
            }

            // slave coords
            uint tNumNodesSlave = tFacet->slave()->number_of_nodes() ;
            Matrix< real > tXs( tNumNodesSlave, 3 );
            for( uint k=0; k<tFacet->slave()->number_of_nodes(); ++k )
            {
                tXs( k, 0 ) = tFacet->slave()->node( k )->x() ;
                tXs( k, 1 ) = tFacet->slave()->node( k )->y() ;
                tXs( k, 2 ) = tFacet->slave()->node( k )->z() ;
            }

            Vector< real > tP(3) ;
            Vector< real > tQ(3) ;
            Vector< real > tR(3) ;

            for( uint k=0; k<tSurface->number_of_integration_points(); ++k )
            {
                tP( 0 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 0 )  );
                tP( 1 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 1 )  );
                tP( 2 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 2 )  );
            }

            for( uint k=0; k<tSurface->number_of_integration_points(); ++k )
            {
                tQ( 0 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 0 )  );
                tQ( 1 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 1 )  );
                tQ( 2 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 2 )  );
                if( norm( tP - tQ ) > 1e-12 )
                {
                    return false ;
                }
            }

            for( uint k=0; k<tSurface->number_of_integration_points(); ++k )
            {
                tR( 0 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 0 )  );
                tR( 1 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 1 )  );
                tR( 2 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 2 )  );

                if ( norm( tP - tR ) > 1e-12 )
                {
                    return false ;
                }
            }




            delete tSurface ;
            delete tMaster ;
            delete tSlave ;
        }

        delete tMesh;
    }
    return true ;
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
    //const string tPath ="/home/christian/codes/belfem/test/fem/test_geometry.hdf5" ;
    const string tPath ="/Users/christian/codes/belfem_tests/geometry/test_geometry.hdf5" ;
    // load the file
    HDF5 tFile( tPath, FileMode::OPEN_RDONLY );

    //EXPECT_TRUE( test_interface( & tFile, ElementType::TRI3 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TRI6 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::QUAD4 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::QUAD9 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TET4 ) );
    //EXPECT_TRUE( test_interface( & tFile, ElementType::TET10 ) );
    test_interface( & tFile, ElementType::HEX8 );

    tFile.close() ;


    // close communicator
    return gComm.finalize();
}
