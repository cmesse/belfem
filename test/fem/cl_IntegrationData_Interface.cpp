//
// Created by christian on 3/30/23.
//

#include <gtest/gtest.h>
#include "typedefs.hpp"

#include "cl_HDF5.hpp"
#include "cl_FEM_SideSet.hpp"
#include "FEM_geometry.hpp"
#include "Mesh_Enums.hpp"
#include "fn_norm.hpp"
#include "fn_det.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "cl_Pipette.hpp"
#include "cl_Element_Factory.hpp"


using namespace belfem ;
using namespace fem ;
using namespace geometry ;

extern belfem::HDF5 * gDatabase;

bool
test_interface( HDF5 * aDatabase , ElementType aType )
{
    // number of facets per element
    uint tNumFacets = mesh::number_of_facets( aType );

    // number of permutations
    uint tNumPermutations = 0 ;
    for( uint f=0; f<tNumFacets; ++f )
    {
        tNumPermutations += mesh::number_of_orientations( aType, f );
    }

    uint tOrder = auto_integration_order( aType );
    tOrder = 7 ;

    for( uint p=0; p<tNumPermutations; ++p )
    {
        // load mesh
        Mesh * tMesh = geometry::create_test_mesh( aDatabase,
                                                   "interfaces",
                                                   aType,
                                                   p );

        for( uint f=0; f<tNumFacets; ++f )
        {
            mesh::Facet * tFacet = tMesh->facets()( f );

            IntegrationData * tSurface = new IntegrationData( tFacet->element()->type() );
            tSurface->populate( tOrder );

            IntegrationData * tMaster = new IntegrationData( tFacet->master()->type() );
            tMaster->populate_for_master( tFacet->master_index(), tOrder );

            IntegrationData * tSlave = new IntegrationData( tFacet->slave()->type() );
            tSlave->populate_for_slave( tFacet->slave_index(),
                                        tFacet->orientation_on_slave(),
                                        tOrder );

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

            Vector< real > tP(3 ) ;
            Vector< real > tQ(3 ) ;
            Vector< real > tR(3 ) ;

            for( uint k=0; k<tSurface->number_of_integration_points(); ++k )
            {
                // point on surface
                tP( 0 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 0 )  );
                tP( 1 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 1 )  );
                tP( 2 ) = dot(  tSurface->phi( k ).vector_data(), tX.col( 2 )  );

                // point on master
                tQ( 0 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 0 )  );
                tQ( 1 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 1 )  );
                tQ( 2 ) = dot(  tMaster->phi( k ).vector_data(), tXm.col( 2 )  );

                // point on slave
                tR( 0 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 0 )  );
                tR( 1 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 1 )  );
                tR( 2 ) = dot(  tSlave->phi( k ).vector_data(), tXs.col( 2 )  );

                if( norm( tP - tQ ) > 1e-12 )
                {
                    return false ;
                }

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

TEST( GEOMETRY, interface_tri3 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::TRI3 ) ) ;
}

TEST( GEOMETRY, interface_tri6 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::TRI6 ) ) ;
}

TEST( GEOMETRY, interface_quad4 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::QUAD4 ) ) ;
}

TEST( GEOMETRY, interface_quad9 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::QUAD9 ) ) ;
}

TEST( GEOMETRY, interface_tet4 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::TET4 ) ) ;
}

TEST( GEOMETRY, interface_tet10 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::TET10 ) ) ;
}

TEST( GEOMETRY, interface_hex8 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::HEX8 ) ) ;
}

TEST( GEOMETRY, interface_hex27 )
{
    EXPECT_TRUE( test_interface( gDatabase, ElementType::HEX27 ) ) ;
}