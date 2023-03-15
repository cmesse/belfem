//
// Created by christian on 3/14/23.
//

//
// Created by christian on 3/14/23.
//

#include <gtest/gtest.h>
#include "typedefs.hpp"

#include "cl_HDF5.hpp"
#include "cl_FEM_SideSet.hpp"
#include "FEM_geometry.hpp"
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

TEST( GEOMETRY, normal_tri6 )
{
    // create a new group
    Group * tGroup = new SideSet( ElementType::LINE3,
                                  ElementType::TRI6,
                                  ElementType::TRI6 );

    // work vector for X-Coordinates
    Matrix< real > & tX = tGroup->work_Xm() ;

    // work vector for normals
    Matrix< real  > tNormals ;

    // work vector for expected value
    Vector< real > tExpect( 2 );

    // work vector for edge lengths
    Vector< real > tEdges( 3 );

    // surface
    real tSurface ;

    // load data from HDF5 file
    gDatabase->select_group("normals" );
    gDatabase->select_group("tri6");
    gDatabase->load_data("Nodes", tX );
    gDatabase->load_data( "Normals", tNormals );
    gDatabase->load_data( "Surface", tSurface );
    gDatabase->load_data( "Edges", tEdges );


    // integration points
    Matrix< real > tPoints ;
    gDatabase->load_data( "Points", tPoints );

    // integration weights
    Vector< real > tWeights ;
    gDatabase->load_data( "Weights", tWeights  );

    // close file
    gDatabase->close_tree();

    // make sure that integration points are in correct order
    EXPECT_TRUE( test_intpoints( tGroup, tPoints, tWeights ) );

    uint tNumIntpoints = tWeights.length() ;

    index_t tCount = 0 ;

    // loop over all edges and comput normals and edge lengths
    for( uint e=0; e<3; ++e )
    {
        // edge length
        real tLength = 0 ;

        // loop over all integration points
        for( uint k = 0; k<tNumIntpoints; ++k )
        {
            const Vector< real > & tN = normal_tri( tGroup, e, k );

            real tError = norm( tN - tNormals.col( tCount++ ) ) ;

            EXPECT_NEAR( tError, 0 , 1E-11 );
            tLength += tWeights( k ) * tGroup->dx() ;
        }
        EXPECT_NEAR( tLength , tEdges( e ) , 1E-11 );
    }

    // compute volume
    EXPECT_NEAR( test_pipette( tGroup ) , tSurface , 1E-11 );

    // delete the group
    delete tGroup ;
}