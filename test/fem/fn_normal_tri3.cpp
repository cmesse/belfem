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

TEST( GEOMETRY, normal_tri3 )
{

    // create a new group
    Group * tGroup = new SideSet( ElementType::LINE2,
                                  ElementType::TRI3,
                                  ElementType::TRI3 );

    // matrix for Jacobian is also 2d!
    Matrix< real > & tJ = tGroup->work_J() ;

    // work vector for X-Coordinates
    Matrix< real > & tX = tGroup->work_Xm() ;

    // work vector for normals
    Matrix< real  > tNormals( 3, 2 );

    // work vector for expected value
    Vector< real > tExpect( 2 );

    // work vector for edge lengths
    Vector< real > tEdges( 3 );

    // surface
    real tSurface ;

    // load data from HDF5 file
    gDatabase->select_group("normals" );
    gDatabase->select_group("tri3" );
    gDatabase->load_data("Nodes", tX );
    gDatabase->load_data( "Normals", tNormals );
    gDatabase->load_data( "Surface", tSurface );
    gDatabase->load_data( "Edges", tEdges );
    gDatabase->close_tree();

    // test computation of normals
    for( uint e=0; e<3; ++e )
    {
        // compute normal
        const Vector< real > & tCompute = normal_tri( tGroup, e );

        // check value from datafile
        tExpect = tNormals.col( e );

        // test normal
        EXPECT_NEAR( norm( tCompute - tExpect ), 0, 1E-11 );

        // test edge length
        EXPECT_NEAR( 2.0 * tGroup->dx(), tEdges( e ), 1E-11 );

        // here we know that J is constant
        tJ = tGroup->master_integration( 0 )->dNdXi( 0 ) * tX ;

        // perform pipette test
        EXPECT_NEAR( test_pipette( tGroup ) , tSurface , 1E-11 );
    }

    // delete the group
    delete tGroup ;

}