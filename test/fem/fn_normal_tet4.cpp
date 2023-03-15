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

TEST( GEOMETRY, normal_tet4 )
{
    // create a new group
    Group * tGroup = new SideSet( ElementType::TRI3,
                                  ElementType::TET4,
                                  ElementType::TET4 );

    // vector with node coordinates
    Matrix< real > & tX = tGroup->work_Xm() ;
    Matrix< real > & tJ = tGroup->work_J() ;


    // matrix with normals
    Matrix< real > tNormals ;

    // vector with surfaces
    Vector< real > tSurfaces( 4 );

    // expected value
    Vector< real > tExpect( 3 );

    real tVolume = 0.0 ;

    // load data from HDF5 file
    gDatabase->select_group("normals");
    gDatabase->select_group("tet4" );
    gDatabase->load_data("Nodes", tX );
    gDatabase->load_data( "Normals", tNormals );
    gDatabase->load_data( "Surfaces", tSurfaces );
    gDatabase->load_data( "Volume", tVolume );
    gDatabase->close_tree() ;

    // Jacobian is constant for linear element
    tJ = tGroup->master_integration( 0 )->dNdXi( 0 ) * tX ;

    for( uint f=0; f<4; ++f )
    {
        const Vector< real > & tCompute = normal_tet( tGroup, f );

        // check value from datafile
        tExpect = tNormals.col( f );

        // test normal
        EXPECT_NEAR( norm( tCompute - tExpect ), 0, 1E-11 );

        // test surface

        EXPECT_NEAR( tGroup->dx(), tSurfaces( f ) , 1E-11 );
    }

    // compute volume
    EXPECT_NEAR( test_pipette( tGroup ), tVolume, 1E-11 );

    // delete the group
    delete tGroup ;
}