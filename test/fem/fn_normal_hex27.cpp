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

TEST( GEOMETRY, normal_hex27 )
{
    // create a new group
    Group * tGroup = new SideSet( ElementType::QUAD9,
                                  ElementType::HEX27,
                                  ElementType::HEX27 );
    
    // vector with node coordinates
    Matrix< real > & tX = tGroup->work_Xm();
    
    // vector with surfaces
    Vector< real > tSurfaces( 4 );
    
    // expected value
    Vector< real > tExpect( 3 );
    
    // matrix with normals
    Matrix< real > tNormals ;
    
    // volume of hex
    real tVolume ;
    
    // matrix with normals
    Matrix< real > tPoints ;
    
    // integration weights
    Vector< real > tWeights ;

    gDatabase->select_group( "normals" );
    gDatabase->select_group( "hex27" );
    gDatabase->load_data( "Nodes", tX );
    gDatabase->load_data( "Normals", tNormals );
    gDatabase->load_data( "Surfaces", tSurfaces );
    gDatabase->load_data( "Volume", tVolume );
    gDatabase->load_data( "Points", tPoints );
    gDatabase->load_data( "Weights", tWeights );
    gDatabase->close_tree();

    EXPECT_TRUE( test_intpoints( tGroup, tPoints, tWeights ) );
    
    // initialize counter
    index_t tCount = 0 ;
    
    // number of integration pointstPoints = {belfem::Matrix<double>}
    uint tNumPoints = tWeights.length() ;
    
    // loop over all surfaces
    for( uint f=0; f<6; ++f )
    {
        // reset surface
        real tS = 0.0 ;
        
        // loop over all integraition points
        for( uint k=0 ; k<tNumPoints; ++k )
        {
            // compute the normal
            const Vector< real > & tCompute = normal_hex( tGroup, f, k );
            
            // grab the expected value
            tExpect = tNormals.col( tCount++ );
            
            // test normal
            EXPECT_NEAR( norm( tCompute - tExpect ), 0, 1E-11 );
            
            // compute surface
            tS += tWeights( k ) * tGroup->dx() ;
        }
        
        EXPECT_NEAR( tS, tSurfaces( f ) , 1E-11 );
    
    }
    
    // compute volume
    EXPECT_NEAR( test_pipette( tGroup ), tVolume, 1E-11 );
    
    delete tGroup ;
}