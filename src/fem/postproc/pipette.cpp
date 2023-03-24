//
// Created by christian on 8/23/22.
//
#include "cl_Communicator.hpp"
#include "cl_Pipette.hpp"
#include "cl_Element_Factory.hpp"
#include "random.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_rotation_matrix.hpp"
#include "constants.hpp"

using namespace belfem;
using namespace mesh;

Communicator gComm;

//------------------------------------------------------------------------------

real
triangle_surface( const Node * aA, const Node * aB, const Node * aC )
{
    return 0.5*std::sqrt(
            std::pow(-(-aA->x()+aB->x())*(-aA->y()+aC->y())+(-aA->x()+aC->x())*(-aA->y()+aB->y()),2)
             +std::pow((-aA->x()+aB->x())*(-aA->z()+aC->z())-(-aA->x()+aC->x())*(-aA->z()+aB->z()),2)
             +std::pow(-(-aA->y()+aB->y())*(-aA->z()+aC->z())+(-aA->y()+aC->y())*(-aA->z()+aB->z()),2));

}

//------------------------------------------------------------------------------

real
distance( const Node * aA, const Node * aB, const Node * aC, const Node * aD )
{
    Vector< real > tA( 3 );
    tA( 0 ) = aA->x() ;
    tA( 1 ) = aA->y() ;
    tA( 2 ) = aA->z() ;


    Vector< real > tB( 3 );
    tB( 0 ) = aB->x() ;
    tB( 1 ) = aB->y() ;
    tB( 2 ) = aB->z() ;

    Vector< real > tC( 3 );
    tC( 0 ) = aC->x() ;
    tC( 1 ) = aC->y() ;
    tC( 2 ) = aC->z() ;

    Vector< real > tD( 3 );
    tD( 0 ) = aD->x() ;
    tD( 1 ) = aD->y() ;
    tD( 2 ) = aD->z() ;

    // compute the normal of the plane
    Vector< real > tN( cross( tB-tA,tC-tA) );
    tN /= norm( tN );

    // center of plane
    Vector< real > tP( tA + tB + tC );
    tP /= 3.0 ;

    // distance of plane from origin
    real tH = dot( tC, tN );

    // the distance of the point p to the plane is
    return std::abs( dot( tN, tD ) - tH );
}

//------------------------------------------------------------------------------

real
tet_volume( const Node * aA, const Node * aB, const Node * aC, const Node * aD )
{
    Vector< real > tA( 3 );
    tA( 0 ) = aA->x() ;
    tA( 1 ) = aA->y() ;
    tA( 2 ) = aA->z() ;


    Vector< real > tB( 3 );
    tB( 0 ) = aB->x() ;
    tB( 1 ) = aB->y() ;
    tB( 2 ) = aB->z() ;

    Vector< real > tC( 3 );
    tC( 0 ) = aC->x() ;
    tC( 1 ) = aC->y() ;
    tC( 2 ) = aC->z() ;

    Vector< real > tD( 3 );
    tD( 0 ) = aD->x() ;
    tD( 1 ) = aD->y() ;
    tD( 2 ) = aD->z() ;

    Vector< real > tAB( tB - tA );
    Vector< real > tAC( tC - tA );
    Vector< real > tAD( tD - tA );

    return std::abs( dot( tAB, cross( tAC, tAD ) ) )/6.0 ;
}

//------------------------------------------------------------------------------
bool
test_tri3()
{
    random_seed() ;

    ElementType tType = ElementType::TRI3 ;

    // create the nodes
    mesh::Node * tA = new mesh::Node( 1, -belfem::rand(), -belfem::rand(), 0 );
    mesh::Node * tB = new mesh::Node( 2, belfem::rand(), -belfem::rand(), 0 );
    mesh::Node * tC = new mesh::Node( 3, belfem::rand(), belfem::rand(), 0 );


    mesh::ElementFactory tFactory ;

    Element * tElement = tFactory.create_element( tType, 1 );
    tElement->insert_node( tA, 0 );
    tElement->insert_node( tB, 1 );
    tElement->insert_node( tC, 2 );


    // reference solution
    real tRef = triangle_surface( tA, tB, tC );

    Pipette tPipette;

    tPipette.set_element_type( tType );

    real tValue = tPipette.measure( tElement );

    delete tElement ;
    delete tA ;
    delete tB ;
    delete tC ;

    return std::abs( tValue - tRef ) < BELFEM_MESH_EPSILON ;
}

//------------------------------------------------------------------------------

bool
test_quad4()
{
    random_seed() ;

    ElementType tType = ElementType::QUAD4 ;

    // create the nodes
    real tEPS = 1e-6 ;

    mesh::Node * tA = new mesh::Node( 1, -belfem::rand()-tEPS, -belfem::rand()-tEPS, 0 );
    mesh::Node * tB = new mesh::Node( 2, belfem::rand()+tEPS, -belfem::rand()-tEPS, 0 );
    mesh::Node * tC = new mesh::Node( 3, belfem::rand()+tEPS, belfem::rand()+tEPS, 0 );
    mesh::Node * tD = new mesh::Node( 4, -belfem::rand()-tEPS, belfem::rand()+tEPS, 0 );

    mesh::ElementFactory tFactory ;

    Element * tElement = tFactory.create_element( tType, 1 );
    tElement->insert_node( tA, 0 );
    tElement->insert_node( tB, 1 );
    tElement->insert_node( tC, 2 );
    tElement->insert_node( tD, 3 );

    // reference solution
    real tRef =   triangle_surface( tA, tB, tC )
                + triangle_surface( tA, tC, tD );

    Pipette tPipette;

    tPipette.set_element_type( tType );

    real tValue = tPipette.measure( tElement );

    delete tElement ;
    delete tA ;
    delete tB ;
    delete tC ;
    delete tD ;

    return std::abs( tValue - tRef ) < BELFEM_MESH_EPSILON ;
}

bool
test_tet4()
{
    random_seed() ;

    ElementType tType = ElementType::TET4 ;

    // create the nodes
    mesh::Node * tA = new mesh::Node( 1, belfem::rand(), belfem::rand(), -belfem::rand() );
    mesh::Node * tB = new mesh::Node( 2, belfem::rand(), belfem::rand(), -belfem::rand() );
    mesh::Node * tC = new mesh::Node( 3, belfem::rand(), belfem::rand(), -belfem::rand() );
    mesh::Node * tD = new mesh::Node( 4, belfem::rand(), belfem::rand(), belfem::rand() );

    mesh::ElementFactory tFactory ;

    Element * tElement = tFactory.create_element( tType, 1 );
    tElement->insert_node( tA, 0 );
    tElement->insert_node( tB, 1 );
    tElement->insert_node( tC, 2 );
    tElement->insert_node( tD, 3 );


    // reference solution
    real tRef = tet_volume( tA, tB, tC, tD );

    Pipette tPipette;

    tPipette.set_element_type( tType );

    real tValue = tPipette.measure( tElement );


    delete tElement ;
    delete tA ;
    delete tB ;
    delete tC ;
    delete tD ;

    return std::abs( tValue - tRef ) < BELFEM_MESH_EPSILON ;
}

bool
test_penta6()
{
    random_seed() ;

    ElementType tType = ElementType::PENTA6 ;

    // baseline of arbitrary triangle
    real tG = belfem::rand() ;

    // height of arbitrary triangle
    real tH = belfem::rand() ;

    // coordinate of perpendicular point
    real tXi = belfem::rand() ;

    // length of penta
    real tL = belfem::rand() ;

    // create the nodes
    mesh::Node * tA = new mesh::Node( 1, 0, 0, 0 );
    mesh::Node * tB = new mesh::Node( 2, tG, 0, 0 );
    mesh::Node * tC = new mesh::Node( 3, tG*tXi, tH, 0 );
    mesh::Node * tD = new mesh::Node( 4, 0, 0, tL );
    mesh::Node * tE = new mesh::Node( 5, tG, 0, tL );
    mesh::Node * tF = new mesh::Node( 6, tG*tXi, tH, tL );

    // now we create an arbitraty rotation matrix
    Matrix< real > tR( 3, 3 );

    real tYaw   = (belfem::rand()-0.5) * 2.0 * constant::pi ;
    real tPitch = (belfem::rand()-0.5) * 2.0 * constant::pi ;
    real tRoll  = (belfem::rand()-0.5) * 2.0 * constant::pi ;
    rotation_matrix( tYaw, tPitch, tRoll, tR );

    Vector< real > tX( 3 );

    // rotate the node coordinates
    tX = tR * tA->coords() ;
    tA->set_coords( tX );

    tX = tR * tB->coords() ;
    tB->set_coords( tX );

    tX = tR * tC->coords() ;
    tC->set_coords( tX );

    tX = tR * tD->coords() ;
    tD->set_coords( tX );

    tX = tR * tE->coords() ;
    tE->set_coords( tX );

    tX = tR * tF->coords() ;
    tF->set_coords( tX );

    mesh::ElementFactory tFactory ;

    Element * tElement = tFactory.create_element( tType, 1 );
    tElement->insert_node( tA, 0 );
    tElement->insert_node( tB, 1 );
    tElement->insert_node( tC, 2 );
    tElement->insert_node( tD, 3 );
    tElement->insert_node( tE, 4 );
    tElement->insert_node( tF, 5 );

    // reference solution
    real tRef =   tet_volume( tA, tB, tC, tF )
                + tet_volume( tD, tE, tF, tA )
                + tet_volume( tA, tB, tE, tF );

    Pipette tPipette;

    tPipette.set_element_type( tType );

    real tValue = tPipette.measure( tElement );

    delete tElement ;
    delete tA ;
    delete tB ;
    delete tC ;
    delete tD ;
    delete tE ;
    delete tF ;

    return std::abs( tValue - tRef ) < BELFEM_MESH_EPSILON ;
}

bool
test_hex8()
{
    random_seed() ;

    ElementType tType = ElementType::HEX8 ;


    // create the nodes
    // create the nodes
    real tEPS = 1e-6 ;
    mesh::Node * tA = new mesh::Node( -belfem::rand(), -belfem::rand()-tEPS, -belfem::rand()-tEPS );
    mesh::Node * tB = new mesh::Node(  belfem::rand(), -belfem::rand()-tEPS, -belfem::rand()-tEPS );
    mesh::Node * tC = new mesh::Node(  belfem::rand(),  belfem::rand()+tEPS, -belfem::rand()-tEPS );
    mesh::Node * tD = new mesh::Node( -belfem::rand(),  belfem::rand()+tEPS, -belfem::rand()-tEPS );
    mesh::Node * tE = new mesh::Node( -belfem::rand(), -belfem::rand()-tEPS, belfem::rand()+tEPS );
    mesh::Node * tF = new mesh::Node(  belfem::rand(), -belfem::rand()-tEPS, belfem::rand()+tEPS );
    mesh::Node * tG = new mesh::Node(  belfem::rand(),  belfem::rand()+tEPS, belfem::rand()+tEPS );
    mesh::Node * tH = new mesh::Node( -belfem::rand(),  belfem::rand()+tEPS, belfem::rand()+tEPS );

    mesh::ElementFactory tFactory ;

    Element * tElement = tFactory.create_element( tType, 1 );
    tElement->insert_node( tA, 0 );
    tElement->insert_node( tB, 1 );
    tElement->insert_node( tC, 2 );
    tElement->insert_node( tD, 3 );
    tElement->insert_node( tE, 4 );
    tElement->insert_node( tF, 5 );
    tElement->insert_node( tG, 6 );
    tElement->insert_node( tH, 7 );

    // reference solution
    real tRef =     tet_volume( tA, tB, tC, tF )
                  + tet_volume( tA, tC, tD, tH )
                  + tet_volume( tA, tF, tE, tH )
                  + tet_volume( tF, tC, tG, tH )
                  + tet_volume( tD, tC, tH, tA );


    Pipette tPipette;

    tPipette.set_element_type( tType );

    real tValue = tPipette.measure( tElement );

    delete tElement ;
    delete tA ;
    delete tB ;
    delete tC ;
    delete tD ;
    delete tE ;
    delete tF ;
    delete tG ;
    delete tH ;

    return std::abs( tValue - tRef ) < BELFEM_MESH_EPSILON ;
}

//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );


    std::cout << "tri3  " << test_tri3() << std::endl ;
    std::cout << "quad4 " << test_quad4() << std::endl ;
    std::cout << "tet4 " << test_tet4() << std::endl ;
    std::cout << "penta6 " << test_penta6() << std::endl ;
    std::cout << "hex8 " << test_hex8() << std::endl ;
    return gComm.finalize();
}
