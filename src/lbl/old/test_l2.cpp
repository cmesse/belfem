//
// Created by Christian Messe on 03.02.22.
//

//
// Created by Christian Messe on 10.01.22.
//

#include <iostream>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Mesh.hpp"
#include "cl_Logger.hpp"
#include "cl_Element_Factory.hpp"
#include "en_IWGs.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell_L2_old.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

using namespace belfem;
using namespace fem ;

Communicator gComm;
Logger       gLog( 4 );

void init_2d( Mesh & aMesh )
{
    Vector< real > & tBx = aMesh.create_field("bx");
    Vector< real > & tBy = aMesh.create_field( "by");
    Vector< real > & tIz = aMesh.create_field( "iz");

    for( mesh::Node * tNode : aMesh.nodes() )
    {
        real tX = tNode->x() ;
        real tY = tNode->y() ;

        real tR = std::sqrt( tX * tX + tY * tY ) ;

        if ( abs( tR ) < 1e-12 )
        {
            tBx( tNode->index() ) = 0.0;
            tBy( tNode->index() ) = 0.0;
        }
        else
        {
            tBx( tNode->index() ) = -tY / tR;
            tBy( tNode->index() ) = tX / tR;
            tIz( tNode->index() ) = 1 / tR;
        }
    }

    Vector< real > & tBx0 = aMesh.create_field( "b0x");
    Vector< real > & tBy0 = aMesh.create_field( "b0y");

    tBx0 = tBx ;
    tBy0 = tBy ;
}

void init_3d( Mesh & aMesh )
{
    Vector< real > & tBx = aMesh.create_field("bx");
    Vector< real > & tBy = aMesh.create_field( "by");
    Vector< real > & tBz = aMesh.create_field( "bz");
    tBz.fill( 0.0 );

    Vector< real > & tIx = aMesh.create_field("ix");
    Vector< real > & tIy = aMesh.create_field( "iy");
    Vector< real > & tIz = aMesh.create_field( "iz");
    tIx.fill( 0.0 );
    tIy.fill( 0.0 );
    tIz.fill( 0.0 );

    for( mesh::Node * tNode : aMesh.nodes() )
    {
        real tX = tNode->x() + 2.0 ;
        real tY = tNode->y() + 1.0 ;
        real tR = std::sqrt( tX * tX + tY * tY )  ;

        if ( abs( tR ) < 1e-12 )
        {
            tBx( tNode->index() ) = 0.0;
            tBy( tNode->index() ) = 0.0;
        }
        else
        {
            tBx( tNode->index() ) = -tY / tR;
            tBy( tNode->index() ) = tX / tR;
            tIz( tNode->index() ) = 1 / tR ;
        }
    }

    Vector< real > & tBx0 = aMesh.create_field( "gx");
    Vector< real > & tBy0 = aMesh.create_field( "gy");
    Vector< real > & tBz0 = aMesh.create_field( "gz");
    tBx0 = tBx ;
    tBy0 = tBy ;
    tBz0 = tBz ;
}

void
check_mesh( Mesh * aMesh )
{
    Cell< mesh::Element * > & tElements = aMesh->elements();

    Vector< real > & tCheck = aMesh->create_field( "check", EntityType::ELEMENT );
    tCheck.fill( 1.0 );


    Vector< real > tA( 3 );
    Vector< real > tB( 3 );
    Vector< real > tC( 3 );
    Vector< real > tD( 3 );
    Vector< real > tM( 3 );
    Vector< real > tN( 3 );
    Vector< real > tP( 3 );
    Vector< real > tQ( 3 );

    for( mesh::Element * tElement : tElements )
    {
        if ( mesh::geometry_type( tElement->type() ) == GeometryType::TET )
        {
            // loop over all surfaces
            for( uint f=0; f<4; ++f )
            {
                uint i ;
                uint j ;
                uint k ;
                uint l ;

                // get nodes of surface
                switch( f )
                {
                    case( 0 ) :
                    {
                        i = 0 ;
                        j = 1 ;
                        k = 3 ;
                        l = 2 ;
                        break ;
                    }
                    case( 1 ) :
                    {
                        i = 1 ;
                        j = 2 ;
                        k = 3 ;
                        l = 0 ;
                        break ;
                    }
                    case( 2 ) :
                    {
                        i = 2 ;
                        j = 0 ;
                        k = 3 ;
                        l = 1 ;
                        break ;
                    }
                    default :
                    {
                        i = 0 ;
                        j = 2 ;
                        k = 1 ;
                        l = 3 ;
                        break ;
                    }
                }

                tA(0) = tElement->node( i )->x() ;
                tA(1) = tElement->node( i )->y() ;
                tA(2) = tElement->node( i )->z() ;

                tB(0) = tElement->node( j )->x() ;
                tB(1) = tElement->node( j )->y() ;
                tB(2) = tElement->node( j )->z() ;

                tC(0) = tElement->node( k )->x() ;
                tC(1) = tElement->node( k )->y() ;
                tC(2) = tElement->node( k )->z() ;

                tD(0) = tElement->node( l )->x() ;
                tD(1) = tElement->node( l )->y() ;
                tD(2) = tElement->node( l )->z() ;

                tM = 1/3 * ( tA + tB + tC );
                tA -= tM ;
                tB -= tM ;
                tC -= tM ;
                tD -= tM ;

                tP = tB - tA ;
                tQ = tC - tA ;

                tN = cross( tP, tQ );
                tN /= norm( tN );

                if( dot( tM, tN ) > 0.0 )
                {
                    tCheck( tElement->index() ) = 0.0 ;
                }
            }
        }
    }
}

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

     // Mesh tMesh( "mesh_tet10.msh" );
    //Mesh tMesh( "quadrangle.msh" );
    //Mesh tMesh( "twotets3.msh" );
    Mesh tMesh( "box.msh" );
    tMesh.create_edges();
    tMesh.create_faces();
    tMesh.flag_curved_elements();

    check_mesh( &tMesh );

    init_3d( tMesh );


    KernelParameters tParams( tMesh );

    Kernel tKernel( &tParams );

    ElementType tType = tMesh.block( 1 )->element_type() ;

    IWG_Maxwell_L2 * tL2a = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::B2H );
    IWG_Maxwell_L2 * tL2b = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2B );
    IWG_Maxwell_L2 * tL2c = new IWG_Maxwell_L2( tType, Maxwell_L2_Mode::H2J );


    tL2a->set_block( 1, DomainType::Conductor );
    tL2b->set_block( 1, DomainType::Conductor );
    tL2c->set_block( 1, DomainType::Conductor );

    DofManager * tFa = tKernel.create_field( tL2a );
    DofManager * tFb = tKernel.create_field( tL2b );
    DofManager * tFc = tKernel.create_field( tL2c );

    // create materials
    MaxwellMaterial * tMat = new MaxwellMaterial( "superconductor" );

    // set an arbitrary prime number so that it is not zero
    tMat->set_rho_el_const( 1.0 );
    tMat->set_mu_r( 1.0 );  // todo: enforce 1.0 for this iwg

    // set the material type for the superconductor ( it will be deleted by the Dof Manager)
    tFa->block( 1 )->set_material( tMat );
    tFb->block( 1 )->set_material( tMat );
    tFc->block( 1 )->set_material( tMat );
#ifdef BELFEM_MUMPS
    tFa->set_solver( SolverType::MUMPS );
    tFb->set_solver( SolverType::MUMPS );
    tFc->set_solver( SolverType::MUMPS );
#else
    tFa->set_solver( SolverType::UMFPACK );
    tFb->set_solver( SolverType::UMFPACK );
    tFc->set_solver( SolverType::UMFPACK );
#endif

    tFa->block( 1 )->set_domain_type( DomainType::Conductor );
    tFb->block( 1 )->set_domain_type( DomainType::Conductor );
    tFc->block( 1 )->set_domain_type( DomainType::Conductor );

    tL2a->set_field( tFa );
    tL2b->set_field( tFb );
    tL2c->set_field( tFc );

    tFa->compute_jacobian(  );
    tFa->compute_rhs();
    tFa->solve() ;
    tFa->save_system("matrix.hdf5");

    tMesh.field_data("bx").fill( 0.0 );
    tMesh.field_data("by").fill( 0.0 );

    tFb->compute_jacobian();
    tFb->compute_rhs();
    tFb->solve();

    /*for ( uint f=0; f<tMesh.number_of_fields(); ++f )
    {
        std::cout << f << " " << tMesh.field( f )->label() << " " << tMesh.field(f )->data().length() << std::endl ;
    }*/


    // tMesh.save("mesh.hdf5");
    tMesh.create_field("ElJx", EntityType::ELEMENT );
    tMesh.create_field("ElJy", EntityType::ELEMENT );
    tMesh.create_field("ElJz", EntityType::ELEMENT );

    tFc->compute_jacobian();
    tFc->compute_rhs();
    tFc->solve();

    /*for ( uint f=0; f<tMesh.number_of_fields(); ++f )
    {
        std::cout << f << " " << tMesh.field( f )->label() << " " << tMesh.field(f )->data().length() << std::endl ;
    } */


    tMesh.save( "mesh.exo");


    // close communicator
    return gComm.finalize();
}