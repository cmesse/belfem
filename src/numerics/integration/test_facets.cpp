/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"


#include "cl_Mesh.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_trans.hpp"
#include "fn_intpoints.hpp"
#include "en_IntegrationScheme.hpp"

// todo: the integration points for the facets are corrently in the interpolation library
//       it would be cleaner to move them to integration
#include "../../fem/interpolation/fn_IF_initialize_integration_points_on_facet.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );


void
tri3( const real xi, const real eta, Vector < real > & N )
{
    N =  { xi, eta, 1.0 - xi - eta } ;
}

void
quad4( const real xi, const real eta, Vector < real > & N )
{
    N = { 0.25*(1.-xi)*(1.-eta), 0.25*(1.+xi)*(1.-eta), 0.25*(1.+xi)*(1.+eta),0.25*(1.-xi)*(1.+eta) } ;
}


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    mesh::ElementFactory tFactory;

    auto tRef = tFactory.create_reference_element( ElementType::TET4 );

    mesh::Element * tElement = tRef->element() ;

    Matrix< uint > tTable ;
    tFactory.create_orientation_table( tElement->type(), tTable );

    uint tRow = 0 ;

    for ( uint f = 0 ; f<tElement->number_of_facets(); ++f )
    {
        // get the nodes
        Cell< mesh::Node * > tNodes ;
        tElement->get_nodes_of_facet(  f, tNodes );

        uint nn = tNodes.size();
        Matrix< real > Xm( nn, 3 );
        Matrix< real > Xs( nn, 3 );

        Vector< real > Nm( nn );
        Vector< real > Ns( nn );

        // determine type of face
        ElementType tType = mesh::element_type_from_numnodes( 2, nn);

        uint nnc = mesh::number_of_corner_nodes( tType );

        // get the interpolation function
        void ( *interp )( const real , const real , Vector < real > & ) = nnc == 3 ? tri3 : quad4 ;

        // get the node coordinates for the slave facet (always on element)
        for ( uint i = 0 ; i<nn; ++i )
        {
            Xs( i, 0 ) = tNodes( i )->x();
            Xs( i, 1 ) = tNodes( i )->y();
            Xs( i, 2 ) = tNodes( i )->z();

        }

        Matrix< real > xi_m ;
        Matrix< real > xi_s ;
        Vector< real > p( 3 ) ;
        Vector< real > q( 3 ) ;
        Vector< real > r( 3 ) ;

        uint tOrder = 5 ;
        Vector< real > w ;

        intpoints(
            IntegrationScheme::GAUSS,
            mesh::geometry_type( tType ),
            tOrder, w, xi_m );

        // loop over all orientations
        for ( uint o=0; o<nnc; ++o )
        {
            // get the coordinates for the master element
            for ( uint i = 0 ; i<nn; ++i )
            {
                uint j = tTable( i, tRow );
                Xm( i, 0 ) = tRef->nodes()(j)->x();
                Xm( i, 1 ) = tRef->nodes()(j)->y();
                Xm( i, 2 ) = tRef->nodes()(j)->z();
            }

            switch ( mesh::geometry_type( tElement->type() ) )
            {
                case GeometryType::TET :
                {
                    fem::facetintpoints::intpoints_tet( f, o, w, xi_s, tOrder );
                    break;
                }
                case GeometryType::HEX :
                {
                    fem::facetintpoints::intpoints_hex( f, o, w, xi_s, tOrder );
                    break ;
                }
                case GeometryType::PENTA :
                {
                    fem::facetintpoints::intpoints_penta( f, o, w, xi_s, tOrder );
                    break ;
                }

                // note: pyra still missing
                default:
                {
                    BELFEM_ERROR( false, "Invalid element type" );
                }
            }

            uint ng = w.length() ;
            for ( uint k=0; k<ng; ++k )
            {
                // master shape function
                ( * interp )( xi_m( 0, k ), xi_m( 1, k ), Nm ) ;

                // slave shape function
                ( * interp )( xi_s( 0, k ), xi_s( 1, k ), Ns ) ;

                p( 0 ) = dot( Nm, Xm.col( 0 ) );
                p( 1 ) = dot( Nm, Xm.col( 1 ) );
                p( 2 ) = dot( Nm, Xm.col( 2 ) );

                q( 0 ) = dot( Ns, Xs.col( 0 ) );
                q( 1 ) = dot( Ns, Xs.col( 1 ) );
                q( 2 ) = dot( Ns, Xs.col( 2 ) );

                r = p - q ;

                std::cout << "#test f : " << f << " o: " << o << " eps: " << norm( r ) << std::endl ;
            }
            ++tRow ;
        }
    }


    delete tRef ;

    return  gComm.finalize();
}