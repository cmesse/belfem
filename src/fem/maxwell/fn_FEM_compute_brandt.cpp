//
// Created by Christian Messe on 22.02.22.
//

#include "fn_FEM_compute_brandt.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_MaxwellMaterial.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        real
        compute_brandt( DofManager * aField, const real aHaHp )
        {
            // applied field, this is the return value
            real aB = BELFEM_QUIET_NAN ;

            // identify width and height
            if ( aField->parent()->is_master() )
            {


                // grab the first block
                Block * tBlock = aField->block( 1 );

                BELFEM_ASSERT( tBlock->domain_type() == DomainType::Conductor,
                              "Block 1 must be of type superconductor" );

                // get the material from that block
                const MaxwellMaterial * tMat = reinterpret_cast< const MaxwellMaterial * > (
                        tBlock->material());


                // get the field
                Vector< real > & tField = aField->mesh()->field_exists( "brandt" ) ?
                                          aField->mesh()->field_data( "brandt" )
                                                                                   : aField->mesh()->create_field(
                                "brandt" );

                // reset the data
                tField.fill( BELFEM_QUIET_NAN );


                // the center should be zero, but we compute it anyways
                real tX0 = 0.25 * ( aField->mesh()->node( 1 )->x()
                                    + aField->mesh()->node( 2 )->x()
                                    + aField->mesh()->node( 3 )->x()
                                    + aField->mesh()->node( 4 )->x());

                real tY0 = 0.25 * ( aField->mesh()->node( 1 )->y()
                                    + aField->mesh()->node( 2 )->y()
                                    + aField->mesh()->node( 3 )->y()
                                    + aField->mesh()->node( 4 )->y());

                // unflag all nodes on the mesh
                aField->mesh()->unflag_all_nodes();

                // flag all nodes on the block
                tBlock->block()->flag_nodes();

                // get the node container
                Cell< mesh::Node * > & tNodes = aField->mesh()->nodes();

                // get the critical current density
                real Jc = tMat->jc();

                // half width of superconductor
                real a = std::abs( 0.5 * ( aField->mesh()->node( 2 )->x()
                                           - aField->mesh()->node( 1 )->x()));

                // half height of superconductor
                real b = std::abs( 0.5 * ( aField->mesh()->node( 3 )->y()
                                           - aField->mesh()->node( 2 )->y()));

                // compute the penetration field Eq. ( 65 )
                real Hp = Jc * b / constant::pi * (
                        2. * a / b * std::atan( b / a )
                        + std::log( 1. + a * a / ( b * b )));

                // compute the applied field
                real Ha = aHaHp * Hp;

                std::cout << "brandt " << a << " " << b << " " << Jc << " " << Ha << " " << Hp << std::endl ;
                std::cout << "org" << tX0 << " " << tY0 << std::endl ;

                // Eq. ( 70, part 1 )
                real x0 = a/std::cosh( constant::pi * Ha / ( 2.0 * Jc * b ) );

                // help factors
                real phi = 4.0 * Jc * b / constant::pi ;
                real psi = std::sqrt( a*a-x0*x0 )/a ;

                // first sector
                for ( mesh::Node * tNode: tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        // grab node coordinates
                        real x = std::abs( tNode->x() - tX0 );

                        // compute sign
                        real s = tNode->x() < tX0 ? -1.0 : 1.0 ;

                        if( x > x0 - BELFEM_EPSILON )
                        {
                            tField( tNode->index() ) = s*Jc ;
                        }
                        else
                        {
                            // Eq. ( 70 )
                            real Js = phi * std::atan( psi * x / std::sqrt( x0*x0-x*x ) );

                            // Eq. ( 71 )
                            real yc = b*( 1.-0.5*Js/(Jc * b) );

                            // get y-coordinate of node
                            real y = std::abs( tNode->y() - tY0 );
                            if( y < yc )
                            {
                                tField( tNode->index() ) = 0.0 ;
                            }
                            else
                            {
                                tField( tNode->index() ) = s*Jc ;
                            }
                        }
                    }
                }

                // compute applied field
                aB = constant::mu0 * Ha ;
            }

            // wait for other procs
            comm_barrier() ;

            // sync values
            broadcast( aField->parent()->master(), aB );

            return aB ;
        }

//------------------------------------------------------------------------------
    }
}