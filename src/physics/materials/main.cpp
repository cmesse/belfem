//
// Created by Christian Messe on 27.07.20.
//

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Material.hpp"
#include "cl_Arguments.hpp"
#include "en_Materials.hpp"
#include "cl_MaterialFactory.hpp"
#include "fn_linspace.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    // create Arguments
    Arguments tArguments( argc, argv );

    string tLabel = tArguments.data( 1 ) ;

    // get material enum
    MaterialType tMaterialType =  string_to_material_type( tArguments.data( 1 ) ) ;

    // material factory
    MaterialFactory tFactory ;

    // create the material
    Material * tMatierial = tFactory.create_material( tMaterialType );

    if( tArguments.data().size() > 2 )
    {
        real tRRR = std::stod( tArguments.data( 2 ) ) ;
        tMatierial->set_rrr( tRRR );
    }

    if( gComm.rank() == gComm.size() -1 )
    {
        print_banner();

        if( tMatierial->number().length() > 0 )
        {
            std::cout << "Material: " << tMatierial->label() << " " << tMatierial->number() << std::endl;
        }
        else
        {
            std::cout << "Material: " << tMatierial->label() << std::endl;
        }


       /*Vector< real > tT = {
                4.0,
                7.5,
                10.0,
                15.0,
                20.0,
                25.0,
                30.0,
                40.0,
                50.0,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
                125.0,
                150.0,
                200.0,
                225.0,
                250.0,
                273.15,
                298.15,
                350.0,
                400.0,
                500.0,
                600.0,
                700.0,
                800.0,
                900.0,
                1000.0,
                1200.0,
                1400.0,
                1600.0,
                1800.0,
                2000.0,
                2200.0,
                2400.0,
                2400.0,
                2600.0,
                2800.0,
                3000.0 }; */


        Vector< real > tT = linspace( 4.0, 300.0, 297 );
              
        // scan material list
        uint n = tT.length();
        for ( uint k = 0; k < tT.length(); ++k )
        {
            if ( tT( k ) > tMatierial->T_max())
            {
                tT( k ) = tMatierial->T_max();
                n = k + 1;
                break;
            }
        }

        // compute case
        uint tCase = tMatierial->has_thermal_expansion() ? 1 : 0;
        tCase += tMatierial->has_mechanical() ? 2 : 0;
        tCase += tMatierial->has_electric_resistivity() ? 4 : 0;
        tCase += tMatierial->has_thermal() ? 8 : 0;

        // magnetic field
        real tB = 0.0 ;

        // check if a b-field was given
        if( tArguments.data().size() > 3 )
        {
            tB = std::stod( tArguments.data( 3 ) ) ;
        }

        switch ( tCase )
        {
            case ( 1 ) :
            {
                // alpha
                printf( "      T    alpha\n" );
                printf( "      K    1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f\n",
                            tT( k ),
                            tMatierial->alpha( tT( k )) * 1.e6 );
                }
                break;
            }
            case ( 2 ) :
            {
                printf( "      T    E        nu\n" );
                printf( "      K    GPa      -\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %8.4f\n",
                            tT( k ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k ))
                    );
                }
                break;
            }
            case ( 3 ) :
            {
                //  E&nu, alpha
                printf( "      T    E        nu    alpha\n" );
                printf( "      K    GPa      -     1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %8.4f %10.3f\n",
                            tT( k ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )),
                            tMatierial->alpha( tT( k )) * 1.e6 );
                }
                break;
            }
            case ( 4 ) :
            {
                // rho_el
                printf( "      T    rho_e0\n" );
                printf( "      K    10-9 A/m²\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f\n",
                            tT( k ),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ));
                }
                break;
            }
            case ( 5 ) :
            {
                // rho_el, alpha
                printf( "      T    rho_e0       alpha\n" );
                printf( "      K    1e-9 A/m²     1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f\n",
                            tT( k ),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ),
                            1e6 * tMatierial->alpha( tT( k )));
                }
                break;
            }
            case ( 6 ) :
            {
                // rho_el, E&nu
                // cp&lambda, rho_el, E^nu
                printf( "      T          rho_e0       E        nu\n" );
                printf( "      K       10-9 A/m²     GPa        - \n" );
                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %8.4f %10.3f\n",
                            tT( k ),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )));
                }
                break;
            }
            case ( 7 ) :
            {
                // rho_el, E&nu, alpha
                printf( "      T          rho_e0       E        nu       alpha\n" );
                printf( "      K       10-9 A/m²     GPa        -         1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %8.4f %10.3f %10.3f\n",
                            tT( k ),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )),
                            1e6 * tMatierial->alpha( tT( k )));
                }
                break;
            }
            case ( 8 ) :
            {
                // cp&lambda
                printf( "      T          c    lambda\n" );
                printf( "      K    J/(kg*K)   W/(m*K)\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k ))
                    );
                }
                break;
            }
            case ( 9 ) :
            {
                // cp&lambda, alpha
                printf( "      T          c    lambda      alpha\n" );
                printf( "      K    J/(kg*K)   W/(m*K)     1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            tMatierial->alpha( tT( k )) * 1.e6
                    );
                }
                break;
            }
            case ( 10 ) :
            {
                // cp&lambda, E&nu,
                printf( "      T          c    lambda      E        nu\n" );
                printf( "      K    J/(kg*K)   W/(m*K)     GPa      -\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %10.3f %10.3f %8.4f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k ))
                    );
                }
                break;
            }
            case ( 11 ) :
            {
                // cp&lambda, E&nu, alpha
                printf( "      T          c    lambda      E        nu        alpha\n" );
                printf( "      K    J/(kg*K)   W/(m*K)     GPa      -         1e-6/K\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            1.e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )),
                            tMatierial->alpha( tT( k )) * 1.e6
                    );
                }
                break;
            }
            case ( 12 ) :
            {
                // cp&lambda, rho_el
                printf( "      T          c    lambda    rho_b0\n" );
                printf( "      K    J/(kg*K)   W/(m*K)    1e-9 A/m²\n" );
                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            tMatierial->rho_el( 0.0, tT( k ), tB ) * 1e9
                    );
                }
                break;
            }
            case ( 13 ) :
            {
                // cp&lambda, rho_el, alpha
                break;
            }
            case ( 14 ) :
            {
                // cp&lambda, rho_el, E^nu
                printf( "      T          c    lambda       rho_e0       E        nu\n" );
                printf( "      K    J/(kg*K)   W/(m*K)   1e-9 A/m²     GPa        - \n" );
                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )));
                }
                break;
            }
            case ( 15 ) :
            {
                // cp&lambda, rho_el, E^nu, alpha
                printf( "      T          c    lambda      rho_e0       E        nu        alpha\n" );
                printf( "      K    J/(kg*K)   W/(m*K)  1e-9 A/m²     GPa        -         1e-6/K\n" );
                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), tB ),
                            1e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )),
                            1e6 * tMatierial->alpha( tT( k )));
                }
                break;
            }
        }
    }
    printf( "density @ %3.2f K : %10.3f kg/m²\n", BELFEM_TREF, tMatierial->rho( BELFEM_TREF ) );

    delete tMatierial ;

    return gComm.finalize();
}
