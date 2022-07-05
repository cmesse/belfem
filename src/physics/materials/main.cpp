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

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    // create Arguments
    Arguments tArguments( argc, argv );

    string tLabel = tArguments.data( 1 ) ;

    // get material enum
    MaterialType tMaterialType =  string_to_material_type( tArguments.data( 1 ) ) ;

    // material factory
    MaterialFactory tFactory ;

    // create the material
    Material * tMatierial = tFactory.create_material( tMaterialType );

    print_banner();

    std::cout << "Material: " << tMatierial->label() << std::endl ;

    Vector< real > tT = {
            4.0,
            10.0,
            20.0,
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
            3000.0 };

    // scan material list
    uint n=tT.length() ;

    for( uint k=0; k<tT.length(); ++k )
    {
        if( tT( k ) > tMatierial->T_max() )
        {
            tT( k ) = tMatierial->T_max();
            n = k + 1 ;
            break ;
        }
    }

    if( tMatierial->has_thermal() )
    {
        if( ! tMatierial->has_mechanical() )
        {
            if( ! tMatierial->has_electric_resistivity() )
            {
                printf( "      T          c    lambda\n");
                printf( "      K    J/(kg*K)   W/(m*K)\n");


                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k ))
                    );
                }
            }
            else
            {
                printf( "      T          c    lambda     ho_b0\n");
                printf( "      K    J/(kg*K)   W/(m*K)    10-9 A/m²\n");
                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {
                    printf( " %6.0f %10.3f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            tMatierial->rho_el( 0.0, tT( k ), 0.0 )*1e9
                    );
                }
            }
        }
        else if( ! tMatierial->has_thermal_expansion() )
        {
            if( ! tMatierial->has_electric_resistivity() )
            {
                if( ! tMatierial->has_electric_resistivity() )
                {
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
                }
                else
                {
                    printf( "      T          c    lambda      E        nu     rho_e0\n" );
                    printf( "      K    J/(kg*K)   W/(m*K)     GPa      -      10-9 A/m²\n" );
                                                               printf( "T          c" );
                    // loop over all entries
                    for ( uint k = 0; k < n; ++k )
                    {

                        printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f\n",
                                tT( k ),
                                tMatierial->c( tT( k )),
                                tMatierial->lambda( tT( k )),
                                1e-9 * tMatierial->E( tT( k )),
                                tMatierial->nu( tT( k )),
                                1e9 * tMatierial->rho_el( 0.0, tT( k ), 0.0 ) );
                    }
                }
            }
        }
        else
        {
            if( ! tMatierial->has_electric_resistivity() )
            {
                printf( "      T          c    lambda      E        nu        alpha\n" );
                printf( "      K    J/(kg*K)   W/(m*K)     GPa      -         10-6\n" );

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
            }
            else
            {
                printf( "      T          c    lambda      E        nu        alpha     rho_e0\n" );
                printf( "      K    J/(kg*K)   W/(m*K)     GPa      -         10-6      10-9 A/m²\n" );

                // loop over all entries
                for ( uint k = 0; k < n; ++k )
                {

                    printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f %10.3f\n",
                            tT( k ),
                            tMatierial->c( tT( k )),
                            tMatierial->lambda( tT( k )),
                            1.e-9 * tMatierial->E( tT( k )),
                            tMatierial->nu( tT( k )),
                            tMatierial->alpha( tT( k )) * 1.e6,
                            1e9 * tMatierial->rho_el( 0.0, tT( k ), 0.0 ) );
                }
            }
        }
    }


    delete tMatierial ;

    return gComm.finalize();
}