/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Material.hpp"
#include "cl_Arguments.hpp"
#include "cl_MaterialFactory.hpp"
#include "fn_linspace.hpp"
#include "fn_create_database_mesh.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 3 );

void
print_help()
{
    std::cout << "Usage: material <material_label> [OPTION]" << std::endl;
    std::cout << "  -a, --angle <field_angle>                 field angle in degree" << std::endl;
    std::cout << "  -b, --field <field_strength>              magnetic field strength in Tesla" << std::endl;
    std::cout << "  -c, --create                              create hdf5 lookup table" << std::endl;
    std::cout << "  -m, --mesh                                save exodus mesh for visualization" << std::endl;
    std::cout << "  -r, --rrr   <value>                       set RRR value" << std::endl;
    std::cout << "  -t, --temperatures <tmin> <tmax> <step>   temperature range in K" << std::endl;
    std::cout << "  -h, --help                                print this help message" << std::endl;
    std::cout << "  -l, --list                                list builtin materials" << std::endl ;
}

void
print_tables( const Material * aMaterial, Vector< real > & aT, real tB, real tBeta )
{
    if( aMaterial->number().length() > 0 )
    {
        std::cout << "Material: " << aMaterial->label() << " " << aMaterial->number();
    }
    else
    {
        std::cout << "Material: " << aMaterial->label() ;
    }

    if ( aMaterial->type() == MaterialType::PureMetal || aMaterial->type() == MaterialType::CompositeAlloy )
    {
        std::cout << " RRR: " << aMaterial->constant_property( MaterialProperty::RRR ) ;
        std::cout << sprint( " B : %4.2f T   angle(B∠j): %5.2f °", tB, tBeta / constant::deg ) ;
    }
    std::cout << std::endl ;

    uint n = aT.length();
    for ( uint k = 0; k < aT.length(); ++k )
    {
        real aTmax = aMaterial->constant_property( MaterialProperty::T_max );
        if ( aT( k ) > aTmax )
        {
            aT( k ) = aTmax;
            n = k + 1;
            break;
        }
    }

    // compute case
    uint aCase = aMaterial->have( MaterialProperty::alpha )? 1 : 0;
    aCase += aMaterial->have( MaterialProperty::E) ? 2 : 0;
    aCase += aMaterial->have( MaterialProperty::rho )? 4 : 0;
    aCase += aMaterial->have( MaterialProperty::cp ) ? 8 : 0;

    // Field-aware accessors: only invoke the 3-arg overloads when the
    // material actually declares a dependency on |B|. Materials such as
    // HastelloyC276 declare have(rho)/have(lambda) but never populate the
    // magnetoresistance database, so the 3-arg overload would assert.
    auto rho_at = [&]( const real T ) -> real
    {
        return aMaterial->depends( MaterialProperty::rho, MaterialDependency::normB )
             ? aMaterial->rho( T, tB, tBeta )
             : aMaterial->rho( T );
    };
    auto lambda_at = [&]( const real T ) -> real
    {
        return aMaterial->depends( MaterialProperty::lambda, MaterialDependency::normB )
             ? aMaterial->lambda( T, tB, tBeta )
             : aMaterial->lambda( T );
    };

    switch ( aCase )
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
                        aT( k ),
                        aMaterial->alpha( aT( k )) * 1.e6 );
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
                        aT( k ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k ))
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
                        aT( k ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )),
                        aMaterial->alpha( aT( k )) * 1.e6 );
            }
            break;
        }
        case ( 4 ) :
        {
            // rho
            printf( "      T         rho\n" );
            printf( "      K    1e-8 Ω·m\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {

                printf( " %6.0f %10.3f\n",
                        aT( k ),
                        1e8 *  rho_at( aT( k ) ) );
            }
            break;
        }
        case ( 5 ) :
        {
            // rho, alpha
            printf( "      T    rho_0        alpha\n" );
            printf( "      K    1e-8 Ω·m      1e-6/K\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %10.3f\n",
                        aT( k ),
                        1e8 *  rho_at( aT( k ) ),
                        1e6 * aMaterial->alpha( aT( k )));
            }
            break;
        }
        case ( 6 ) :
        {
            // rho, E&nu
            // cp&lambda, rho, E^nu
            printf( "      T            rho        E        nu\n" );
            printf( "      K       1e-8 Ω·m      GPa        - \n" );
            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {

                printf( " %6.0f %10.3f %8.4f %10.3f\n",
                        aT( k ),
                        1e8 *  rho_at( aT( k ) ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )));
            }
            break;
        }
        case ( 7 ) :
        {
            // rho, E&nu, alpha
            printf( "      T            rho        E        nu       alpha\n" );
            printf( "      K       1e-8 Ω·m      GPa        -         1e-6/K\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %8.4f %10.3f %10.3f\n",
                        aT( k ),
                        1e8 *  rho_at( aT( k ) ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )),
                        1e6 * aMaterial->alpha( aT( k )));
            }
            break;
        }
        case ( 8 ) :
        {
            // cp&lambda
            printf( "      T         cp    lambda\n" );
            printf( "      K    J/(kg*K)   W/(m*K)\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        aMaterial->lambda( aT( k ))
                );
            }
            break;
        }
        case ( 9 ) :
        {
            // cp&lambda, alpha
            printf( "      T         cp    lambda      alpha\n" );
            printf( "      K    J/(kg*K)   W/(m*K)     1e-6/K\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {

                printf( " %6.0f %10.3f %10.3f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        aMaterial->lambda( aT( k )),
                        aMaterial->alpha( aT( k )) * 1.e6
                );
            }
            break;
        }
        case ( 10 ) :
        {
            // cp&lambda, E&nu,
            printf( "      T         cp    lambda      E        nu\n" );
            printf( "      K    J/(kg*K)   W/(m*K)     GPa      -\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {

                printf( " %6.0f %10.3f %10.3f %10.3f %8.4f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        lambda_at( aT( k ) ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k ))
                );
            }
            break;
        }
        case ( 11 ) :
        {
            // cp&lambda, E&nu, alpha
            printf( "      T         cp    lambda      E        nu        alpha\n" );
            printf( "      K    J/(kg*K)   W/(m*K)     GPa      -         1e-6/K\n" );

            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {

                printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        lambda_at( aT( k ) ),
                        1.e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )),
                        aMaterial->alpha( aT( k )) * 1.e6
                );
            }
            break;
        }
        case ( 12 ) :
        {
            // cp&lambda, rho
            printf( "      T         cp    lambda    rho_b0\n" );
            printf( "      K    J/(kg*K)   W/(m*K)    1e-8 Ω·m\n" );
            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %10.3f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        lambda_at( aT( k ) ),
                         rho_at( aT( k ) ) * 1e8
                );
            }
            break;
        }
        case ( 13 ) :
        {
            // cp&lambda, rho, alpha
            break;
        }
        case ( 14 ) :
        {
            // cp&lambda, rho, E^nu
            printf( "      T         cp    lambda         rho        E        nu\n" );
            printf( "      K    J/(kg*K)   W/(m*K)   1e-8 Ω·m      GPa        - \n" );
            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %10.3f %10.3f %8.4f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        lambda_at( aT( k ) ),
                        1e8 *  rho_at( aT( k ) ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )));
            }
            break;
        }
        case ( 15 ) :
        {
            // cp&lambda, rho, E, nu, alpha
            printf( "      T         cp    lambda         rho       E        nu        alpha\n" );
            printf( "      K    J/(kg*K)   W/(m*K)  1e-8 Ω·m      GPa        -         1e-6/K\n" );
            // loop over all entries
            for ( uint k = 0; k < n; ++k )
            {
                printf( " %6.0f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
                        aT( k ),
                        aMaterial->cp( aT( k )),
                        lambda_at( aT( k ) ),
                        1e8 *  rho_at( aT( k ) ),
                        1e-9 * aMaterial->E( aT( k )),
                        aMaterial->nu( aT( k )),
                        1e6 * aMaterial->alpha( aT( k )));
            }
            break;
        }
    }

    printf( "density @ %3.2f K : %10.3f kg/m³\n", BELFEM_TREF, aMaterial->density( BELFEM_TREF ) ) ;

    if ( aMaterial->is_constant( MaterialProperty::grueneisen ) )
    {
        printf( "Grüneisen param. γ : %5.4f\n", aMaterial->constant_property( MaterialProperty::grueneisen ) );
    }
}

void
save_mesh( Material * aMaterial )
{
    Mesh * tMesh = create_database_mesh(1);

    std::string tLabel = aMaterial->label() + ".exo";

    Cell< mesh::Node * > & tNodes = tMesh->nodes();

    if ( aMaterial->have( MaterialProperty::rho ) )
    {
        Vector< real > & rho = tMesh->create_field( "rho" );

        const bool tRhoB = aMaterial->depends(
                MaterialProperty::rho, MaterialDependency::normB );

        index_t tCount = 0 ;
        for ( mesh::Node * tNode : tNodes )
        {
            real T = tNode->x();
            real B = std::pow( 10, tNode->y() );
            real beta = tNode->z();
            rho( tCount++ ) = ( tRhoB ? aMaterial->rho( T, B, beta )
                                      : aMaterial->rho( T ) ) * 1e8 ;
        }
    }
    if ( aMaterial->have( MaterialProperty::lambda ) )
    {
        Vector< real > & lambda = tMesh->create_field( "lambda" );

        const bool tLambdaB = aMaterial->depends(
                MaterialProperty::lambda, MaterialDependency::normB );

        index_t tCount = 0 ;
        for ( mesh::Node * tNode : tNodes )
        {
            real T = tNode->x();
            real B = std::pow( 10, tNode->y() );
            real beta = tNode->z();
            lambda( tCount++ ) = tLambdaB ? aMaterial->lambda( T, B, beta )
                                          : aMaterial->lambda( T );
        }
    }

    // translate node coordinates
    for ( mesh::Node * tNode : tNodes )
    {
        real T = tNode->x();
        real B = std::pow( 10, tNode->y() );
        real beta = tNode->z() / constant::deg;
        tNode->set_coords( T, B, beta );
    }

    tMesh->save( tLabel );

    delete tMesh ;
}

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );

    // create Arguments
    Arguments tArguments( argc, argv );

    if ( argc < 2 )
    {
        print_help();
        return gComm.finalize();
    }

    string tLabel = tArguments.data( 1 ) ;

    if ( tLabel == "-h" || tLabel == "--help" )
    {
        print_help();
        return gComm.finalize();
    }

    uint na = tArguments.data().size() ;
    uint a = 1 ;

    // magnetic field
    real tB = 0.0 ;

    // field angle
    real tBeta = 0.0 ;

    Vector< real > tT = { 0, 4,5,6,7,8,9, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 77, 80, 90, 100, 150, 200, 250, 273.15 , 293.15 };

    // material factory
    MaterialFactory tFactory ;

    // create the material


    real Tmin  = BELFEM_QUIET_NAN ;
    real Tstep = BELFEM_QUIET_NAN ;

    // we test if this is an alloy, if so, we use a more conservative
    // default value of RRR 10, otherwise RRR 100
    Cell< std::pair< string, real > > tComposition ;
    to_pair( tLabel, tComposition );
    real tRRR = tComposition.size() > 1 ? 10. : 100. ;

    real tTMaxUser = 400.0 ;

    bool tSaveMesh = false ;
    bool tCreateTable = false ;

    while ( a < na )
    {
        if( tArguments.data( a ) == "-h" || tArguments.data( a ) == "--help" )
        {
            print_help();
            return gComm.finalize();
        }
        if( tArguments.data( a ) == "-l" || tArguments.data( a ) == "--list" )
        {
            tFactory.print_material_list( std::cout );
            return gComm.finalize();
        }
        else if( tArguments.data( a ) == "-a" || tArguments.data( a ) == "-A" || tArguments.data( a ) == "--angle" )
        {
            tBeta = std::stod( tArguments.data( ++a ) ) * constant::deg ;
        }
        else if( tArguments.data( a ) == "-b" || tArguments.data( a ) == "-B" || tArguments.data( a ) == "--field" )
        {
            tB = std::stod( tArguments.data( ++a ) ) ;
        }
        else if( tArguments.data( a ) == "-r" ||  tArguments.data( a ) == "-R" ||  tArguments.data( a ) == "--RRR" || tArguments.data( a ) == "--rrr" )
        {
            tRRR = std::stod( tArguments.data( ++a ) ) ;
        }
        else if( tArguments.data( a ) == "-m" ||  tArguments.data( a ) == "--mesh" )
        {
            tSaveMesh = true ;
            tCreateTable = true ;
        }
        else if( tArguments.data( a ) == "-c" ||  tArguments.data( a ) == "--create" )
        {
            tCreateTable = true ;
        }
        else if( tArguments.data( a ) == "-t" || tArguments.data( a ) == "-T" || tArguments.data( a ) == "--temperatures" )
        {
            try
            {
                Tmin  = std::stod( tArguments.data( a + 1 ) ) ;
            }
            catch ( std::exception & e )
            {
                Tmin = BELFEM_QUIET_NAN ;
            }

            if ( ! std::isnan( Tmin ) )
            {
                try
                {
                    tTMaxUser  = std::stod( tArguments.data( a + 2 ) ) ;
                }
                catch ( std::exception & e )
                {
                    tTMaxUser = 0 ;
                }

                try
                {
                    if ( ! std::isnan( tTMaxUser ) )
                    {
                        Tstep = std::stod( tArguments.data( a + 3 ) ) ;
                    }
                }
                catch ( std::exception & e )
                {
                    Tstep = 4.0 ;
                }

                uint n = (uint) ( ( tTMaxUser - Tmin ) / Tstep ) + 1 ;
                tT.set_size( n );
                real T = Tmin ;
                for ( uint i = 0; i < n; ++i)
                {
                    tT( i ) = T ;
                    T += Tstep ;
                }
                tT( n-1 ) = tTMaxUser ;
                a += 3 ;
            }
        }
        ++a ;
    }

    Material * tMaterial = tFactory.create_material( tLabel, tRRR, tCreateTable );

    if( gComm.rank() == 0 )
    {
        if ( tSaveMesh )
        {
            save_mesh( tMaterial );
        }
        else
        {
            print_tables( tMaterial, tT, tB, tBeta );
        }
    }

    delete tMaterial ;

    return gComm.finalize();
}
