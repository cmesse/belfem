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
#include "commtools.hpp"

#include "fn_sum.hpp"
#include "fn_embed_python_guide.hpp"
#include "fn_dot.hpp"
#include "cl_Tensor.hpp"
#include "fn_identity.hpp"

#include "cl_Material_Alloy.hpp"
#include "cl_MaterialFactory.hpp"

#include "fn_invert_symmetric.hpp"
#include "fn_hust.hpp"
#include "cl_Mesh_Distributor.hpp"
#include "cl_HDF5.hpp"
#include "stringtools.hpp"
#include "fn_create_database_mesh.hpp"
#include "fn_rho_database_is_current.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace material
    {
        Alloy::Alloy( const string & aName, const bool aBuildTables ) :
            SplineLookupTable( MaterialType::PureMetal )
        {
            mComputeTables = aBuildTables ;
            this->set_label( aName.length() == 0 ? "Alloy" : aName  );

            Cell< std::pair< string, real > > tPairs ;
            to_pair( aName, tPairs );
            this->set_components( tPairs );

        }

        Alloy::Alloy( const string & aName,
                       const Cell< std::pair< string, real > > & aComponents,
                       const real aRRR,
                       const bool aBuildTables ):
            SplineLookupTable( MaterialType::PureMetal )
        {
            mComputeTables = aBuildTables ;
            this->set_label( aName.length() == 0 ? "Alloy" : aName  );

            // check if components contain an rrr value
            bool tHaveRRR = false ;
            Cell< std::pair< string, real > > tComponents ;
            for ( const std::pair< string, real > & tPair : aComponents )
            {
                if ( string_to_lower( tPair.first ) == "rrr" )
                {
                    tHaveRRR = true ;
                }
                tComponents.push( tPair );
            }

            if ( ! tHaveRRR && ! std::isnan( aRRR ) )
            {
                tComponents.push( std::make_pair( "RRR", aRRR ) );
            }
            this->set_components( tComponents );
        }

        Alloy::~Alloy()
        {
            this->delete_components();

            if ( mRhoData != nullptr )
            {
                delete mRhoData;
            }
        }

        void
        Alloy::set_components( const Cell< std::pair< string, real > > & aComponents )
        {
            if ( aComponents.size() > 0 )
            {
                Cell< string > tComponents ;
                Cell< real > tMassFractions ;

                for ( const std::pair< string, real > & tPair : aComponents )
                {
                    if ( string_to_lower( tPair.first ) == "rrr" )
                    {
                        this->set_constant( MaterialProperty::RRR, tPair.second );
                    }
                    else
                    {
                        tComponents.push(  tPair.first );
                        tMassFractions.push( tPair.second );
                    }
                }


                Vector< real > tX( tMassFractions.size() );
                for ( size_t k=0; k<tMassFractions.size(); ++k )
                {
                    tX( k ) = tMassFractions( k ) * 0.01 ;
                }

                this->set_components( tComponents, tX );
            }
        }

        void
        Alloy::delete_components()
        {
            for ( Metal * tMetal : mComponents ) delete tMetal ;
            mComponents.clear();
        }


        void
        Alloy::set_components(
               const Cell< string > & aComponents,
               const Vector< real > & aMassFractions,
               const real RRR  )
        {
            BELFEM_ERROR( aMassFractions.length() == aComponents.size(),
                "Mass fractions must be the same size as the number of components" );

            this->delete_components();

            uint n = aComponents.size() ;
            mComponents.set_size( n, nullptr );

            Vector< real > & M = mMolarMasses ;

            Vector< real > & Y = mMassFractions ;
            Vector< real > & X = mMolarFractions ;
            Vector< real > & V = mVolumeFractions ;

            Y = aMassFractions ;

            // last entry is always balance
            Y( n - 1) = 0.0 ;
            Y( n - 1 ) = 1.0 - sum( Y );

            X.set_size( n );
            M.set_size( n );
            V.set_size( n );

            real tTmax = BELFEM_REAL_MAX ;

            for ( uint i = 0; i < n; ++i )
            {
                mComponents(i) = this->create_component( aComponents(i) );
                M( i )   = mComponents(i)->constant_property( MaterialProperty::M );
                V( i ) = Y( i ) /  mComponents(i)->density( gTroom );
                X( i ) = Y( i ) / M( i );

                if ( mComponents(i)->constant_property( MaterialProperty::T_max ) < tTmax )
                {
                    tTmax = mComponents(i)->constant_property( MaterialProperty::T_max );
                }
            }
            V /= sum( V );
            X /= sum( X );

            this->set_constant( MaterialProperty::T_max, tTmax );
            this->set_constant( MaterialProperty::M, dot( X, M ) );

            // check if the RRR value has already been set
            real tRRR = RRR ;

            if ( std::isnan( tRRR ) && this->is_constant( MaterialProperty::RRR ) )
            {
                tRRR = this->constant_property( MaterialProperty::RRR );
            }

            // compute the alloy's residual resistivity from the user-specified RRR
            BELFEM_ERROR( !std::isnan( tRRR ) && tRRR > 1.0,
                "Alloy '%s' requires a valid RRR > 1 (got %f)",
                this->label().c_str(), tRRR );

            real tRhoI = 0.0 ;
            for ( uint i = 0; i < n; ++i )
            {
                tRhoI += V( i ) * mComponents( i )->rho_i_custom( gTroom );
            }
            this->set_constant( MaterialProperty::rho_0, tRhoI / ( tRRR - 1.0 ) );
            this->set_constant( MaterialProperty::RRR, tRRR );
            this->create_tables();
            if ( mComputeTables )
            {
                this->populate_rho_database();
            }

            bool tIsMu0 = true ;
            for ( Metal * tComp : mComponents )
            {
                if ( ! tComp->have( MaterialProperty::mu ) )
                {
                    tIsMu0 = false ;
                    break ;
                }
                if ( std::abs( tComp->constant_property( MaterialProperty::mu ) - constant::mu0 ) > BELFEM_EPSILON )
                {
                    tIsMu0 = false ;
                    break ;
                }
            }
            if ( tIsMu0 )
            {
                this->set_constant( MaterialProperty::mu, constant::mu0 );
            }
        }

        void
        Alloy::update_volume_fractions( const real T )
        {
            Vector< real > & V = mVolumeFractions ;
            Vector< real > & Y = mMassFractions ;

            uint n = mComponents.size() ;
            for ( uint i=0; i<n; ++i )
            {
                V( i ) = Y( i ) / mComponents(i)->density( T );
            }
            V/=sum( V );
        }

        Metal *
        Alloy::create_component( const string & aName )
        {
            string tName = string_to_lower( aName );

            MaterialFactory tFactory;
            Material * tMaterial = tFactory.create_material( tName, BELFEM_QUIET_NAN, false );

            BELFEM_ERROR( tMaterial->type() == MaterialType::PureMetal , "%s is not a pure metal", aName.c_str() );

            Metal * aMetal = reinterpret_cast< Metal * >( tMaterial );

            real rho_0 = aMetal->constant_property( MaterialProperty::rho0_pure );

            real rho_i = aMetal->constant_property( MaterialProperty::rho_i_ref );
            real RRR = (rho_i + rho_0) / rho_0 ;
            aMetal->set_constant( MaterialProperty::rho_0, rho_0 );
            aMetal->set_constant( MaterialProperty::RRR, RRR );

            // we don't call set_RRR because this would initialize lookup tables
            // which we don't need
            //aMetal->set_RRR( RRR );

            return aMetal ;
        }

        void
        Alloy::create_tables()
        {
            const real dT = 4.0 ; // temperature step

            size_t tNumPoints = std::ceil( this->constant_property( MaterialProperty::T_max ) / dT )  + 1 ;

            Vector< real > tT( tNumPoints );
            Vector< real > tE( tNumPoints );
            Vector< real > tNu( tNumPoints );
            Vector< real > tAlpha( tNumPoints );
            Vector< real > tCp( tNumPoints );
            Vector< real > tRho( tNumPoints );
            Vector< real > tLambda( tNumPoints );

            real T = 0.0 ;

            uint n = mComponents.size() ;

            Vector< real > K( n );
            Vector< real > G( n );
            Vector< real > cp( n );

            Vector< real > rho( n );    // means electric resistivity
            Vector< real > w( n );

            Cell< Tensor< real > * > tTensors( n+8, nullptr );
            for ( index_t k=0; k<n+8; ++k )
            {
                tTensors( k ) = new Tensor< real > ( 3, 3, 3, 3 );
            }
            tensor::identity( *tTensors( n+6 ) );
            Vector< real > tWork( 72 );
            Vector< int_t >  tPivot( 36 );

            // compute elastic properties of individual components
            real Km = 0.0 ;
            real Gm = 0.0 ;

            real rho_0 = this->constant_property( MaterialProperty::rho_0 );

            const Vector< real > & v = mVolumeFractions ;

            real rho_i = 0.0 ;
            real w_i = 0.0 ;

            real tInvRefDensity = 0.0 ;
            for ( uint k=0; k<n; ++k )
            {
                tInvRefDensity += mMassFractions( k ) / mComponents( k )->density( gTroom );
            }
            real ref_density = 1.0 / tInvRefDensity ;

            this->set_constant( MaterialProperty::ref_density, ref_density );
            this->set_constant( MaterialProperty::T_ref_density, gTroom );

            for ( uint i=0; i<tNumPoints; ++i )
            {
                real K_voigt = 0.0 ;
                real G_voigt = 0.0 ;
                real alpha_voigt = 0.0 ;
                real K_reuss = 0.0 ;
                real G_reuss = 0.0 ;
                real alpha_reuss = 0.0 ;

                this->update_volume_fractions( T );

                for ( index_t k=0; k<n; ++k )
                {
                    real nu = mComponents( k )->nu( T );
                    real E  = mComponents( k )->E( T ) * 1e-9 ; // converting to GPa for better scaling
                    real alpha = mComponents( k )->alpha( T ) ;

                    cp( k ) = mComponents( k )->cp( T ); // specific heat

                    // bulk modulus
                    K(k)  = E / ( 3. * ( 1.0 - 2.0 * nu ) );

                    // shear modulus
                    G(k)  = E / ( 2.0 * ( 1.0 + nu ) );

                    // upper bound
                    K_voigt += v( k ) * K( k ) ;
                    G_voigt += v( k ) * G( k ) ;
                    alpha_voigt += v( k ) * alpha ;

                    // lower bound
                    K_reuss += v( k ) / K( k );
                    G_reuss += v( k ) / G( k );
                    alpha_reuss += v( k ) * alpha / K( k );

                    // populate tensors
                    tTensors( k )->fill_isotropic_elasticity( E, nu );

                }

                K_reuss = 1.0/K_reuss ;
                G_reuss = 1.0/G_reuss ;



                alpha_reuss *= K_reuss ;

                // now we can create an initial guess for the mixed property
                if ( i == 0 )
                {
                    Km = 0.5 * ( K_voigt + K_reuss );
                    Gm = 0.5 * ( G_voigt + G_reuss );
                }

                // iterate self-consistent elasticity properties
                this->self_consistent_KG( Km, Gm, tTensors, tWork, tPivot );

                real Em  = 9. * Km * Gm / ( 3. * Km + Gm );
                real num = ( 3. * Km - 2. * Gm ) / ( 6. * Km + 2. * Gm );


                tT( i ) = T ;
                tE( i ) = Em * 1e9 ;
                tNu( i ) = num ;

                if ( i == 0 )
                {
                    tAlpha( 0 ) = 0.0 ;
                    tCp( 0 ) = 0.0 ;
                    tRho( 0 ) = rho_0 ;
                    tLambda( 0 ) = 0.0 ;
                }
                else
                {
                    // Schapery with self-consistent Km
                    tAlpha( i ) = std::abs( K_reuss - K_voigt ) < BELFEM_EPSILON ? 0.5 * ( alpha_reuss + alpha_voigt ) :
                               alpha_voigt + ( 1./Km - 1./K_voigt) / ( 1./K_reuss - 1./K_voigt ) * ( alpha_reuss - alpha_voigt );

                    tCp( i ) = dot( mMassFractions, cp );


                    for ( index_t k=0; k<n; ++k )
                    {
                        rho( k ) = mComponents( k )->rho_i_custom( T ) ;
                        w( k ) = hust( mComponents( k )->lambda_coefficients(), T ) ;
                    }
                    rho_i = dot( v, rho ) ;
                    w_i = dot( v, w );

                    tRho( i ) = rho_0 + rho_i ;

                    real w0 = rho_0 / ( constant::L0 * T );
                    tLambda( i ) = 1.0 / ( w0 + w_i );
                }

                T+=dT;
            }


            for ( Tensor< real > * tTensor : tTensors )
            {
                delete tTensor ;
            }
            tTensors.clear();

            spline::SplineBC tStartBC = spline::SplineBC::Tangent ;
            spline::SplineBC tEndBC = spline::SplineBC::Parabolic  ;

            SpMatrix tA ;
            spline::create_helpmatrix( tNumPoints, dT, tA, tStartBC, tEndBC );

            this->set_spline( MaterialProperty::E, new Spline( tT, tE, tA, tStartBC, tEndBC, 0.0 ) );

            this->set_spline( MaterialProperty::nu, new Spline( tT, tNu, tA, tStartBC, tEndBC, 0.0 ) );

            real gamma = 0.0 ;
            for ( uint k=0; k<n; ++k )
            {
                gamma += mMassFractions( k ) * mComponents( k )->constant_property( MaterialProperty::gamma );
            }
            this->set_constant( MaterialProperty::gamma, gamma );
            this->set_spline( MaterialProperty::cp, new Spline( tT, tCp, tA, tStartBC, tEndBC, gamma ) );

            this->set_spline( MaterialProperty::rho, new Spline( tT, tRho, tA, tStartBC, tEndBC, 0.0 ) );

            tStartBC = spline::SplineBC::Parabolic  ;
            spline::create_helpmatrix( tNumPoints, dT, tA, tStartBC, tEndBC );
            this->set_spline( MaterialProperty::alpha, new Spline( tT, tAlpha, tA, tStartBC, tEndBC ) );

            this->spline( MaterialProperty::alpha )->create_integral( this->constant_property( MaterialProperty::T_ref_density ), 0.0 );
            this->set_spline( MaterialProperty::lambda, new Spline( tT, tLambda, tA, tStartBC, tEndBC ) );

            this->set_custom( MaterialProperty::density );
        }

        void
        Alloy::populate_rho_database()
        {
            if ( mRhoData != nullptr )
            {
                delete mRhoData ;
            }

            string tFile = sprint( "%s_RRR%u.hdf5" ,
                this->label().c_str(),
                ( uint ) this->constant_property( MaterialProperty::RRR ) );


            // see Metal::populate_rho_database: a cached file carries no format
            // version, so probe for the RRR marker and rebuild an older one rather
            // than failing inside load_rho_database. The verdict must be the same
            // on every rank, hence probe-on-master plus broadcast.
            bool tUsable = file_exists( tFile );

            if ( tUsable && comm_rank() == 0 )
            {
                tUsable = rho_database_is_current( tFile );

                if ( ! tUsable )
                {
                    message( InfoLevel::Default,
                        "    Warning: %s predates the current rho-database format "
                        "( no RRR entry ) and is being rebuilt",
                        tFile.c_str() );
                }
            }
            if ( comm_size() > 1 )
            {
                broadcast( tUsable );
            }

            if ( ! tUsable )
            {
                // all ranks build in lockstep; the workers hold empty grid
                // copies and receive the projected values inside the Database
                // constructor -- see Metal::populate_rho_database
                this->populate_rho_database_serial();

                this->save_rho_database( tFile );
            }
            else
            {
                this->load_rho_database( tFile );
            }

            this->set_have( MaterialProperty::rho );
            this->set_have( MaterialProperty::lambda );

            // the field database is valid from here on: expose the
            // rho( T, B, beta ) and lambda( T, B, beta ) overrides to the
            // dependency-based dispatch ( calculator::MaxwellData keys on
            // depends( . , normB ), not on the material type ), same
            // pattern as Copper
            this->set_dependency( MaterialProperty::rho, MaterialDependency::T );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::normB );
            this->set_dependency( MaterialProperty::rho, MaterialDependency::angleBxJ );

            this->set_dependency( MaterialProperty::lambda, MaterialDependency::T );
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::normB );
            this->set_dependency( MaterialProperty::lambda, MaterialDependency::angleBxJ );

            comm_barrier() ;
        }

        void
        Alloy::populate_rho_database_serial()
        {
            Mesh * tMesh ;

            // create a tensor mesh
            tMesh = create_database_mesh();

            mDatabaseTmin = tMesh->tensorconf()->min( 0 ) ;
            mDatabaseTmax = tMesh->tensorconf()->max( 0 ) ;
            mDatabaseBmin = std::pow( 10, tMesh->tensorconf()->min( 1 )  );
            mDatabaseBmax = std::pow( 10, tMesh->tensorconf()->max( 1 )  );


            Vector< real > & tRho = tMesh->create_field( "rho" );
            this->populate_rho_database( tMesh, tRho );
            string tDisplayLabel = sprint( "%s RRR %u",
                this->label().c_str(),
                ( uint ) this->constant_property( MaterialProperty::RRR ) );
            mRhoData = new Database( tMesh, "rho", true, tDisplayLabel );

            // the Database keeps its own copy of the values and config;
            // the work grid is ours to free
            delete tMesh ;
        }

        void
        Alloy::populate_rho_database( Mesh * aMesh, Vector< real > & aRho )
        {
            proc_t tRank = comm_rank() ;

            real T0 = mDatabaseTmin ;

            Cell< mesh::Node * > & tNodes = aMesh->nodes();

            // populate dataset
            index_t tCount = 0 ;

            // number of components
            uint n = mComponents.size() ;
            Vector< real > rho_i( n );
            Vector< real > rho_0( n );
            Vector< real > rho_ref( n );
            Vector< real > delta_rho( n );
            Vector< real > w_i( n );

            this->update_volume_fractions( T0 );

            const Vector< real > & v = mVolumeFractions ;
            for ( uint k=0; k<n; ++k )
            {
                rho_0( k ) = mComponents( k ) ->constant_property( MaterialProperty::rho_0 );
                rho_i( k ) = mComponents( k ) ->constant_property( MaterialProperty::rho_i_ref );
            }
            rho_ref = rho_0 + rho_i ;
            for ( uint k=0; k<n; ++k )
            {
                rho_i( k ) = mComponents( k ) ->rho_i_custom( T0 );
            }

            real rho_0K = this->constant_property( MaterialProperty::rho_0 );
            real rho_0T  = rho_0K + dot( v, rho_i );

            real rho ;

            for ( mesh::Node * tNode : tNodes )
            {
                if ( tNode->owner() != tRank ) continue ;

                real T = tNode->x();

                if ( T != T0 )
                {
                    this->update_volume_fractions( T );

                    // reference at 0 T
                    for ( uint k=0; k<n; ++k )
                    {
                        rho_i( k ) = mComponents( k ) ->rho_i_custom( T );
                    }

                    rho_0T = rho_0K + dot( v, rho_i );
                    T0 = T ;
                }

                real B = std::pow( 10., tNode->y() );
                real beta = tNode->z();

                for ( uint k=0; k<n; ++k )
                {
                    real S = rho_ref( k ) / ( rho_0K + rho_i( k ));
                    delta_rho( k ) = mComponents( k ) ->kohler( B, S, beta );
                }

                // projected value for temperature, field and angle
                if ( B < 0.001 )
                {
                    rho = rho_0T ;
                }
                else
                {
                    rho = rho_0T * ( 1. + dot( v, delta_rho ) ) ;
                }

                // store data in field
                aRho( tCount++ )    = std::log(rho) ;
            }
        }

        void
        Alloy::self_consistent_KG( real & Km, real & Gm, Cell< Tensor< real > * > & Tensors, Vector< real > & Work, Vector< int_t > & Pivot )
        {
            // relaxation factor
            const real omega = 0.99 ;


            uint n = mComponents.size() ;

            // elasticity
            Tensor< real > & C = *Tensors( n   ) ;

            // compliance
            Tensor< real > & M = *Tensors( n+1 ) ;

            // eshelby
            Tensor< real > & S = *Tensors( n+2 ) ;

            // hill polarization
            Tensor< real > & P = *Tensors( n+3 ) ;

            // influence tensor
            Tensor< real > & A = *Tensors( n+4 ) ;

            // compliance from previous timestep
            Tensor< real > & C0 = *Tensors( n+5 ) ;

            // identity tensor
            Tensor< real > & I = *Tensors( n+6 ) ;

            // work tensor
            Tensor< real > & B = *Tensors( n+7 ) ;

            // volume fractions
            const Vector< real > & v = mVolumeFractions ;

            uint tCount = 0 ;

            real f = 1 ;

            while ( abs( f ) > 1e-8 )
            {
                real K0 = Km ;
                real G0 = Gm ;

                // mean elasticity tensor
                C.fill( Km, Gm );

                // mean compliance tensor
                M.fill( 1./(9.*Km), 1./(4.*Gm) );

                // compute the factors for the eshelby tensor
                // Qu, page 53, Eq. (4.3.40)
                real gamma = Km / ( 3. * Km + 4. * Gm );
                real delta = 3. * ( Km + 2. * Gm ) / ( 5. * ( 3. * Km + 4. * Gm ) );

                S.fill( gamma, delta );

                // compute the hill tensor
                P = S % M ;

                C0 = C ;
                C.fill( 0.0 );


                for ( index_t k=0; k<n; ++k )
                {
                    Tensor< real > & Ck = *Tensors( k );

                    // compute the compliance tensor
                    B = Ck - C0 ;
                    A = P % B ;
                    A += I ;
                    tensor::invert_symmetric( A, Work, Pivot );

                    B = Ck % A ;
                    B *= v( k ) ;
                    C += B ;
                }

                // recover new mixing values
                real G1 = C(1,0,1,0);
                real K1 = C(0,0,0,0) - 4./3.*G1 ;

                Gm = ( 1. - omega ) * G0 + omega * G1 ;
                Km = ( 1. - omega ) * K0 + omega * K1 ;

                real kerr = ( Km - K0 ) / K0 ;
                real gerr = ( Gm - G0 ) / G0 ;
                f = std::sqrt( kerr*kerr + gerr*gerr ) ;

                BELFEM_ERROR( ++tCount < 1000, "Failed to converge" );
            }
        }


        void
        Alloy::save_rho_database( const std::string & aPath )
        {
            if ( comm_rank() != 0 ) return ;

            HDF5 tFile( aPath, FileMode::NEW );
            tFile.save_data( "label" , this->label() );
            tFile.save_data( "RRR",  this->constant_property( MaterialProperty::RRR ) );
            tFile.save_data( "Tmin", mDatabaseTmin );
            tFile.save_data( "Tmax", mDatabaseTmax );
            tFile.save_data( "Bmin", mDatabaseBmin );
            tFile.save_data( "Bmax", mDatabaseBmax );

            tFile.create_group( "rho" );
            mRhoData->save( tFile.active_group() );
            tFile.close_active_group();

            // ship the reference reader with the table, if it is available
            material::embed_python_guide( tFile );

            tFile.close();
        }

        void
        Alloy::load_rho_database( const std::string & aPath )
        {
            if ( comm_rank() == 0 )
            {
                HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

                string tLabel ;
                tFile.load_data( "label", tLabel );

                BELFEM_ERROR( tLabel == this->label(), "Found material %s but expect %s",
                    tLabel.c_str(), this->label().c_str() );

                real tRRR ;
                tFile.load_data( "RRR", tRRR );
                BELFEM_ERROR( tRRR == this->constant_property( MaterialProperty::RRR ), "Found RRR %g for material %s but expect %g",
                   ( double ) tRRR, this->label().c_str(), ( double ) this->constant_property( MaterialProperty::RRR ) );

                tFile.load_data( "Tmin", mDatabaseTmin );
                tFile.load_data( "Tmax", mDatabaseTmax );
                tFile.load_data( "Bmin", mDatabaseBmin );
                tFile.load_data( "Bmax", mDatabaseBmax );

                tFile.select_group( "rho" );
                if ( mRhoData != nullptr ) delete mRhoData ;
                mRhoData = new Database( tFile.active_group() , "rho");
                tFile.close_active_group();
                tFile.close();

                Vector< real > tLimits( 4 );
                tLimits(0) = mDatabaseTmin ;
                tLimits(1) = mDatabaseTmax ;
                tLimits(2) = mDatabaseBmin ;
                tLimits(3) = mDatabaseBmax ;
                share( tLimits );
            }
            else
            {
                mRhoData    = new Database( ( hid_t ) 0 , "rho");

                Vector< real > tLimits( 4 );
                receive( tLimits );
                mDatabaseTmin = tLimits(0) ;
                mDatabaseTmax = tLimits(1) ;
                mDatabaseBmin = tLimits(2) ;
                mDatabaseBmax = tLimits(3) ;
            }
        }
    }
}
