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

#ifndef BELFEM_CL_MATERIAL_ALLOY_HPP
#define BELFEM_CL_MATERIAL_ALLOY_HPP

#include "cl_Material_SplineLookupTable.hpp"

#include "cl_Material_Metal.hpp"
#include "cl_Tensor.hpp"

namespace belfem
{
    namespace material
    {

        class Alloy : public SplineLookupTable
        {
            Vector< real > mMolarMasses ;
            Vector< real > mMolarFractions ;
            Vector< real > mVolumeFractions ;
            Vector< real > mMassFractions ;

            Cell< Metal * > mComponents ;

            Database * mRhoData = nullptr ;

            real mDatabaseTmin = BELFEM_QUIET_NAN ;
            real mDatabaseTmax = BELFEM_QUIET_NAN ;
            real mDatabaseBmin = BELFEM_QUIET_NAN ;
            real mDatabaseBmax = BELFEM_QUIET_NAN ;

            const real mInvLog10 = 1.0 / std::log(10.0);

            bool mComputeTables = true ;
        public:

            Alloy( const string & aName,
                   const bool aBuildTables = true );

            Alloy( const string & aName,
                   const Cell< std::pair< string, real > > & aComponents,
                   const real aRRR = BELFEM_QUIET_NAN,
                   const bool aBuildTables = true );

            ~Alloy() override;

            void
            set_components(
                const Cell< string > & aComponents,
                const Vector< real > & aMassFractions,
                const real RRR = BELFEM_QUIET_NAN  );

            // Bring base class single-parameter versions into scope
            // (needed because C++ name hiding hides Material::rho(T) when we override rho(T,B,beta))
            using Material::rho;
            using Material::lambda;
            using Material::drhodT;
            using Material::dlambdadT;

            /**
             * @brief Electrical resistivity with magnetoresistance (Kohler's rule)
             * @param T Temperature [K]
             * @param B Magnetic field magnitude [T]
             * @param beta Angle between current and field [rad]
             * @return Electrical resistivity [Ω·m]
             *
             * Applies Kohler's rule: ρ(B,T) = ρ₀(T) · [1 + K(B·S, β)]
             * where S = ρ_ref/ρ₀(T) is the similarity parameter.
             */
            real
            rho( const real T, const real B, const real beta ) const override ;

            real
            drhodT( const real T, const real B, const real beta ) const override ;

            real
            drhodB( const real T, const real B, const real beta ) const override ;

            real
            drhodbeta( const real T, const real B, const real beta ) const override ;

            /**
             * @brief Thermal conductivity with magnetoresistance
             * @param T Temperature [K]
             * @param B Magnetic field magnitude [T]
             * @param beta Angle between current and field [rad]

             * @return Thermal conductivity [W/(m·K)]
             *
             * Uses Wiedemann-Franz relation corrected for magnetoresistance.
             */
            real
            lambda( const real T, const real B, const real beta ) const override ;

            real
            dlambdadT( const real T, const real B, const real beta ) const override ;

        private:

            void
            set_components( const Cell< std::pair< string, real > > & aComponents );

            void
            update_volume_fractions( const real T );

            Metal *
            create_component( const string & aName );

            void
            delete_components();

            void
            create_tables();

            void
            self_consistent_KG( real & Km, real & Gm, Cell< Tensor< real > * > & Tensors, Vector< real > & Work, Vector< int_t > & Pivot );

            void
            populate_rho_database();

            void
            populate_rho_database( Mesh * aMesh, Vector< real > & aRho );

            void
            populate_rho_database_serial() ;

            void
            save_rho_database( const std::string & aPath );

            void
            load_rho_database( const std::string & aPath );

        };

        inline real
        Alloy::rho( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real log10B = std::log(std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            return std::exp(mRhoData->evaluate( theta, log10B, angle )) ;
        }

        inline real
        Alloy::drhodT( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real log10B = std::log(std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y    = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydT = mRhoData->evaluate_derivx( theta, log10B, angle ) ;
            return std::exp( y ) * dydT ;
        }

        // field derivatives: same table chain rules as Metal::drhodB_table /
        // drhodbeta_table ( rho = exp(y), y on the (T, log10B, angle) grid )
        inline real
        Alloy::drhodB( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real Bc = std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) ;

            real log10B = std::log( Bc )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y       = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydlogB = mRhoData->evaluate_derivy( theta, log10B, angle ) ;

            return std::exp( y ) * dydlogB * mInvLog10 / Bc ;
        }

        inline real
        Alloy::drhodbeta( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real Bc = std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) ;

            real log10B = std::log( Bc )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y       = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydbeta = mRhoData->evaluate_derivz( theta, log10B, angle ) ;

            return std::exp( y ) * dydbeta ;
        }

        inline real
        Alloy::lambda( real T, const real B, const real beta ) const
        {
            return this->lambda( T ) * this->rho( T ) / this->rho( T, B, beta );
        }

        inline real
        Alloy::dlambdadT( real T, const real B, const real beta ) const
        {
            real a = this->lambda( T ) ;
            real b = this->rho( T ) ;
            real c = this->rho( T, B, beta );

            real da = this->dlambdadT( T );
            real db = this->drhodT( T );
            real dc = this->drhodT( T, B, beta );

            return ( c * ( da*b + a * db ) - a * b * dc ) /( c * c );
        }

    }
}
#endif //BELFEM_CL_MATERIAL_ALLOY_HPP