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

#ifndef BELFEM_CL_GT_REFGAS_HPP
#define BELFEM_CL_GT_REFGAS_HPP

#include "typedefs.hpp"

#include "cl_Map.hpp"
#include "cl_Cell.hpp"

#include "cl_Spline.hpp"

#include "cl_GT_GasData.hpp"
#include "cl_GT_HeatPoly.hpp"
#include "cl_GT_TransportPoly.hpp"

namespace belfem
{
    namespace gastables
    {

        class InputThermo;
        class InputTransport;
        class RefGasFactory;

        /**
         * how a RefGas evaluates its properties
         *
         * POLY   : directly from the piecewise NASA-9 polynomials, exact but
         *          with a linear search for the interval that holds T
         * SPLINE : from cubic splines sampled off those polynomials, which is
         *          what the factory leaves every gas in
         */
        enum class RefGasMode
        {
            POLY      = 0,
            SPLINE    = 1,
            UNDEFINED = 2
        };

        /**
         * A single chemical species with its caloric and transport properties.
         *
         * The caloric side is the NASA-9 formulation of CEA ( NASA RP-1311 ):
         * a set of temperature intervals, each carrying a seven coefficient
         * polynomial for cp/R plus two integration constants for H and S.
         * Transport uses the CEA correlation ln( eta ) = A ln T + B/T + C/T^2 + D
         * over its own intervals. Both are read from the shipped tables by
         * RefGasFactory, which also synthesizes the intervals the tables do not
         * cover: glue polynomials across interval junctions, a cryogenic
         * extrapolation below the lowest interval, and viscosity and thermal
         * conductivity from the Lucas and Chung correlations for species with
         * critical data but no transport record.
         *
         * Everything is molar unless the accessor is lower case: H() is J/mol
         * and h() is J/kg, the same convention as Cp()/cp(). Temperatures are
         * in K throughout.
         *
         * The object owns its polynomials and deletes them; mData is held by
         * value. Construction is through RefGasFactory, not directly, because
         * a usable object needs the data files.
         *
         * @ingroup grp_physics_gastables
         * @see @ref physics_gastables_gastables_usage_guide
         */
        class RefGas
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // flag telling if this is a noble gas
            const bool mAmNoble;

            GasData mData;

            const string & mLabel;

            // heat polynomials
            Cell< HeatPoly * > mHeatPolys;

            // transport polynomials
            Cell< TransportPoly * > mViscosityPolys;
            Cell< TransportPoly * > mConductivityPolys;

            // flag telling if this gas has thermo data
            bool mHaveThermo = false;

            // flag telling if this gas has viscosity data
            bool mHaveViscosity = false;

            // flag telling if this gas has conductivity data
            bool mHaveConductivity = false;

            // flag telling if this is a liquid polynomial according to CEA
            bool mLiquidFlag = false;

            // flag telling if intermediate polynomials have been calculated
            bool mFinalizedFlag = false;

            bool mHaveComponents = false;

            Spline mHeatSpline;
            Spline mViscositySpline;
            Spline mConductivitySpline;

            real
            ( RefGas:: * mFunctionCp )          ( const real T ) const ;

            real
            ( RefGas:: * mFunctiondCpdT )       ( const real T ) const ;

            real
            ( RefGas:: * mFunctiond2CpdT2 )     ( const real T ) const ;

            real
            ( RefGas:: * mFunctionH )           ( const real T ) const ;

            real
            ( RefGas:: * mFunctionS )           ( const real T ) const ;

            real
            ( RefGas:: * mFunctiondSdT )        ( const real T ) const ;

            real
            ( RefGas:: * mFunctionMu )          ( const real T ) const ;

            real
            ( RefGas:: * mFunctiondMudT )       ( const real T ) const ;

            real
            ( RefGas:: * mFunctiond2MudT2 )     ( const real T ) const ;

            real
            ( RefGas:: * mFunctionLambda )      ( const real T ) const ;

            real
            ( RefGas:: * mFunctiondLambdadT )   ( const real T ) const ;

            real
            ( RefGas:: * mFunctiond2LambdadT2 ) ( const real T ) const ;

            friend InputThermo;
            friend InputTransport;
            friend RefGasFactory;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            RefGas( const string & aLabel );

//------------------------------------------------------------------------------

            ~RefGas();

//------------------------------------------------------------------------------

            // owns raw polynomial pointers and holds a reference into mData
            RefGas( const RefGas & ) = delete;
            RefGas & operator=( const RefGas & ) = delete;

//------------------------------------------------------------------------------

            inline const string &
            label() const;

//------------------------------------------------------------------------------

            /**
             * return the composition if the element exists
             */
            inline real
            component_multiplicity( const string & aLabel ) const;

//------------------------------------------------------------------------------

            /**
             * molar mass in kg/Mol
             */
            inline const real &
            M() const;

//------------------------------------------------------------------------------

            /**
             * reference formation enthalpy
             */
            inline real
            reference_formation_enthalpy() const;

//------------------------------------------------------------------------------

            inline bool
            is_liquid() const;
//------------------------------------------------------------------------------

            /**
             * @name Caloric properties
             *
             * Upper case is molar, lower case is mass specific: Cp is
             * J/(mol K) and cp is J/(kg K), H is J/mol and h is J/kg. The two
             * differ by the molar mass and nothing else. Temperature is in K.
             *
             * The enthalpy scale is set by the factory, not by the source
             * record: it discards the record's integration constant and shifts
             * every interval so that
             *
             *     H( 298.15 K ) = dHf( 298.15 K ) + [ H( 298.15 K ) - H( 0 K ) ]
             *
             * with the formation enthalpy and the sensible term both taken from
             * thermo.inp. Equivalently, H( 0 K ) is the formation enthalpy at
             * 298.15 K. That is a hybrid of the two usual conventions - CEA
             * places the formation enthalpy at 298.15 K, a 0 K scale places it
             * at 0 K - and it is deliberate.
             *
             * No consumer in the tree can see the difference. The combustion
             * solver re-references to the CEA convention itself, subtracting
             * h( 298.15 ) and adding the formation enthalpy back
             * ( cl_CN_Scheme.cpp in the nonfree tree ), and the flow routines only ever evaluate
             * enthalpy differences, where the constant cancels. It would matter
             * only if an absolute enthalpy were compared against an outside
             * table.
             *
             * Entropy is absolute, anchored on the standard state entropy of
             * the record.
             *
             * Which representation answers these calls depends on the mode -
             * the polynomials themselves in POLY, the sampled splines in
             * SPLINE. The factory leaves every gas in SPLINE.
             * @{
             */

            real
            Cp( const real T ) const;

//------------------------------------------------------------------------------

            real
            H( const real T ) const;

//------------------------------------------------------------------------------

            real
            S( const real T ) const;

//------------------------------------------------------------------------------

            real
            dSdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            dCpdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            cp( const real T ) const;

//------------------------------------------------------------------------------

            real
            h( const real T ) const;

//------------------------------------------------------------------------------

            real
            h_ref() const;

//------------------------------------------------------------------------------

            real
            H_ref() const;

//------------------------------------------------------------------------------

            real
            s( const real T ) const;

//------------------------------------------------------------------------------

            real
            dcpdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            d2cpdT2( const real T ) const;

            /** @} */

//------------------------------------------------------------------------------

            /**
             * @name Transport properties
             *
             * Dynamic viscosity in Pa s and thermal conductivity in W/(m K),
             * both SI - the tabulated CEA coefficients are micropoise and
             * microwatt per centimetre kelvin and are converted on the way out.
             * Temperature is in K.
             *
             * These describe the dilute gas. Pressure dependence is not part of
             * this class; the gasmodels layer adds it.
             *
             * A species with no transport record returns zero rather than
             * raising, except where the factory was able to synthesize
             * viscosity and conductivity from the critical point.
             * @{
             */

            real
            mu( const real T ) const;

//------------------------------------------------------------------------------

            real
            dmudT( const real T ) const;

//------------------------------------------------------------------------------

            real
            d2mudT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            lambda( const real T ) const;

//------------------------------------------------------------------------------

            real
            dlambdadT( const real T ) const;

//------------------------------------------------------------------------------

            real
            d2lambdadT2( const real T ) const;

            /** @} */

//------------------------------------------------------------------------------

            /**
             * @name Data availability
             *
             * Which of the source tables actually had a record for this
             * species. Callers use these to decide whether a property is worth
             * asking for - the accessors of a missing set return zero rather
             * than raising, so a silent zero is otherwise indistinguishable
             * from a real one.
             * @{
             */

            inline bool
            has_thermo() const;

//------------------------------------------------------------------------------

            inline bool
            has_conductivity() const;

//------------------------------------------------------------------------------

            inline bool
            has_viscosity() const;


//------------------------------------------------------------------------------

            inline bool
            has_components() const;

//------------------------------------------------------------------------------

            inline bool
            is_noble() const;

            /** @} */

//------------------------------------------------------------------------------

            /**
             * Switch between evaluating the CEA polynomials directly and
             * evaluating the splines sampled off them. The splines are faster
             * because they avoid the interval search, at interpolation
             * accuracy; the factory leaves every gas in SPLINE.
             *
             * Only rebinds function pointers - both modes describe the same
             * ideal gas, so this is not a change of physics.
             */
            void
            set_mode( const RefGasMode & aMode );

//------------------------------------------------------------------------------

            /**
             * expose data object
             */
            inline GasData *
            data();

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * expose data object ( const version )
             */
            inline const GasData *
            data() const;


//------------------------------------------------------------------------------

            /**
             * expose heat spline
             */
            inline Spline *
            heat_spline();

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * expose heat spline ( const version )
             */
            inline const Spline *
            heat_spline() const;

//------------------------------------------------------------------------------

            /**
             * expose viscosity spline
             */
            inline Spline *
            viscosity_spline();

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * expose viscosity spline ( const version )
             */
            inline const Spline *
            viscosity_spline() const;

//------------------------------------------------------------------------------

            /**
             * expose conductivity spline
             */
            inline Spline *
            conductivity_spline();

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * expose conductivity spline ( const version )
             */
            inline const Spline *
            conductivity_spline() const;

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * to be called by RefGasFactory once all polynomials have been read
             */
            void
            finalize();

//------------------------------------------------------------------------------

            void
            create_splines(
                    const Vector< real > & T,
                          SpMatrix       & aHelpMatrix );


//------------------------------------------------------------------------------

            /**
             * to be called by InputThermo
             */
            void
            add_component( const string & aLabel, const real aValue );

//------------------------------------------------------------------------------

            /**
             * set the molar mass
             * @param[in] aMolarMass in kg/Mol, not g/Mol!
             */
            void
            set_molar_mass( const real aMolarMass );

//------------------------------------------------------------------------------

            /**
             * set the formation enthalpy in J/mol
             */
            void
            set_reference_formation_enthalpy( const real aDeltaHf );

//------------------------------------------------------------------------------

            /**
             * set the reference enthalpy H(298.15 K) - H(0 K) in J/mol,
             * as read from the thermo.inp record
             */
            void
            set_reference_enthalpy( const real aHref );

//------------------------------------------------------------------------------

            void
            set_liquid_flag();


//------------------------------------------------------------------------------

            void
            set_component_flag();

//------------------------------------------------------------------------------

            void
            unset_liquid_flag();

//------------------------------------------------------------------------------

            void
            add_heat_poly( HeatPoly * aHeatPoly );

//------------------------------------------------------------------------------

            void
            add_transport_poly( TransportPoly * aTransportPoly );

//------------------------------------------------------------------------------

            void
            delete_heat_polys();

//------------------------------------------------------------------------------

            void
            delete_transport_polys();

//------------------------------------------------------------------------------

            void
            create_glue_polys_heat( const uint & aNumberOfOriginalPolynomials );

//------------------------------------------------------------------------------

            void
            create_cryo_poly_heat();

//------------------------------------------------------------------------------

            void
            create_hot_poly_heat( const uint & aNumberOfOriginalPolynomials );

//------------------------------------------------------------------------------

            void
            fix_reference_points(  const uint & aStart, const uint & aEnd );

//------------------------------------------------------------------------------

            const HeatPoly *
            find_heat_poly( const real T ) const;

//------------------------------------------------------------------------------

            const TransportPoly *
            find_viscosity_poly( const real T ) const;

//------------------------------------------------------------------------------

            const TransportPoly *
            find_conductivity_poly( const real T ) const;

//------------------------------------------------------------------------------

            void
            finalize_transport();

//------------------------------------------------------------------------------

            void
            create_cryo_poly_transport( Cell< TransportPoly * > & aPolys );

//------------------------------------------------------------------------------

            void
            create_glue_polys_transport(
                    Cell< TransportPoly * > & aPolys,
                    const uint & aNumberOfOriginalPolynomials );

//------------------------------------------------------------------------------

            void
            create_hot_poly_transport( Cell< TransportPoly * > & aPolys,
                                       const uint & aNumberOfOriginalPolynomials );

//------------------------------------------------------------------------------

            void
            finalize_thermo();

//------------------------------------------------------------------------------

            real
            zero( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_Cp( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_dCpdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_d2CpdT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_H( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_S( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_dSdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_Mu( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_dMudT( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_d2MudT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_Lambda( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_dLambdadT( const real T ) const;

//------------------------------------------------------------------------------

            real
            poly_d2LambdadT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_Cp( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_dCpdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_H( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_S( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_dSdT( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_Mu( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_dMudT( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_d2MudT2( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_Lambda( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_dLambdadT( const real T ) const;

//------------------------------------------------------------------------------

            real
            spline_d2LambdadT2( const real T ) const;

//------------------------------------------------------------------------------

            void
            fix_switches();

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------

        const string &
        RefGas::label() const
        {
            return mData.label();
        }

// ------------------------------------------------------------------------------

        GasData *
        RefGas::data()
        {
            return &mData;
        }

// ------------------------------------------------------------------------------

        const GasData *
        RefGas::data() const
        {
            return & mData;
        }

// ------------------------------------------------------------------------------

        real
        RefGas::component_multiplicity( const string & aLabel ) const
        {
           return mData.component_multiplicity( aLabel );
        }

// ------------------------------------------------------------------------------

        bool
        RefGas::is_liquid() const
        {
            return mLiquidFlag;
        }

// ------------------------------------------------------------------------------

        const real &
        RefGas::M() const
        {
            return mData.M();
        }

// ------------------------------------------------------------------------------

        real
        RefGas::reference_formation_enthalpy() const
        {
            return mData.Hf();
        }

// ------------------------------------------------------------------------------

        bool
        RefGas::has_thermo() const
        {
            return mHaveThermo;
        }

//------------------------------------------------------------------------------

        bool
        RefGas::has_conductivity() const
        {
            return mHaveConductivity;
        }

//------------------------------------------------------------------------------

        bool
        RefGas::has_viscosity() const
        {
            return mHaveViscosity ;
        }


//------------------------------------------------------------------------------

        bool
        RefGas::has_components() const
        {
            return mHaveComponents;
        }

//------------------------------------------------------------------------------

        bool
        RefGas::is_noble() const
        {
            return mAmNoble;
        }

//------------------------------------------------------------------------------

        Spline *
        RefGas::heat_spline()
        {
            return & mHeatSpline;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Spline *
        RefGas::heat_spline() const
        {
            return & mHeatSpline;
        }

//------------------------------------------------------------------------------

        Spline *
        RefGas::viscosity_spline()
        {
            return & mViscositySpline;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Spline *
        RefGas::viscosity_spline() const
        {
            return & mViscositySpline;
        }

//------------------------------------------------------------------------------

        Spline *
        RefGas::conductivity_spline()
        {
            return & mConductivitySpline;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const Spline *
        RefGas::conductivity_spline() const
        {
            return & mConductivitySpline;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GT_THERMO_HPP
