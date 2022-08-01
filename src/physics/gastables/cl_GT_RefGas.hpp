//
// Created by Christian Messe on 2019-08-18.
//

#ifndef BELFEM_CL_GT_THERMO_HPP
#define BELFEM_CL_GT_THERMO_HPP

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

        enum class RefGasMode
        {
            POLY      = 0,
            SPLINE    = 1,
            UNDEFIED  = 2
        };

        // a cea thermo object
        class RefGas
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            const bool mAmNoble;

            GasData mData;

            const string & mLabel;

            // cell with elements this gas
            Cell< string > mElements;

            // map containing the composition
            Map< string, real > mComposition;

            // heat polynomials
            Cell< HeatPoly * > mHeatPolys;

            // transport polynomials
            Cell< TransportPoly * > mViscosityPolys;
            Cell< TransportPoly * > mConductivityPolys;

            // flag telling if this gas has thermo data
            bool mHaveThermo = false;

            // flag tellinf if this gas has viscosity data
            bool mHaveViscosity = false;

            // flag tellinf if this gas has viscosity data
            bool mHaveConductivity = false;

            // flag telling if this is a liquid polynomial according to CEA
            bool mLiquidFlag = false;

            // flag telling if intermediate polynomials have been calculated
            bool mFinalizedFlag = false;

            // flag telling if data exists in crthermo.inp and crtrans.inp
            bool mHaveCryoThermo = false;
            bool mHaveCryoTransport = false;

            bool mHaveComponents = false;

            // flag telling if this is a noble gas


            Spline mHeatSpline;
            Spline mViscositySpline;
            Spline mConductivitySpline;

            RefGasMode mMode;

            real
            ( RefGas:: * mFunctionCp )          ( const real aT );

            real
            ( RefGas:: * mFunctiondCpdT )       ( const real aT );

            real
            ( RefGas:: * mFunctiond2CpdT2 )     ( const real aT );

            real
            ( RefGas:: * mFunctionH )           ( const real aT );

            real
            ( RefGas:: * mFunctionS )           ( const real aT );

            real
            ( RefGas:: * mFunctiondSdT )        ( const real aT );

            real
            ( RefGas:: * mFunctionMu )          ( const real aT );

            real
            ( RefGas:: * mFunctiondMudT )       ( const real aT );

            real
            ( RefGas:: * mFunctiond2MudT2 )     ( const real aT );

            real
            ( RefGas:: * mFunctionLambda )      ( const real aT );

            real
            ( RefGas:: * mFunctiondLambdadT )   ( const real aT );

            real
            ( RefGas:: * mFunctiond2LambdadT2 ) ( const real aT );

            real
            ( RefGas:: * mFunctionZero) ( const real aT );

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

            real
            Cp( const real aT );

//------------------------------------------------------------------------------

            real
            H( const real aT );

//------------------------------------------------------------------------------

            real
            S( const real aT );

//------------------------------------------------------------------------------

            real
            dSdT( const real aT );

//------------------------------------------------------------------------------

            real
            dCpdT( const real aT );

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real aT );

//------------------------------------------------------------------------------

            real
            cp( const real aT );

//------------------------------------------------------------------------------

            real
            h( const real aT );

//------------------------------------------------------------------------------

            real
            h_ref() const;

//------------------------------------------------------------------------------

            real
            H_ref() const;

//------------------------------------------------------------------------------

            real
            s( const real aT );

//------------------------------------------------------------------------------

            real
            dcpdT( const real aT );

//------------------------------------------------------------------------------

            real
            d2cpdT2( const real aT );

//------------------------------------------------------------------------------

            real
            mu( const real aT );

//------------------------------------------------------------------------------

            real
            dmudT( const real aT );

//------------------------------------------------------------------------------

            real
            d2mudT2( const real aT );

//------------------------------------------------------------------------------

            real
            lambda( const real aT );

//------------------------------------------------------------------------------

            real
            dlambdadT( const real aT );

//------------------------------------------------------------------------------

            real
            d2lambdadT2( const real aT );

//------------------------------------------------------------------------------

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
            has_cryo_thermo() const;

//------------------------------------------------------------------------------

            inline bool
            has_cryo_transport() const;

//------------------------------------------------------------------------------

            inline bool
            has_components() const;

//------------------------------------------------------------------------------

            inline bool
            is_noble() const;

//------------------------------------------------------------------------------

            /**
             * Use either original CEA equations or Splines.
             * Later ones are faster.
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

//------------------------------------------------------------------------------

            /**
             * expose viscosity spline
             */
            inline Spline *
            viscosity_spline();

//------------------------------------------------------------------------------

            /**
             * expose conductivity spline
             */
            inline Spline *
            conductivity_spline();

//------------------------------------------------------------------------------

            // flag telling if this gas is marked for reference enthalpy
            // calculation
            void
            set_href_flag();

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * to be called by constructor
             */
            void
            finalize();

//------------------------------------------------------------------------------

            void
            create_splines(
                    const Vector< real > & aT,
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
             * set the formation enthalpy in J/mol
             */
            void
            set_reference_enthalpy( const real aHref );

//------------------------------------------------------------------------------

            void
            set_liquid_flag();

//------------------------------------------------------------------------------

            void
            set_cryo_transport_flag();

//------------------------------------------------------------------------------

            void
            set_cryo_thermo_flag();

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

            HeatPoly *
            find_heat_poly( const real aT );

//------------------------------------------------------------------------------

            TransportPoly *
            find_viscosity_poly( const real aT );

//------------------------------------------------------------------------------

            TransportPoly *
            find_conductivity_poly( const real aT );

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
            zero( const real aT );

//------------------------------------------------------------------------------

            real
            poly_Cp( const real aT );

//------------------------------------------------------------------------------

            real
            poly_dCpdT( const real aT );

//------------------------------------------------------------------------------

            real
            poly_d2CpdT2( const real aT );

//------------------------------------------------------------------------------

            real
            poly_H( const real aT );

//------------------------------------------------------------------------------

            real
            poly_S( const real aT );

//------------------------------------------------------------------------------

            real
            poly_dSdT( const real aT );

//------------------------------------------------------------------------------

            real
            poly_Mu( const real aT );

//------------------------------------------------------------------------------

            real
            poly_dMudT( const real aT );

//------------------------------------------------------------------------------

            real
            poly_d2MudT2( const real aT );

//------------------------------------------------------------------------------

            real
            poly_Lambda( const real aT );

//------------------------------------------------------------------------------

            real
            poly_dLambdadT( const real aT );

//------------------------------------------------------------------------------

            real
            poly_d2LambdadT2( const real aT );

//------------------------------------------------------------------------------

            real
            spline_Cp( const real aT );

//------------------------------------------------------------------------------

            real
            spline_dCpdT( const real aT );

//------------------------------------------------------------------------------

            real
            spline_H( const real aT );

//------------------------------------------------------------------------------

            real
            spline_S( const real aT );

//------------------------------------------------------------------------------

            real
            spline_dSdT( const real aT );

//------------------------------------------------------------------------------

            real
            spline_Mu( const real aT );

//------------------------------------------------------------------------------

            real
            spline_dMudT( const real aT );

//------------------------------------------------------------------------------

            real
            spline_d2MudT2( const real aT );

//------------------------------------------------------------------------------

            real
            spline_Lambda( const real aT );

//------------------------------------------------------------------------------

            real
            spline_dLambdadT( const real aT );

//------------------------------------------------------------------------------

            real
            spline_d2LambdadT2( const real aT );

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
            return mHaveViscosity;
        }

//------------------------------------------------------------------------------

        bool
        RefGas::has_cryo_thermo() const
        {
            return mHaveCryoThermo;
        }

//------------------------------------------------------------------------------

        bool
        RefGas::has_cryo_transport() const
        {
            return mHaveCryoTransport;
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

//------------------------------------------------------------------------------

        Spline *
        RefGas::viscosity_spline()
        {
            return & mViscositySpline;
        }

//------------------------------------------------------------------------------

        Spline *
        RefGas::conductivity_spline()
        {
            return & mConductivitySpline;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_GT_THERMO_HPP
