//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_GASDATA_HPP
#define BELFEM_CL_GT_GASDATA_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {
        class Helmholtz;
    }

    namespace gastables
    {

        class RefGas;
        class InputData;
        class InputAlpha;

//------------------------------------------------------------------------------

        class GasData
        {
            // the name of the gas as chemical symbol
            string mLabel = "uintitled";

            // the name of the gas as written word
            string mName = "untitled";

            // the cas number of this gas
            string mCasNumber = "none";

            // specific gas constant in J/(kg*K)
            real mR = BELFEM_QUIET_NAN;

            // the molar mass in kg/mol;
            real mM = BELFEM_QUIET_NAN;

            // critical temperature in K
            real mTcrit = BELFEM_QUIET_NAN;

            // critical pressure in Pascal
            real mPcrit= BELFEM_QUIET_NAN;

            // critical density in kg / m^3
            real mRhocrit= BELFEM_QUIET_NAN;

            // real gas factor
            real mZcrit= BELFEM_QUIET_NAN;

            // acentric factor
            real mAcentric = BELFEM_QUIET_NAN;

            // dipole moment
            real mDipole = BELFEM_QUIET_NAN;

            // coeffs for SRK EoS
            Vector< real > mCoeffsSRK;

            // coeffs for PR EOS
            Vector< real > mCoeffsPR;

            // flag telling if data is taken from Mahmoodi and Sedig
            // 10.1016/j.fluid.2016.12.015
            bool mCubicFlag = false;

            // cell with elements this gas
            Cell< string > mElements;

            // map containing the composition
            Map< string, real > mComposition;

            // formation enthalpy at 298.15 K
            real mHf = BELFEM_QUIET_NAN;

            // enthalpy at 298.15 K
            real mHref = BELFEM_QUIET_NAN;

            // entropy at 298.15 K and 1 bar ( calculated from spline )
            real mSref = BELFEM_QUIET_NAN;

            friend RefGas;
            friend InputData;
            friend InputAlpha;
            friend gasmodels::Helmholtz ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            GasData() = default;

            ~GasData() = default;

//------------------------------------------------------------------------------

            inline bool
            has_crit() const;

//------------------------------------------------------------------------------

            inline bool
            has_cubic() const;

//------------------------------------------------------------------------------

            inline const string &
            label() const;

//------------------------------------------------------------------------------

            inline const string &
            name() const;
//------------------------------------------------------------------------------

            inline const string &
            cas() const;

//------------------------------------------------------------------------------

            inline const real &
            R() const;

//------------------------------------------------------------------------------

            inline const real &
            M() const;

//------------------------------------------------------------------------------

            inline const real &
            acentric() const;

//------------------------------------------------------------------------------

            inline const real &
            dipole() const;

//------------------------------------------------------------------------------

            inline const real &
            Hf() const;

//------------------------------------------------------------------------------

            inline const real &
            Href() const;

//------------------------------------------------------------------------------

            inline const real &
            Sref() const;

//------------------------------------------------------------------------------

            // cell with elements this gas
            inline const Cell< string > &
            elements() const;

//------------------------------------------------------------------------------

            // map containing the composition
            inline const Map< string, real > &
            composition() const;

//------------------------------------------------------------------------------

            inline real
            component_multiplicity( const string & aLabel ) const;

//------------------------------------------------------------------------------

            inline const real &
            T_crit() const;

//------------------------------------------------------------------------------

            inline const real &
            p_crit() const;

//------------------------------------------------------------------------------

            inline const real &
            rho_crit() const;

//------------------------------------------------------------------------------

            inline const real &
            Z_crit() const;


//------------------------------------------------------------------------------

            inline const Vector< real > &
            srk() const;

//------------------------------------------------------------------------------

            inline const Vector< real > &
            pr() const;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            set_label( const string & aLabel );

//------------------------------------------------------------------------------

            void
            set_name( const string & aName );

//------------------------------------------------------------------------------

            void
            set_cas_number( const string & aCasNumber );

//------------------------------------------------------------------------------

            void
            set_molar_mass( const real & aM );

//------------------------------------------------------------------------------

            void
            set_t_crit( const real aTcrit );

//------------------------------------------------------------------------------

            void
            set_p_crit( const real aPcrit );

//------------------------------------------------------------------------------

            void
            set_rho_crit( const real & aRhocrit );

//------------------------------------------------------------------------------

            void
            set_z_crit( const real & aZcrit );

//------------------------------------------------------------------------------

            void
            set_acentric_factor( const real & aOmega );

//------------------------------------------------------------------------------

            void
            set_dipole_moment( const real & aDipole );

//------------------------------------------------------------------------------

            void
            set_formation_enthalpy( const real & aHf );

//------------------------------------------------------------------------------

            void
            set_reference_enthalpy( const real & aHref );

//------------------------------------------------------------------------------

            void
            set_reference_entropy( const real & aSref );

//------------------------------------------------------------------------------

            inline Vector< real > &
            srk();

//------------------------------------------------------------------------------

            inline Vector< real > &
            pr();

//------------------------------------------------------------------------------

            inline void
            set_cubic_flag();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        bool
        GasData::has_crit() const
        {

            return         (     ! std::isnan( mTcrit    ) )
                           &&  ( ! std::isnan( mPcrit    ) )
                           &&  ( ! std::isnan( mZcrit    ) )
                           &&  ( ! std::isnan( mAcentric ) )
                           &&  ( ! std::isnan( mDipole   ) );
        }

//------------------------------------------------------------------------------

        bool
        GasData::has_cubic() const
        {
            return mCubicFlag;
        }

//------------------------------------------------------------------------------

        const string &
        GasData::label() const
        {
            return mLabel;
        }

//------------------------------------------------------------------------------

        const string &
        GasData::name() const
        {
            return mName;
        }

//------------------------------------------------------------------------------

        const string &
        GasData::cas() const
        {
            return mCasNumber;
        }

//------------------------------------------------------------------------------

        const real &
        GasData::R() const
        {
            return mR;
        }

//------------------------------------------------------------------------------

        const real &
        GasData::M() const
        {
            return mM;
        }

//------------------------------------------------------------------------------

        inline const real &
        GasData::acentric() const
        {
            return mAcentric;
        }

//------------------------------------------------------------------------------

        inline const real &
        GasData::dipole() const
        {
            return mDipole;
        }

//------------------------------------------------------------------------------

        inline const real &
        GasData::Hf() const
        {
            return mHf;
        }

//------------------------------------------------------------------------------

        inline const real &
        GasData::Href() const
        {
            return mHref;
        }

//------------------------------------------------------------------------------

        inline const real &
        GasData::Sref() const
        {
            return mSref;
        }

//------------------------------------------------------------------------------

        // cell with elements this gas
        const Cell< string > &
        GasData::elements() const
        {
            return mElements;
        }

//------------------------------------------------------------------------------

        // map containing the composition
        const Map< string, real > &
        GasData::composition() const
        {
            return mComposition;
        }

//------------------------------------------------------------------------------

        real
        GasData::component_multiplicity( const string & aLabel ) const
        {
            // test if parameter exists
            if ( mComposition.key_exists( aLabel ) )
            {
                return mComposition( aLabel );
            }
            else
            {
                return  0.0;
            }
        }

//------------------------------------------------------------------------------

        const real &
        GasData::T_crit() const
        {
            return mTcrit;
        }

//------------------------------------------------------------------------------

        const real &
        GasData::p_crit() const
        {
            return mPcrit;
        }

//------------------------------------------------------------------------------

        const real &
        GasData::rho_crit() const
        {
            return mRhocrit;
        }

//------------------------------------------------------------------------------

        const real &
        GasData::Z_crit() const
        {
            return mZcrit;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        GasData::srk()
        {
            return mCoeffsSRK;
        }

        const Vector< real > &
        GasData::srk() const
        {
            return mCoeffsSRK;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        GasData::pr()
        {
            return mCoeffsPR;
        }

        const Vector< real > &
        GasData::pr() const
        {
            return mCoeffsPR;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_cubic_flag()
        {
            mCubicFlag = true;
        }

//------------------------------------------------------------------------------

    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_GASDATA_HPP
