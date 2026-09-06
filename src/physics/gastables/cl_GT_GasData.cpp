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

#include <cmath>
#include "constants.hpp"
#include "assert.hpp"
#include "cl_GT_GasData.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        void
        GasData::set_label( const string & aLabel )
        {
            mLabel = aLabel;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_name( const string & aName )
        {
            mName = aName;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_cas_number( const string & aCasNumber )
        {
            mCasNumber = aCasNumber;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_molar_mass( const real & aM )
        {
            mM = aM;
            mR = constant::Rm / mM;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_t_crit( const real aTcrit )
        {
            mTcrit = aTcrit;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_p_crit( const real aPcrit )
        {
            mPcrit = aPcrit;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_rho_crit( const real & aRhocrit )
        {
            // the derived value needs p_crit, molar mass and t_crit first
            BELFEM_ASSERT( std::isnan( aRhocrit ) ||
                    ! ( std::isnan( mPcrit ) || std::isnan( mR ) || std::isnan( mTcrit ) ),
                    "set_rho_crit called before p_crit, molar mass and t_crit for %s",
                    mLabel.c_str() );

            mRhocrit = aRhocrit;
            mZcrit = mPcrit / ( mR * mTcrit * mRhocrit );
        }

//------------------------------------------------------------------------------

        void
        GasData::set_z_crit( const real & aZcrit )
        {
            // the derived value needs p_crit, molar mass and t_crit first
            BELFEM_ASSERT( std::isnan( aZcrit ) ||
                    ! ( std::isnan( mPcrit ) || std::isnan( mR ) || std::isnan( mTcrit ) ),
                    "set_z_crit called before p_crit, molar mass and t_crit for %s",
                    mLabel.c_str() );

            mZcrit = aZcrit;
            mRhocrit = mPcrit / ( mR * mTcrit * mZcrit );
        }

//------------------------------------------------------------------------------

        void
        GasData::set_acentric_factor( const real & aOmega )
        {
            mAcentric = aOmega;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_dipole_moment( const real & aDipole )
        {
            mDipole = aDipole;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_formation_enthalpy( const real & aHf )
        {
            mHf = aHf;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_reference_enthalpy( const real & aHref )
        {
            mHref = aHref;
        }

//------------------------------------------------------------------------------

        void
        GasData::set_reference_entropy( const real & aSref )
        {
            mSref = aSref;
        }

//------------------------------------------------------------------------------

    }
}