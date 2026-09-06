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

#include "cl_GT_RefGasFactory.hpp"


#include "cl_Spline.hpp"
#include "cl_SpMatrix.hpp"
#include "fn_linspace.hpp"

#include "GT_globals.hpp"
#include "fn_GT_data_path.hpp"
#include "cl_GT_InputThermo.hpp"
#include "cl_GT_InputTransport.hpp"
#include "cl_GT_InputData.hpp"
#include "cl_GT_InputAlpha.hpp"
#include "cl_GT_RefGas.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        RefGasFactory::RefGasFactory( const string aDataPath )
        {
            BELFEM_ERROR( aDataPath.size() > 0,
                    "No gas data directory found. Set BELFEM_DATA or run from a "
                    "directory with share/fluid in reach." );

            // open the thermo file
            mThermo    = new InputThermo( aDataPath + "/thermo.inp" );

            // open the transport file
            mTransport = new InputTransport( aDataPath + "/trans.inp" );

            // open the data file
            mData     = new InputData( aDataPath + "/gasdata.inp" );

            // open the cubic file
            mAlpha    = new InputAlpha( aDataPath + "/cubicalpha.inp" );

            // temperatures used for the splines
            this->create_temperature_steps( mTemperatures );

            // help matrix for spline creation ( actully only needed on master proc )
            this->create_helpmatrix( mHelpMatrix );
        }

//------------------------------------------------------------------------------

        RefGasFactory::~RefGasFactory()
        {
            delete mThermo;
            delete mTransport;
            delete mData;
            delete mAlpha;
        }

//------------------------------------------------------------------------------

        void
        RefGasFactory::create_helpmatrix( SpMatrix & aHelpMatrix )
        {
            spline::create_helpmatrix(
                    gastables::gNumberOfSplinePoints,
                    gastables::gDeltaT,
                    aHelpMatrix );
        }

//------------------------------------------------------------------------------

        void
        RefGasFactory::create_temperature_steps( Vector< real > & aTemperatureSteps )
        {
            linspace(
                    0.0,
                    gastables::gTmax,
                    gastables::gNumberOfSplinePoints,
                    aTemperatureSteps );
        }

//------------------------------------------------------------------------------

        RefGas *
        RefGasFactory::create_refgas( const string & aLabel )
        {
            // create the refgas pointer
            RefGas * aRefGas = new RefGas( aLabel );

            // an unknown label would otherwise produce a silent zero/NaN gas
            BELFEM_ERROR( mThermo->entry_exists( aRefGas ),
                    "Species %s not found in thermo.inp ( labels are case sensitive )",
                    aLabel.c_str() );

            // load thermal data. Low temperature intervals are part of the record
            // in thermo.inp, so there is no second file and no order to respect.
            mThermo->read_data( aRefGas );

            // load transport data
            mTransport->read_data( aRefGas );

            // load cas number, masses and critical data
            mData->read_data( aRefGas->data() );

            // load srk and pr data
            mAlpha->read_data( aRefGas->data() );

            // carefully create interpolating and extrapolating functions
            aRefGas->finalize();

            // create interpolating splines
            aRefGas->create_splines( mTemperatures, mHelpMatrix );

            // return the object
            return aRefGas;
        }

//------------------------------------------------------------------------------

        bool
        RefGasFactory::interaction_viscosity_exists(
                const string & aA,
                const string & aB )
        {
            return mTransport->interaction_parameter_exists( aA, aB );
        }

//------------------------------------------------------------------------------

        RefGas *
        RefGasFactory::create_interaction_viscosity( const string & aA, const string & aB )
        {
            // make sure that this call is legal; must stay active in release,
            // otherwise a missing database entry produces a garbage gas
            BELFEM_ERROR( mTransport->interaction_parameter_exists( aA, aB ),
                    "Interaction parameter of %s and %s does not exist.",
                          aA.c_str(), aB.c_str() );

            // create the refgas pointer
            auto aRefGas = new RefGas( aA + "@" + aB );

            // read data from transport only
            mTransport->read_data( aRefGas );

            // finalize the gas
            aRefGas->finalize();

            // create interpolating splines
            aRefGas->create_splines( mTemperatures, mHelpMatrix );

            // return pointer
            return aRefGas;
        }

    }
}