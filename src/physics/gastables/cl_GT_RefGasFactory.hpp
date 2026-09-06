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

#ifndef BELFEM_CL_GT_REFGASFACTORY_HPP
#define BELFEM_CL_GT_REFGASFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_GT_RefGas.hpp"
#include "fn_GT_data_path.hpp"

namespace belfem
{
    class SpMatrix;

    namespace gastables
    {
        class InputThermo;
        class InputTransport;
        class InputData;
        class InputAlpha;

//------------------------------------------------------------------------------

        /**
         * Builds RefGas objects from the shipped data tables.
         *
         * The factory owns the four readers for the lifetime of the object, so
         * the tables are parsed once however many species are created. Its work
         * is more than lookup: a raw record covers only the temperature range
         * the source data covers, and create_refgas fills the rest in - glue
         * polynomials across interval junctions, a cryogenic interval below the
         * lowest tabulated one, a hot extrapolation above the highest, and, for
         * species that have a critical point but no transport record,
         * viscosity and conductivity synthesized from the Lucas and Chung
         * correlations.
         *
         * The returned RefGas is owned by the caller and must be deleted by it.
         * The factory holds no reference to what it produced and may be
         * destroyed first.
         *
         * The data path defaults to gastables::data_path(), which resolves
         * $BELFEM_DATA and otherwise searches for share/fluid relative to
         * the working directory.
         *
         * @ingroup grp_physics_gastables
         * @see @ref physics_gastables_gastables_usage_guide
         */
        class RefGasFactory
        {
            InputThermo    * mThermo;
            InputTransport * mTransport;
            InputData      * mData;
            InputAlpha     * mAlpha;

            Vector< real >   mTemperatures;
            SpMatrix         mHelpMatrix;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            RefGasFactory( const string aDataPath = gastables::data_path() );

//------------------------------------------------------------------------------

            ~RefGasFactory();

//------------------------------------------------------------------------------

            // owns the four input parsers
            RefGasFactory( const RefGasFactory & ) = delete;
            RefGasFactory & operator=( const RefGasFactory & ) = delete;

//------------------------------------------------------------------------------

            RefGas *
            create_refgas( const string & aLabel );

//------------------------------------------------------------------------------

            void
            create_temperature_steps( Vector< real > & aTemperatureSteps );

//------------------------------------------------------------------------------

            void
            create_helpmatrix( SpMatrix & aHelpMatrix );

//------------------------------------------------------------------------------

            bool
            interaction_viscosity_exists(
                    const string & aA,
                    const string & aB );

//------------------------------------------------------------------------------

            RefGas *
            create_interaction_viscosity( const string & aA, const string & aB );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_REFGASFACTORY_HPP
