/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * The symmetry default of the IWG base constructor. A subclass that omits
 * the argument must come out Unsymmetric: MUMPS is handed the full matrix,
 * and a symmetric mode makes it factorize a different one ( cl_IWG.hpp,
 * the comment on the constructor ). Free of Maxwell headers on purpose:
 * this translation unit is compiled unconditionally.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_IWG.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    // the shape the user guide documents for IWG_Custom: type and
    // dimensionality only, every other argument left to its default
    class IwgProbe : public fem::IWG
    {
        public:
            IwgProbe() :
                fem::IWG( IwgType::UNDEFINED, ModelDimensionality::TwoD )
            {
            }
    };

    // the same subclass with the symmetry stated
    class SymmetricIwgProbe : public fem::IWG
    {
        public:
            SymmetricIwgProbe() :
                fem::IWG( IwgType::UNDEFINED, ModelDimensionality::TwoD,
                          IwgMode::Iterative, SymmetryMode::PositiveDefiniteSymmetric )
            {
            }
    };
}

TEST( IwgDefaults, SubclassOmittingTheArgumentIsUnsymmetric )
{
    IwgProbe tProbe ;
    EXPECT_EQ( tProbe.symmetry_mode(), SymmetryMode::Unsymmetric );
}

TEST( IwgDefaults, ExplicitSymmetricIsHonored )
{
    // the accessor reports what the constructor was given, so the case
    // above pins the default and not a constant accessor
    SymmetricIwgProbe tProbe ;
    EXPECT_EQ( tProbe.symmetry_mode(), SymmetryMode::PositiveDefiniteSymmetric );
}
