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
 * The block domain types the thermal factory admits. One predicate serves
 * the block selection and the dispatch over the selected blocks, so the
 * dispatch's refusal cannot fire on a selected block; this pins the
 * predicate itself over every enumerator, so a type added to the enum
 * lands in a decision rather than in the refused set by accident.
 */

#include <gtest/gtest.h>
#include <set>

#include "typedefs.hpp"
#include "en_DomainType.hpp"
#include "cl_ThermalFactory.hpp"

using namespace belfem ;

TEST( ThermalBlockAdmission, ExactlyTheFourConductingTypesAreSupported )
{
    const std::set< DomainType > tSupported = {
        DomainType::Conductor,
        DomainType::ThinShell,
        DomainType::Ferro,
        DomainType::Buffer };

    const uint tLast = static_cast< uint >( DomainType::UNDEFINED );

    for ( uint k = 0; k <= tLast; ++k )
    {
        const DomainType tType = static_cast< DomainType >( k );
        EXPECT_EQ( fem::ThermalFactory::block_is_supported( tType ),
                   tSupported.count( tType ) == 1 )
            << "domain type " << k << " ( " << to_string( tType ) << " )" ;
    }
}

TEST( ThermalBlockAdmission, AirIsNotAThermalBlock )
{
    EXPECT_FALSE( fem::ThermalFactory::block_is_supported( DomainType::Air ) );
    EXPECT_FALSE( fem::ThermalFactory::block_is_supported( DomainType::Coil ) );
    EXPECT_FALSE( fem::ThermalFactory::block_is_supported( DomainType::Default ) );
}
