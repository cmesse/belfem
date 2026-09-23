/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * A linear solver section read from a deck goes through the same setters a
 * C++ caller uses, so a value the setter refuses is refused from a deck too.
 * The sparse suite has no input parser, which is why these cases live here.
 */

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_InputFile.hpp"
#include "cl_SolverParameters.hpp"

namespace
{
    //! a deck with one linear solver section in the test temp directory.
    //! The path is absolute, because a relative one is resolved against a
    //! directory ctest may have changed; the file is removed again on scope exit
    struct LinearDeck
    {
        std::string mPath ;

        LinearDeck( const std::string & aName, const std::string & aBody ) :
            mPath( ::testing::TempDir() + aName )
        {
            std::ofstream tFile( mPath );
            tFile << "solver\n{\n    linear\n    {\n" << aBody << "\n    }\n}\n";
        }

        ~LinearDeck()
        {
            std::remove( mPath.c_str() );
        }
    };

    //! builds the parameters of the deck's linear section
    belfem::SolverParameters
    read_linear_section( const LinearDeck & aDeck )
    {
        belfem::InputFile tInput( aDeck.mPath );
        return belfem::SolverParameters( tInput.section( "solver" )->section( "linear" ) );
    }
}

TEST( SolverParametersDeck, RelativeToleranceIsValidatedFromADeck )
{
    // Mutation: let the deck arm write the member directly again, bypassing
    // the setter, and the refusals below go red
    for ( const char * tValue : { "0", "-1e-6" } )
    {
        const LinearDeck tDeck( "linear_rtol_bad.conf",
                std::string( "        relative tolerance : " ) + tValue + " ;" );
        EXPECT_THROW( read_linear_section( tDeck ), std::runtime_error )
                << "relative tolerance : " << tValue ;
    }

    const LinearDeck tDeck( "linear_rtol_good.conf", "        relative tolerance : 1e-9 ;" );
    const belfem::SolverParameters tParams = read_linear_section( tDeck );
    EXPECT_DOUBLE_EQ( tParams.relative_tolerance(), 1e-9 );
    EXPECT_TRUE( tParams.have_relative_tolerance() );
}

TEST( SolverParametersDeck, MaxIterationsIsValidatedFromADeck )
{
    // a negative count must be refused, not wrapped into a huge budget
    for ( const char * tValue : { "0", "-1" } )
    {
        const LinearDeck tDeck( "linear_maxit_bad.conf",
                std::string( "        max iterations : " ) + tValue + " ;" );
        EXPECT_THROW( read_linear_section( tDeck ), std::runtime_error )
                << "max iterations : " << tValue ;
    }

    const LinearDeck tDeck( "linear_maxit_good.conf", "        max iterations : 7 ;" );
    EXPECT_EQ( read_linear_section( tDeck ).max_iterations(), 7u );
}
