/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * jc and n have no T-only custom callback.
 *
 * Material::set_custom( jc ) and set_custom( n ) once wrote their callbacks
 * into the rho_i and debye dispatch slots, while the power law evaluates jc
 * and n through a JcFunction and never read them. Both cases now abort, so
 * the resistivity and Debye routing of a material cannot be hijacked.
 */

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_Material_Copper.hpp"

TEST( JcCustomGuard, SetCustomRejectsJcAndN )
{
    // a built-in metal: tables off, analytic Kohler path
    belfem::material::Copper tCopper( 50.0, false );

    EXPECT_THROW( tCopper.set_custom( belfem::MaterialProperty::jc ), std::exception );
    EXPECT_THROW( tCopper.set_custom( belfem::MaterialProperty::n ),  std::exception );

    // the guard fires before the have-flag or any dispatch slot is touched
    EXPECT_FALSE( tCopper.have( belfem::MaterialProperty::jc ) );
    EXPECT_FALSE( tCopper.have( belfem::MaterialProperty::n ) );
    EXPECT_GT( tCopper.rho_i( 77.0 ), 0.0 );
    EXPECT_GT( tCopper.debye( 77.0 ), 0.0 );
}
