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

// The plugin API. A user-defined material, defect or source function needs
// only this header.
//
// NEVER include cl_Vector.hpp or cl_Matrix.hpp here, or any class that uses
// them: the plugin API is backend-free by contract, and every header listed
// below compiles standalone with neither BELFEM_ARMADILLO nor BELFEM_BLAZE
// defined (the three entry headers carry the same note).
//
// gTbulk is only declared here (globals.hpp); a plugin that instantiates the
// power laws takes it as an undefined symbol that the host resolves at load.

#ifndef BELFEM_USER_API_HPP
#define BELFEM_USER_API_HPP
#include "typedefs.hpp"
#include "constants.hpp"
#include "globals.hpp"
#include "fn_sprint.hpp"
#include "assert.hpp"
#include "cl_Bitset.hpp"
#include "cl_Cell.hpp"
#include "cl_Material.hpp"
#include "cl_Material_UserDefined.hpp"
#include "fn_Material_UserDefinedPolynomials.hpp"
#include "cl_JcFunction.hpp"
#include "cl_JcFunction_UserDefined.hpp"
#include "cl_SourceFunction.hpp"
#include "powerlaws.hpp"
#endif //BELFEM_USER_API_HPP
