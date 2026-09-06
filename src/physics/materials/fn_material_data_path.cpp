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

#include "filetools.hpp"
#include "globals.hpp"
#include "fn_material_data_path.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        //! subdirectory of the data path that holds the material databases.
        //! Bare, without a separator: both readers below join explicitly
        constexpr char gMaterialsSubdir[] = "material";

//------------------------------------------------------------------------------

        string
        data_path()
        {
            return gBelfemDataPath.size() > 0 ?
                gBelfemDataPath + "/" + gMaterialsSubdir : "";
        }

//------------------------------------------------------------------------------

        string
        data_file( const string & aFile )
        {
            // the search order itself lives in io/filetools, one layer down,
            // so that the source-function plugin loader in numerics/sources
            // can reach it too -- numerics sits BELOW physics and cannot
            // include this header
            return search_data_file( aFile, gMaterialsSubdir );
        }

//------------------------------------------------------------------------------
    }
}
