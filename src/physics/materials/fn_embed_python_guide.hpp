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

#ifndef BELFEM_FN_EMBED_PYTHON_GUIDE_HPP
#define BELFEM_FN_EMBED_PYTHON_GUIDE_HPP

#include "typedefs.hpp"

namespace belfem
{
    class HDF5 ;

    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Copy the python database reader into a generated table
         *
         * If ``$BELFEM_DATA/python/database/howtoread.py`` exists, its full
         * text is stored in the open database file as group ``python``,
         * dataset ``howtoread.py``. A table written by the material tool then
         * ships its own reference reader: anyone holding the file can recover
         * the properties without a BELFEM installation.
         *
         * The feature is optional by design -- when the guide is not found
         * the function is a silent no-op, so table generation never depends
         * on the data directory being complete.
         *
         * @param aFile   open HDF5 database, positioned at the file root
         * @return        true if the guide was found and embedded
         *
         * The caller is responsible for the rank guard: like the save
         * routines that call it, this runs on rank 0 only.
         */
        bool
        embed_python_guide( HDF5 & aFile );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_EMBED_PYTHON_GUIDE_HPP
