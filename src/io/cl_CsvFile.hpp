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

#ifndef BELFEM_CL_CSVFILE_HPP
#define BELFEM_CL_CSVFILE_HPP
#include "cl_Ascii.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    /**
     * @brief CSV reader for numeric data.
     *
     * @ingroup grp_io
     * @see @ref io_io_usage_guide
     */
    class CsvFile : public Ascii
    {
        Matrix< real > mData ;

        public:

        CsvFile( const string & aPath, const char aDelimiter=',');

        const Matrix< real > & data() const { return mData; }

        Matrix< real > & data() { return mData; }

        ~CsvFile() override;
    };
}
#endif //BELFEM_CL_CSVFILE_HPP
