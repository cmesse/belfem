/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_VTK_CURVE_HPP
#define BELFEM_CL_VTK_CURVE_HPP

#include "vtktypes.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace vtk
    {
        /**
         * @class Curve
         * @brief A polyline VTK actor built from an ordered list of 3D points.
         *
         * Generic and source-agnostic: it knows nothing about meshes, orbits, or
         * kepler — it just turns points into a polyline. Points are a 3 x N
         * matrix (each column one point, contiguous in BELFEM's column-major
         * storage). Pass aClosed = true to join the last point back to the first
         * (e.g. a closed orbit ellipse). Rendered unlit so the curve reads as a
         * flat overlay rather than a shaded tube.
         *
         * @ingroup grp_visualizer
         * @see @ref visualizer_index
         */
        class Curve
        {
            Points    mPoints ;
            PolyData  mPolyData ;
            Mapper    mMapper ;
            Actor     mActor ;

        public:

            //! Build a polyline actor from a 3 x N point matrix.
            Curve( const Matrix< real > & aPoints,
                   real aLineWidth = 2.0,
                   bool aClosed    = false );

            ~Curve() = default;

            //! The VTK actor; ownership stays with this wrapper.
            inline Actor actor() { return mActor; }

            //! Set the line color (RGB components in [0,1]).
            void set_color( real aRed, real aGreen, real aBlue );
        };
    }
}
#endif //BELFEM_CL_VTK_CURVE_HPP
