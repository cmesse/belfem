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

#ifndef BELFEM_CL_DATABASEPROJECTOR_HPP
#define BELFEM_CL_DATABASEPROJECTOR_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Solver.hpp"
#include "cl_SpMatrix.hpp"

namespace belfem
{
    namespace database
    {
        /**
         * @brief Build-time L2 projection of sampled values onto a B-spline basis.
         *
         * @ingroup grp_physics_database
         * @see @ref physics_database_database_usage_guide
         */
        class Projector
        {
            const proc_t mCommRank ;
            const proc_t mCommSize ;
            const TensorMeshConfig * mConfig ;
            Mesh * mMesh ;
            string mLabel ;
            SpMatrix * mA = nullptr;
            Solver * mSolver = nullptr;
            Matrix< real > mT ;
            Matrix< real > mMel ;
            Matrix< real > mAel ;
            Matrix< real > mBel ;

        public:

            Projector( Mesh * aMesh ) ;

            ~Projector();

            void
            project( const string & aField, Vector< real > & aResult );


            void
            set_label( const string & aLabel );

            const string &
            label() const ;

        private:

            void
            allocate_system_matrix( Mesh * aMesh );

            void
            compute_element_matrices();

        };

        inline void
        Projector::set_label( const string & aLabel )
        {
            mLabel = aLabel ;
        }

        inline const string &
        Projector::label() const
        {
            return mLabel ;
        }

    }
}
#endif // BELFEM_CL_DATABASEPROJECTOR_HPP