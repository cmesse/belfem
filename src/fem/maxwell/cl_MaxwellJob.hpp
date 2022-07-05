//
// Created by Christian Messe on 13.01.22.
//

#ifndef BELFEM_CL_MAXWELLJOB_HPP
#define BELFEM_CL_MAXWELLJOB_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_MaxwellMaterial.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * The job class is a meta class that contains the data needed to run a job.
         * It is also used to run the tests.
         */
        class MaxwellJob
        {
            //! the rank of this proc
            const proc_t mMyRank ;

            //! number of procs
            const proc_t mCommSize ;

            //! the mesh. This is the real mesh object for the master,
            //! otherwise it points to the mesh on the kernel
            Mesh * mMesh = nullptr ;

            //! the main equation object (owned and deleted by the Kernel )
            fem::IWG_Maxwell * mEquation = nullptr ;

            //! the kernel parameters for this job
            KernelParameters * mParams = nullptr ;

            //! the kernel itself
            Kernel * mKernel = nullptr ;

            DofManager * mMagneticField = nullptr ;

            string mLabel = "job";

            Cell< Material * > mMaterials ;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * this constructor only sets the comm size etc.
             * Needed for Tests.
             */
             MaxwellJob();

//------------------------------------------------------------------------------

             /**
              * this constructor will be called by the Maxwell factory
              */
             MaxwellJob( Mesh * aMesh, Kernel * aKernel );

//------------------------------------------------------------------------------

             ~MaxwellJob();

//------------------------------------------------------------------------------
            /**
             * this subroutine sets a fully functional reference problem
             * for a non-physical example. The goal is to compute the element
             * matrices and make sure that they are consistent with known solutions.
             */
             void
             initialize_test( const IwgType aIwgType, const bool aCurved = false );

//------------------------------------------------------------------------------

             void
             run_test( const string & aReferenceFilePath, Vector< real > & aResult, const bool aPrint = false );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_test_mesh_tri3() ;

//------------------------------------------------------------------------------

            void
            create_test_mesh_tri6( const bool aCurved=true ) ;

//------------------------------------------------------------------------------

            void
            create_test_mesh_tet4() ;

//------------------------------------------------------------------------------

            void
            create_test_mesh_tet10( const bool aCurved=true ) ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_MAXWELLJOB_HPP
