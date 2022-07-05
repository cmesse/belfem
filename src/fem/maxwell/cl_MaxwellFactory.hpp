//
// Created by christian on 9/20/21.
//

#ifndef BELFEM_CL_MAXWELLFACTORY_HPP
#define BELFEM_CL_MAXWELLFACTORY_HPP

#include "cl_IwgFactory.hpp"
#include "cl_InputFile.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_Kernel.hpp"

#include "cl_IWG_Maxwell.hpp"
#include "cl_FEM_DomainGroup.hpp"
#include "cl_MaxwellMaterial.hpp"
#include "cl_FEM_MaxwellBoundaryConditionCurrent.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"

namespace belfem
{
    class Spline ;

    namespace fem
    {
//------------------------------------------------------------------------------

        struct NonlinearSettings
        {
            //! minumum number of iterations
            const uint minIter;

            //! maximum number of iterations
            const uint maxIter;

            //! relaxation factor
            const real picardOmega;
            const real picardEpsilon;

            const real newtonOmega;
            const real newtonEpsilon;

            NonlinearSettings( const real aMinIter,
                               const real aMaxIter,
                               const real aPicardOmega,
                               const real aPicardEpsilon,
                               const real aNewtonOmega,
                               const real aNewtonEpsilon );

            ~NonlinearSettings() = default;
        };

//------------------------------------------------------------------------------

        class MaxwellFactory : public IwgFactory
        {
            // ref to the input file
            const InputFile & mInputFile  ;

            // parameter object
            KernelParameters * mParameters = nullptr ;

            // kernel object
            Kernel *           mKernel = nullptr ;

            // flag telling if we own the mesh
            // (unset if create_mesh() was called )
            bool mOwnMesh = true ;

            //! flag telling if we own the kernel object
            bool mOwnKernel = true ;

            //! flag telling of we own the boundary conditions
            bool mOwnBoundaryConditions = true ;

            //! flag telling if we want to compute the biot-savart field
            bool mBiotSavartFlag = false ;

            // timestepping info
            real mTimeStep = BELFEM_QUIET_NAN ;
            real mTheta    = BELFEM_QUIET_NAN ;

            // list of selected block ids
            Vector< id_t > mBlockIDs ;
            Vector< id_t > mSideSetIDs ;
            Cell< DomainType >  mBlockTypes ;
            Cell< DomainType >  mSideSetTypes ;

            // list of nedelec blocks
            Vector< id_t > mNedelecBlocks ;
            Vector< id_t > mNedelecSideSets ;

            // written by taperoller
            Vector< id_t > mGhostSideSets ;

            // this one is needed to tell the IWG where the DOFs are picket from
            id_t mGhostMaster = 0 ;

            // domain groups
            Cell< DomainGroup * >        mBlocks ;
            Map< string, DomainGroup * > mBlockMap ;

            Cell< DomainGroup * >        mBoundaries ;
            Map< string, DomainGroup * > mBoundaryMap ;

            Cell< DomainGroup * >        mInterfaces ;

            Cell< DomainGroup * >        mCuts ;
            Map< string, DomainGroup * > mCutMap ;

            Cell< DomainGroup * >        mSideSets ;

            Cell< DomainGroup * >        mTapes ;

            Cell< Material * >        mMaterials ;
            Map< string, Material * > mMaterialMap ; // links materials to block names
            Map< string, Material * > mMaterialLabelMap ;  // links materials to material labels

            // label for formulation
            IwgType mFormulation = IwgType::UNDEFINED ;
            IwgType mProjectJ    = IwgType::UNDEFINED ;
            IwgType mProjectB    = IwgType::UNDEFINED ;

            Cell< BoundaryCondition * > mBoundaryConditions ;
            Map< string, MaxwellBoundaryConditionCurrent * > mCurrentBcMap ;

            Map< string, MaxwellBoundaryConditionMagfield * > mMagfieldBcMap ;

            //! list of blocks for which we compute the element-wise current
            Vector< id_t > mCurrentBlocks ;

            //! flag telling if linear field is enforced for maxwell
            bool mEnforceLinear = false ;

            //! flag telling if J-Field is computed
            bool mHaveJ = false ;

            //! flag telling if B-Field is computed
            bool mHaveB = false ;

            //! projection factor for j
            real mAlphaJ = 0.0 ;

            //! projector factor for b
            real mAlphaB = 0.0 ;

            //! flag telling if thermal field exists
            bool mHaveThermal = false ;

            //! flag telling if mech field exists
            bool mHaveMech = false ;

            //! flag telling if we want to compute biot-savart
            bool mComputeBiotSavart = false ;

            //! flag if we want to compute norm b
            bool mComputeNormB = false ;

            DofManager * mMagneticField = nullptr ;
            DofManager * mThermalField  = nullptr ;

            Map< id_t, id_t > mSideSetToCutMap ;

            //! for tapes
            Cell< string > mTapeMaterialLabels ;
            Vector< real > mTapeThicknesses ;

            //! container with tape materials. Todo: add to sidesets
            Cell< Material * > mTapeMaterials ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            MaxwellFactory( const InputFile & aInputFile ) ;

            ~MaxwellFactory();

//------------------------------------------------------------------------------

            Mesh *
            mesh() ;

            KernelParameters *
            parameters() ;

            Kernel *
            kernel() ;

            DofManager *
            magnetic_field() ;

            /**
             * return the name of the exodus file
             */
            string
            outfile() const ;

            /**
             * return the name of the backup file
             */
            string
            backupfile() const ;

            bool
            restart() const ;

            real
            maxtime() const ;

            NonlinearSettings
            nonlinear_settings() ;

//------------------------------------------------------------------------------

            bool
            compute_biot_savart() const;

//------------------------------------------------------------------------------

            bool
            compute_normb() const;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * read the formulation which determines the IWG type
             * form the input file
             */
             void
             read_formulation();

//------------------------------------------------------------------------------

            /**
             * read the relevant information for the timestepping method
             */
            void
            read_timestep();

//------------------------------------------------------------------------------

            /**
             * secondary fields such as
             * j, b, h
             * are computed in a postprocessing step
             */
            void
            read_field_switches();

//-----------------------------------------------------------------------------

            void
            remove_coils();

//-----------------------------------------------------------------------------

            void
            set_sheet_types();

//------------------------------------------------------------------------------

            void
            create_bcs();

//------------------------------------------------------------------------------

            void
            create_materials();

//------------------------------------------------------------------------------

            Material *
            create_maxwell_material( const input::Section * aInput ) ;

            void
            create_layers();

//------------------------------------------------------------------------------

            void
            link_materials( DofManager * aDofManager );

//------------------------------------------------------------------------------

            /**
             * a subroutine for the create_domain_groups() routine.
             * Builds a new group based on a section in the input file.
             */
            DomainGroup *
            create_domain_group( const input::Section * aSection );

//------------------------------------------------------------------------------

            // called by constructor
            Mesh *
            load_mesh( const InputFile & aInputFile  );

//------------------------------------------------------------------------------

            IWG_Maxwell *
            create_iwg( const IwgType aType );

//------------------------------------------------------------------------------

            void
            set_tape_materials( DofManager * aDofManager );

//------------------------------------------------------------------------------

            void
            create_kernel();

//------------------------------------------------------------------------------

            void
            create_topology();

//------------------------------------------------------------------------------

            void
            create_magnetic_field();

// -----------------------------------------------------------------------------

            void
            set_solver( DofManager * aDofMaganer, const input::Section * aSection );

// -----------------------------------------------------------------------------

            void
            create_current_projector();

// -----------------------------------------------------------------------------

            void
            create_magfield_projector();

// -----------------------------------------------------------------------------

            void
            add_postprocessor_to_kernel( IWG_Maxwell * aEquation );

// -----------------------------------------------------------------------------

            // if a sideset has a dirichlet bc, we must remove it from the IWG
            void
            blacklist_dirichlet_sidesets();

// -----------------------------------------------------------------------------
// for topology
// -----------------------------------------------------------------------------

            void
            create_domain_groups();

// -----------------------------------------------------------------------------

            void
            select_blocks();

// -----------------------------------------------------------------------------

            // populate physical tags in mesh
            void
            write_block_phystags_to_mesh();

// -----------------------------------------------------------------------------

            void
            write_sideset_phystags_to_mesh();

// -----------------------------------------------------------------------------

            void
            fix_facet_masters();

// -----------------------------------------------------------------------------

            void
            create_interfaces();

// -----------------------------------------------------------------------------

            void
            select_bcs_and_cuts();

// -----------------------------------------------------------------------------

            void
            select_sidesets();

// -----------------------------------------------------------------------------

            Spline *
            read_bhfile(
                    const string & aPath,
                    const string & aUnitB,
                    const string & aUnitH,
                    const value    aMaxB,
                          real   & aM );

// -----------------------------------------------------------------------------

            /**
             * help function to determine the element type of a block
             */
             ElementType
             element_type( const id_t aBlockID );

// -----------------------------------------------------------------------------

            void
            set_domain_types() ;

// -----------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

        inline Mesh *
        MaxwellFactory::mesh()
        {
            mOwnMesh = false ;
            return mMesh ;
        }

//------------------------------------------------------------------------------

        inline KernelParameters *
        MaxwellFactory::parameters()
        {
            // if this is called, it's the user's responsibility to destroy
            mKernel->claim_parameter_ownership( false );
            return mParameters ;
        }

//------------------------------------------------------------------------------

        inline Kernel *
        MaxwellFactory::kernel()
        {
            mOwnKernel = false ;
            return mKernel ;
        }

//------------------------------------------------------------------------------

        inline DofManager *
        MaxwellFactory::magnetic_field()
        {
            return mMagneticField ;
        }

//------------------------------------------------------------------------------

        inline bool
        MaxwellFactory::compute_biot_savart() const
        {
            return mComputeBiotSavart ;
        }

//------------------------------------------------------------------------------

        inline bool
        MaxwellFactory::compute_normb() const
        {
            return mComputeNormB && mHaveB ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MAXWELLFACTORY_HPP