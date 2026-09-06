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
#ifndef CL_MAXWELLPOSTPROCESSOR_HPP
#define CL_MAXWELLPOSTPROCESSOR_HPP

#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Postprocessor.hpp"
#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        enum class MaxwellPostprocessorType
        {
            Air,
            Ferro,
            Conductor,      // <- no j/jc
            SuperConductor, // <- has j/jc
            ThinShellConductor,       // <- no j/jc
            ThinShellSuperConductor,  // <- has j/jc
            SideConnector,            // <- edge coating walls: full h = ht + hb + hn
            UNDEFINED
        };

        /**
         * @brief Field recovery for B, H, J and J/Jc.
         *
         * @ingroup grp_fem_maxwell
         * @see @ref fem_maxwell_postprocessor_recovery_theory
         */
        class MaxwellPostprocessor : public Postprocessor
        {
            const MaxwellPostprocessorType mType ;
            const bool mCreateElementFields ;

            Cell< index_t > mMyOwnedElementIndices ;
            Cell< Vector< index_t > > mAllOwnedElementIndices ;

            Vector< real > mX0 ; // coordinates for node
            Vector< real > mX ; // coordinates for integration point
            Vector< real > mY ;
            Vector< real > mZ ;

            Vector< real > mH ;   // magnetic field vector
            Vector< real > mB ;   // magnetic flux density vector
            Vector< real > mJ ;   // current vector
            Vector< real > mJJc ; // current vector, normalized to critical current

            Vector< real > mDOFs ; // dof vector for current element

            uint mRecoveryDepth = 0 ;

            Cell< string > mElementTargetFields ;

            void
            ( MaxwellPostprocessor::*mUpdateFunction )();

            const Vector< real > &
            ( MaxwellPostprocessor::*mComputeFunction )( const uint aK );


            Map< id_t, Material * > mMaterialMap ;
            Matrix< real > mElementData ;

            bool mComputeElementFields = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            MaxwellPostprocessor(
                Kernel * aKernel,
                const Map< id_t, DomainType > & aBlockTypes,
                const Map< id_t, string >     & aMaterialMap,
                const MaxwellPostprocessorType aType,
                const bool aCreateElementFields = false );

            ~MaxwellPostprocessor() override;

            void
            run() override ;

            void
            initialize() override ;

        protected:

            void
            update_element_dofs() override;

        private:

            void
            create_other_fields();

            void
            set_type( const MaxwellPostprocessorType aType );

            void
            create_element_fields();

            void
            select_blocks_and_materials( const Map< id_t, DomainType > & aBlockTypes, const Map< id_t, string >  & aMaterialMap );

            void
            select_owned_elements();

            void
            update_dofs_lagrange();

            void
            update_dofs_nedelec();

            const Vector< real > &
            compute( const uint aK ) override;

            const Vector< real > &
            compute_air( const uint aK );

            const Vector< real > &
            compute_ferro( const uint aK );

            const Vector< real > &
            compute_conductor( const uint aK );

            //! thin-shell layer variant: adds the normal flux recovered
            //! from the volume traces on both sides ( -grad phi on a
            //! phi-region, the conductor's own E*q on an h-conductor,
            //! compute_h_trace ), which the layer's E*q cannot see
            const Vector< real > &
            compute_conductor_ts( const uint aK );

            //! edge coating wall: h = ht + hb + hn via the MaxwellData
            //! recovery, same route as the assembly kernel
            const Vector< real > &
            compute_side_connector( const uint aK );

            //! edge coating wall visualization: copy the tape seam values
            //! of T, H and B onto the wall nodes' own field entries ( the
            //! strip is thin pure metal, so the seam field IS the wall
            //! field to sub-mT accuracy; J stays the wall's own C*q ).
            //! Master only — assembly-side writes are rank-local and the
            //! save-time field gather collects dof-carrying entities only
            void
            copy_seam_fields();

            const Vector< real > &
            compute_superconductor( const uint aK );

            const Vector< real > &
            compute_superconductor_ts( const uint aK );

            void
            compute_element_data();

            void
            collect_element_data();

            DomainType
            select_domain_type( const MaxwellPostprocessorType aType );

        };

        inline void
        MaxwellPostprocessor::update_element_dofs()
        {
            (this->*mUpdateFunction)();
        }


    }
}
#endif //CL_MAXWELLPOSTPROCESSOR_HPP
