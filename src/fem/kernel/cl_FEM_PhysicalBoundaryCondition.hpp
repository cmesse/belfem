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

#ifndef BELFEM_CL_FEM_BOUNDARYCONDITION_HPP
#define BELFEM_CL_FEM_BOUNDARYCONDITION_HPP

#include "cl_InputFile.hpp"
#include "cl_SourceFunction.hpp"
#include "en_FEM_BoundaryConditionType.hpp"

namespace belfem
{
    namespace fem
    {
//-----------------------------------------------------------------------------

        // pointer-only members below: the .cpp includes
        // cl_FEM_DofManager.hpp ( which itself pulls cl_FEM_Dof.hpp ).
        // Keeping the fat DofManager header out of this one keeps the
        // circuit library ( which includes this file ) off the
        // DofManager → Mesh/IWG/Material include graph — the fact that
        // lets circuit's CMake drop those include dirs
        class DofManager ;
        class Dof ;

//-----------------------------------------------------------------------------

        class PhysicalBoundaryCondition
        {
//-----------------------------------------------------------------------------
        protected :
//-----------------------------------------------------------------------------

            const index_t mIndex ;

            BoundaryConditionType mType ;

            //! Exodus/global-variable name of this condition, assigned by the
            //! creating factory from the input-file group ( section header
            //! label, else type + occurrence ordinal ). Empty = this BC does
            //! not publish a global ( e.g. Bearing, which imposes no scalar )
            string mLabel = "" ;

            SourceFunction * mFunction = nullptr ;

            real mValue = BELFEM_EPS ;

            //Scaling
            real mScale = 1.0 ;

            Vector < real > mDirection = Vector < real >(3,0.0) ;

            unit mUnits;

            //! the deck stated the amplitude as a flux density ( tesla family )
            //! rather than as a field strength, so mScale carries the nu0 = 1/mu0
            //! conversion. Read by MaxwellFactory to warn about a ferro-bounded
            //! sideset, where B = mu0*H does not hold. Explicitly initialised:
            //! mUnits above is NOT, and cannot be used for this
            bool mAmplitudeIsFluxDensity = false ;

            Vector < id_t > mDomains;

            bool mIsThinShell ;

            DofManager * mDofMngr = nullptr;

            Dof * mAbstractDof = nullptr;

//-----------------------------------------------------------------------------

        public :
//-----------------------------------------------------------------------------

            PhysicalBoundaryCondition( const index_t aIndex ) ;

            ~PhysicalBoundaryCondition() ;

//-----------------------------------------------------------------------------

            index_t
            index() const ;

//-----------------------------------------------------------------------------

            void
            set_function( SourceFunction * aBCFunction ) ;

//-----------------------------------------------------------------------------

            void
            set_direction( real x, real y, real z ) ;

//-----------------------------------------------------------------------------

            void
            set_type( BoundaryConditionType aType ) ;

//-----------------------------------------------------------------------------

            void
            set_units( unit aUnits ) ;

//-----------------------------------------------------------------------------

            void
            set_amplitude_is_flux_density( const bool aFlag ) ;

//-----------------------------------------------------------------------------

            bool
            amplitude_is_flux_density() const ;

//-----------------------------------------------------------------------------

            void
            set_domains( Vector < id_t > aDomains, const bool aIsThinShell = false ) ;

//-----------------------------------------------------------------------------

            void
            set_field( DofManager * aField ) ;

//-----------------------------------------------------------------------------

            void
            set_scale( real aScale ) ;

//-----------------------------------------------------------------------------

            real &
            scale() ;

//-----------------------------------------------------------------------------

            SourceFunction *
            function() ;

//-----------------------------------------------------------------------------

            Vector< id_t > &
            domains() ;

//-----------------------------------------------------------------------------

            BoundaryConditionType
            type() ;

//-----------------------------------------------------------------------------

            bool
            is_thinshell() ;

//-----------------------------------------------------------------------------

            real
            value() ;

//-----------------------------------------------------------------------------

            void
            fix( real aValue ) ;

//-----------------------------------------------------------------------------

            Vector < real > &
            direction() ;

//-----------------------------------------------------------------------------

            unit
            units() ;

//-----------------------------------------------------------------------------

            void
            impose_bc( real aTime ) ;

//-----------------------------------------------------------------------------

            //! set the global-variable name ( factories only; empty = no
            //! global is published for this condition )
            void
            set_label( const string & aLabel ) ;

//-----------------------------------------------------------------------------

            const string &
            label() const ;

//-----------------------------------------------------------------------------

            //! write mValue into this condition's mesh global variable, if a
            //! label is set and the variable exists. The variables live on
            //! RANK 0 ONLY ( created there by the factories; Exodus and the
            //! memdump read them there ), so on workers the exists-guard
            //! makes this a no-op — for function-evaluated types called from
            //! impose_bc() on all ranks, and for circuit types written from
            //! fix() at circuit-solve time on rank 0. No worker code reads
            //! mesh globals; anyone adding such a consumer must add an
            //! explicit synch first. Deliberately a no-op before the
            //! factories have created the variables
            void
            update_global() ;

//-----------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_FEM_BOUNDARYCONDITION_HPP
