//
// Created by christian on 10/8/21.
//

#ifndef BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONMAGFIELD_HPP
#define BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONMAGFIELD_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_BoundaryCondition.hpp"
#include "en_FEM_MagfieldBcType.hpp"

namespace belfem
{
    namespace fem
    {

//-----------------------------------------------------------------------------

        class MaxwellBoundaryConditionMagfield : public BoundaryCondition
        {

            MagfieldBcType mSubType = MagfieldBcType::UNDEFINED ;

            // the amplitude as vector value
            Vector< real > mAmplitude ;

            // indices with nodes on which vector field is imposed
            Vector< index_t > mNodeIndices ;

            Cell< string > mBoundaryFields ;
            Cell< string > mDofFields ;

            // flag telling if this is a h-phi formulation
            bool mIsHPhi = false ;

            // the dofs are owned and destroyed by the dof manager
            // this cell is needed for the h-phi formulation
            Cell< Dof * > mFixedDofs ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            MaxwellBoundaryConditionMagfield(
                    const BoundaryConditionImposing aImposing,
                    const string             & aLabel );

//-----------------------------------------------------------------------------

            ~MaxwellBoundaryConditionMagfield() = default ;

//-----------------------------------------------------------------------------

            void
            set( const BoundaryConditionScaling aType,
                      const Vector< real > & aAmplitude,
                      const real aPeriod=BELFEM_REAL_MAX,
                      const real aPhase=0.0 );

//-----------------------------------------------------------------------------

            void
            set( const MagfieldBcType aType );

//-----------------------------------------------------------------------------

            const Vector< real > &
            values() const ;

//-----------------------------------------------------------------------------

            void
            compute( const real aTime );

//-----------------------------------------------------------------------------

            MagfieldBcType
            subtype() const ;

//-----------------------------------------------------------------------------
        protected:
//-----------------------------------------------------------------------------

            void
            link_sidesets(  DofManager * aDofManager, const Vector< id_t > & aIDs );

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------

        inline const Vector< real > &
        MaxwellBoundaryConditionMagfield::values() const
        {
            return mAmplitude ;
        }

//-----------------------------------------------------------------------------

        MagfieldBcType
        inline MaxwellBoundaryConditionMagfield::subtype() const
        {
            return mSubType ;
        }

//-----------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace belfem */
#endif //BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONMAGFIELD_HPP
