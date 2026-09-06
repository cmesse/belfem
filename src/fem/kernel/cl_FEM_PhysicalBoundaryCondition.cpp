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
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_Node.hpp"

namespace belfem
{
    namespace fem
    {

        PhysicalBoundaryCondition::PhysicalBoundaryCondition( const index_t aIndex ) :
            mIndex( aIndex )
        {

        }
//-----------------------------------------------------------------------------

        PhysicalBoundaryCondition::~PhysicalBoundaryCondition()
        {
            delete mFunction ;
        }

//-----------------------------------------------------------------------------

        index_t
        PhysicalBoundaryCondition::index() const
        {
            return mIndex ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_function( SourceFunction * aBCFunction )
        {
            mFunction = aBCFunction ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_direction(real x, real y, real z)
        {
            mDirection(0) = x;
            mDirection(1) = y;
            mDirection(2) = z;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_type( BoundaryConditionType aType)
        {
            mType = aType ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_units( unit aUnits)
        {
            mUnits = aUnits ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_amplitude_is_flux_density( const bool aFlag )
        {
            mAmplitudeIsFluxDensity = aFlag ;
        }

//-----------------------------------------------------------------------------

        bool
        PhysicalBoundaryCondition::amplitude_is_flux_density() const
        {
            return mAmplitudeIsFluxDensity ;
        }

        //-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_domains( Vector < id_t > aDomains, const bool aIsThinShell )
        {
            mIsThinShell = aIsThinShell ;
            mDomains = aDomains ;
        }

//-----------------------------------------------------------------------------

        SourceFunction *
        PhysicalBoundaryCondition::function()
        {
            return mFunction ;
        }

//-----------------------------------------------------------------------------

        Vector < id_t > &
        PhysicalBoundaryCondition::domains()
        {
            return mDomains ;
        }

//-----------------------------------------------------------------------------

        BoundaryConditionType
        PhysicalBoundaryCondition::type()
        {
            return mType ;
        }

//-----------------------------------------------------------------------------

        bool
        PhysicalBoundaryCondition::is_thinshell()
        {
            return mIsThinShell ;
        }

//-----------------------------------------------------------------------------

        real
        PhysicalBoundaryCondition::value()
        {
            return mValue ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::fix( real aValue )
        {
            mValue = aValue*mScale ;

            // circuit types receive their value through fix() from the
            // circuit loop ( rank 0 ); publish it to the mesh global here,
            // since impose_bc() skips these types
            this->update_global() ;
        }

//-----------------------------------------------------------------------------

        Vector < real > &
        PhysicalBoundaryCondition::direction()
        {
            return mDirection ;
        }

//-----------------------------------------------------------------------------

        unit
        PhysicalBoundaryCondition::units()
        {
            return mUnits ;
        }

//-----------------------------------------------------------------------------


        void
        PhysicalBoundaryCondition::set_field( DofManager * aField)
        {
            mDofMngr = aField ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_scale( real aScale )
        {
            mScale = aScale ;
        }

//-----------------------------------------------------------------------------

        real &
        PhysicalBoundaryCondition::scale()
        {
            return mScale ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::set_label( const string & aLabel )
        {
            mLabel = aLabel ;
        }

//-----------------------------------------------------------------------------

        const string &
        PhysicalBoundaryCondition::label() const
        {
            return mLabel ;
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::update_global()
        {
            if ( mLabel.size() > 0
                && mDofMngr != nullptr
                && mDofMngr->mesh()->global_variable_exists( mLabel ) )
            {
                mDofMngr->mesh()->global_variable_data( mLabel ) = mValue ;
            }
        }

//-----------------------------------------------------------------------------

        void
        PhysicalBoundaryCondition::impose_bc( real aTime )
        {
            BELFEM_ERROR(mDofMngr!= nullptr,"The boundary condition does not have a Dof Manager") ;

            if (mType !=  BoundaryConditionType::Bearing
                && mType != BoundaryConditionType::CircuitCurrent
                && mType != BoundaryConditionType::CircuitVoltage) mValue = mFunction->compute(aTime) * mScale ;

            // publish the imposed value to this condition's mesh global.
            // MUST sit before the switch: every handled case returns early.
            // Circuit types are published from fix() instead ( their mValue
            // here is whatever the circuit loop last set on this rank )
            if ( mType != BoundaryConditionType::CircuitCurrent
                && mType != BoundaryConditionType::CircuitVoltage )
            {
                this->update_global() ;
            }

            switch (mType)
            {
                case ( BoundaryConditionType::Bearing ) :
                case ( BoundaryConditionType::Current ) :
                case ( BoundaryConditionType::Voltage ) :
                case ( BoundaryConditionType::CircuitCurrent ) :
                case ( BoundaryConditionType::CircuitVoltage ) :
                {
                    //cases already handled in the Maxwell Factory or in the script
                    return ;
                }
                case ( BoundaryConditionType::Gauge ) :
                {
                    for(id_t tID : mDomains)
                    {
                        mDofMngr->bearing( tID )->impose_dirichlet( mFunction->compute(aTime) ) ;
                    }
                    return ;
                }
                case ( BoundaryConditionType::Dirichlet ) :
                {
                    for(id_t tID : mDomains)
                    {
                        mDofMngr->sideset( tID )->impose_dirichlet( mFunction->compute(aTime) ) ;
                    }
                    return ;
                }
                case ( BoundaryConditionType::Background ) : //Maxwell Specific
                {
                    // this BC is imposed weakly with its own stiffness matrix
                    return ;
                }
                case ( BoundaryConditionType::BackgroundDirichlet ) : //Maxwell Specific
                {
                    // NOTE, if this type is ever revived: unlike every branch
                    // above, the value below is taken straight from the source
                    // function and mScale is NOT applied. A background stated in
                    // tesla carries its nu0 = 1/mu0 conversion in mScale, so it
                    // would be silently dropped here and the field imposed a
                    // factor 8e5 too small. The type is refused at creation today
                    // ( MaxwellBoundaryConditionFactory, default case ), which is
                    // the only reason this is a comment and not a bug
                    for(id_t tID : mDomains)
                    {
                        for (mesh::Node * tNode : mDofMngr->sideset(tID)->nodes())
                        {
                            real tValue = -1.0*(mFunction->compute(aTime)*mDirection(0)*tNode->x()+
                                          mFunction->compute(aTime)*mDirection(1)*tNode->y()+
                                          mFunction->compute(aTime)*mDirection(2)*tNode->z());
                            if(!tNode->is_duplicate())
                            {
                                mDofMngr->dof(mDofMngr->calculate_dof_id(tNode,0))->fix(tValue) ;
                            }
                        }
                    }
                    return ;
                }
                default:
                {
                    BELFEM_ERROR(false, "Boundary condition imposition not implemented yet") ;
                }
            }
        }

//-----------------------------------------------------------------------------
    }
}
