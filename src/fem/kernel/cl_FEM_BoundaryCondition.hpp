//
// Created by christian on 9/29/21.
//

#ifndef BELFEM_CL_FEM_BOUNDARYCONDITION_HPP
#define BELFEM_CL_FEM_BOUNDARYCONDITION_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_FEM_BoundaryConditionImposing.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_SideSet.hpp"

#include "cl_Node.hpp"
#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
        enum class BoundaryConditionScaling
        {
            Constant    = 0,
            Ramp        = 1,
            Sigmoid     = 2,
            Sine        = 3,
            Square      = 4,
            Triangle    = 5,
            Sawtooth    = 6,
            UNDEFINED   = 7
        };

        enum class BoundaryConditionPhysics
        {
            Default     = 0,
            Current     = 1,
            Magfield    = 2,
            UNDEFINED   = 3
        };



//-----------------------------------------------------------------------------

        BoundaryConditionScaling
        bc_scaling_type( const string & aString );

        class DofManager ;
//-----------------------------------------------------------------------------

        class BoundaryCondition
        {
//-----------------------------------------------------------------------------
        protected:
//-----------------------------------------------------------------------------

            const BoundaryConditionPhysics     mPhysics  ;
            const BoundaryConditionImposing    mImposing ;

            const string mLabel ;

            DofManager * mParent = nullptr ;

            // fields
            Cell< string > mFields ;

            // - - - - - - - - - - - - - - - - - - - -
            // cell with blocks, if selected
            Cell< Block * > mBlocks ;

            // cell with sidesets, if selected
            Cell< SideSet * > mSideSets ;

            // cell with dofs this BC is related to

            // list with nodes that are associated with this bc
            Cell< mesh::Node * > mNodes ;

            // list with dof types
            Vector< uint > mDofTypes ;

            // list to be used in preprocessing
            Vector< id_t > mSideSetIDs ;
            Vector< id_t > mBlockIDs ;

            // the scaling function for the amplitude
            real
            ( BoundaryCondition::*mFunScale )( const real aTime ) const ;


            // this is a data container that can be accessed by
            // FEM Blocks and Sidesets
            Vector< real > mData ;

            // penalty factor. makes only sense for a weak BC
            real mPenalty = 1.0 ;

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------



            BoundaryConditionScaling mScaling
                = BoundaryConditionScaling::UNDEFINED ;

            // time, if this is a ramp
            real mPeriod    = BELFEM_QUIET_NAN ;

            // the phase, if this is a sine
            real mPhase    = 0.0 ;

            // the time offset for the phase
            real mTimeOffset   = 0.0 ;

            // the frequency = inverse of period
            real mFrequency = BELFEM_QUIET_NAN ;

            // the angular frequency if this is a sine
            real mOmega = BELFEM_QUIET_NAN ;

            // fuziness factor, only relevant for sigmoid function
            real mFuziness = BELFEM_QUIET_NAN;


//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            BoundaryCondition(
                    const BoundaryConditionPhysics   aPhysics,
                    const BoundaryConditionImposing  aImposing,
                    const string                  aLabel="boundaryCondition" );

            virtual
            ~BoundaryCondition() = default ;

//-----------------------------------------------------------------------------

            virtual void
            set_dof_manager( DofManager * aDofManager );


//-----------------------------------------------------------------------------

            /**
             * return the imposing type
             */
            BoundaryConditionImposing
            imposing() const ;

//-----------------------------------------------------------------------------

            /**
             * return the physics type
             */
            BoundaryConditionPhysics
            physics() const ;

//-----------------------------------------------------------------------------

            /**
             * return the scaling type
             */
            BoundaryConditionScaling
            scaling() const ;

//-----------------------------------------------------------------------------

            /**
             * return the parent
             */
            DofManager *
            parent() ;

//-----------------------------------------------------------------------------

            /**
             * to be used if dof manager has not been set
             */
            void
            set_sidesets( const Vector< id_t > & aSideSetIDs );

//-----------------------------------------------------------------------------

            /**
             * set the penalty factor. It is only used for a weak BC
             */
            void
            set_penalty( const real aValue );

//-----------------------------------------------------------------------------

            /**
             * to be used if dof manager has not been set
             */
            void
            set_blocks( const Vector< id_t > & aBlockIDs );

//-----------------------------------------------------------------------------

            /**
             * calls the iwg and runs the computation
             * @return
             */
            virtual void
            compute( const real aTimeStamp = 0.0 );

//-----------------------------------------------------------------------------

            /**
             * called by dof manager
             */
            void
            update_dof_types();

//-----------------------------------------------------------------------------

            /**
             * return the label of this bc
             */
             const string &
             label() const ;

//-----------------------------------------------------------------------------

            /**
             *
             * @param aType
             * @param aPeriod
             * @param aPhase      time offset
             * @param aFuzzyness  only for sigmoid function
             */
            void
            set_scaling(
                    const BoundaryConditionScaling aType,
                    const real        aPeriod=BELFEM_REAL_MAX,
                    const real        aPhase=0.0,
                    const real        aFuzzyness=50 );

//-----------------------------------------------------------------------------

            /**
             * return the list of blocks
             */
             const Vector< id_t > &
             blocks() const ;

//-----------------------------------------------------------------------------

            /**
             * return the list of sidesets
             */
            const Vector< id_t > &
            sidesets() const ;

//-----------------------------------------------------------------------------

            /**
             * access the data container
             */
            const Vector< real > &
            data() const ;

//-----------------------------------------------------------------------------

            /**
             * access one value of the data container
             */
            real
            data( const uint aIndex ) const ;

//-----------------------------------------------------------------------------

            /**
             * return the penalty factor
             */
            real
            penalty() const ;

//-----------------------------------------------------------------------------
        protected:
//-----------------------------------------------------------------------------

            /**
             * check if dof manager has been linked
             */
            void
            check_if_linked_to_dof_manager();

//-----------------------------------------------------------------------------

            virtual void
            link_sidesets(  DofManager * aDofManager, const Vector< id_t > & aIDs );

//-----------------------------------------------------------------------------

            virtual void
            link_blocks(  DofManager * aDofManager, const Vector< id_t > & aIDs );

//-----------------------------------------------------------------------------
            /**
             * called by constructor
             */
            void
            collect_nodes();

//-----------------------------------------------------------------------------

            real
            scale( const real aTime ) const ;

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            real
            scale_constant( const real aTime ) const ;

//-----------------------------------------------------------------------------

            real
            scale_ramp( const real aTime ) const ;

//-----------------------------------------------------------------------------

            // see doi:10.1016/j.fss.2005.02.016
            real
            scale_sigmoid( const real aTime ) const ;

//-----------------------------------------------------------------------------

            real
            scale_sine( const real aTime ) const ;

//-----------------------------------------------------------------------------

            real
            scale_sawtooth( const real aTime ) const ;

//----------------------------------------------------------------------------

            real
            scale_square( const real aTime ) const ;

//----------------------------------------------------------------------------

            real
            scale_triangle( const real aTime ) const ;

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------

        inline DofManager *
        BoundaryCondition::parent()
        {
            return mParent ;
        }

//-----------------------------------------------------------------------------

        inline BoundaryConditionImposing
        BoundaryCondition::imposing() const
        {
            return mImposing ;
        }

//-----------------------------------------------------------------------------

        inline const string &
        BoundaryCondition::label() const
        {
            return mLabel ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale( const real aTime ) const
        {
            return (this->*mFunScale)( aTime );
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_constant( const real aTime ) const
        {
            return 1.0 ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_sigmoid( const real aTime ) const
        {

            real a = 0.5*mPeriod ;  // upper limit
            real b = mFuziness/mPeriod ;  // fuziness parameter
            real c = 0.5*mPeriod ;  // slope
            real d = 1.0/mFuziness ;
            real x = aTime - mTimeOffset ;

            return d*std::log(((1.+std::exp(b*(x-c+a)))/
                (1.+std::exp(b*(x-c-a)))));
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_ramp( const real aTime ) const
        {
            return aTime < mPeriod ? aTime / mPeriod : 1.0 ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_sine( const real aTime ) const
        {
            return std::sin( mOmega * aTime + mPhase );
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_sawtooth( const real aTime ) const
        {
            return ( aTime + mTimeOffset ) * mFrequency
                - std::floor( ( aTime + mTimeOffset ) * mFrequency );
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_square( const real aTime ) const
        {
            return this->scale_sawtooth( aTime ) < 0.5 ? 1.0 : -1.0 ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::scale_triangle( const real aTime ) const
        {
            real tX = 4.0 * this->scale_sawtooth( aTime ) ;
            return tX < 2.0 ? tX - 1.0 : 3.0 - tX ;
        }

//-----------------------------------------------------------------------------

        inline BoundaryConditionScaling
        BoundaryCondition::scaling() const
        {
            return mScaling ;
        }

//-----------------------------------------------------------------------------

        inline BoundaryConditionPhysics
        BoundaryCondition::physics() const
        {
            return mPhysics ;
        }

//-----------------------------------------------------------------------------

        inline const Vector< id_t > &
        BoundaryCondition::blocks() const
        {
            return mBlockIDs ;
        }

//-----------------------------------------------------------------------------

        inline const Vector< id_t > &
        BoundaryCondition::sidesets() const
        {
            return mSideSetIDs ;
        }

//-----------------------------------------------------------------------------

        inline const Vector< real > &
        BoundaryCondition::data() const
        {
            return mData ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::data( const uint aIndex ) const
        {
            return mData( aIndex ) ;
        }

//-----------------------------------------------------------------------------

        inline real
        BoundaryCondition::penalty() const
        {
            return mPenalty ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_BOUNDARYCONDITION_HPP
