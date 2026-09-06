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

#include "typedefs.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "fn_dot.hpp"
#include <algorithm>
#include <cmath>
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_DofManager.hpp"

#include "fn_entity_type.hpp"
#include "fn_norm.hpp"

    namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Timestep::IWG_Timestep(
                const IwgType aType,
                const ModelDimensionality aDimensionality,
                const IwgMode aMode,
                const SymmetryMode aSymmetryMode,
                const DofMode      aDofMode,
                const SideSetDofLinkMode aSideSetDofLinkMode ) :
                IWG( aType, aDimensionality, aMode, aSymmetryMode, aDofMode, aSideSetDofLinkMode )
        {
            mNumberOfRhsCols = 1 ;
            mAlpha = 1.0 ;
            mBeta.set_size( 5, BELFEM_QUIET_NAN );
            mH.set_size( 5, 0.0 );
            mTimeStepMatrices = new TimestepMatrices() ;

            // safe default so the dispatch pointer is never null; drivers and
            // the Controller may reconfigure this later ( the old-timestep
            // storage is then recreated for the new order )
            this->set_timestepping_method( EulerMethod::BackwardDifference1 );
        }

        IWG_Timestep::~IWG_Timestep()
        {
            delete mTimeStepMatrices ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Timestep::set_timestepping_method( const EulerMethod aMethod, const bool aHaveStiffness )
        {
            // remember the method
            mMethod = aMethod ;
            mHaveStiffness = aHaveStiffness ;

            if( aHaveStiffness )
            {
                switch( aMethod )
                {
                    case( EulerMethod::Static ) :
                    {
                        mTimestep = & IWG_Timestep::static_subfield ;
                        mOrder = 0 ;
                        break ;
                    }
                    case( EulerMethod::ForwardExplicit ) :
                    {
                        mTimestep = & IWG_Timestep::explicit_euler ;
                        mOrder = 1 ;
                        break;
                    }
                    case( EulerMethod::CrankNicolson ) :
                    case( EulerMethod::Galerkin ) :
                    {
                        // the Newton correction assemble_dJdx is exact for the
                        // BDF family only ( J = alpha*M + dt*K ); these theta
                        // schemes would get a wrong tangent, so they are
                        // disabled until dJdx learns their contract
                        mTimestep = nullptr ;
                        BELFEM_ERROR( false,
                            "Crank-Nicolson and Galerkin timestepping are disabled: "
                            "the Newton tangent assembly supports the BDF family only" );
                        break ;
                    }
                    case( EulerMethod::BackwardDifference1 ) :
                    {
                        mTimestep = & IWG_Timestep::bdf1 ;
                        mOrder = 1 ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference2 ) :
                    {
                        mOrder = 2 ;
                        mTimestep = & IWG_Timestep::bdf2 ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference3 ) :
                    {
                        mOrder = 3 ;
                        mTimestep = & IWG_Timestep::bdf3 ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference4 ) :
                    {
                        mOrder = 4 ;
                        mTimestep = & IWG_Timestep::bdf4 ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference5 ) :
                    {
                        mOrder = 5 ;
                        mTimestep = & IWG_Timestep::bdf5 ;
                        break ;
                    }
                    case( EulerMethod::Derivative ) :
                    {
                        mTimestep = & IWG_Timestep::derivative ;
                        break ;
                    }
                    case( EulerMethod::StiffnessOnly ) :
                    {
                        mTimestep = & IWG_Timestep::stiffness_matrix_only ;
                        break ;
                    }
                    case( EulerMethod::MassOnly ) :
                    {
                        mTimestep = & IWG_Timestep::mass_matrix_only ;
                        break ;
                    }

                    default:
                    {
                        mTimestep = nullptr ;
                        BELFEM_ERROR( false, "Invalid euler method");
                        break ;
                    }
                }
            }
            else
            {
                switch( aMethod )
                {
                    case( EulerMethod::Static ) :
                    {
                        mTimestep = nullptr ;
                        BELFEM_ERROR( false,
                                      "Invalid euler method: must have a stiffness if field is static");
                        break ;
                    }
                    case( EulerMethod::CrankNicolson ) :
                    case( EulerMethod::Galerkin ) :
                    {
                        // do not silently alias onto bdf1_nok — these schemes
                        // are disabled, see the stiffness branch above
                        mTimestep = nullptr ;
                        BELFEM_ERROR( false,
                            "Crank-Nicolson and Galerkin timestepping are disabled: "
                            "the Newton tangent assembly supports the BDF family only" );
                        break ;
                    }
                    case( EulerMethod::ForwardExplicit ) :
                    case( EulerMethod::BackwardDifference1 ) :
                    {
                        mOrder = 1 ;
                        mTimestep = & IWG_Timestep::bdf1_nok ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference2 ) :
                    {
                        mOrder = 2 ;
                        mTimestep = & IWG_Timestep::bdf2_nok ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference3 ) :
                    {
                        mOrder = 3 ;
                        mTimestep = & IWG_Timestep::bdf3_nok ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference4 ) :
                    {
                        mOrder = 4 ;
                        mTimestep = & IWG_Timestep::bdf4_nok ;
                        break ;
                    }
                    case( EulerMethod::BackwardDifference5 ) :
                    {
                        mOrder = 5 ;
                        mTimestep = & IWG_Timestep::bdf5_nok ;
                        break ;
                    }
                    case( EulerMethod::Derivative ) :
                    {
                        mTimestep = & IWG_Timestep::derivative_no_stiffness ;
                        break ;
                    }
                    default:
                    {
                        mTimestep = nullptr ;
                        BELFEM_ERROR( false, "Invalid euler method");
                        break ;
                    }
                }
            }

            // default dispatch; compute_bdf_coefficients may override this
            // during the startup ramp of a BDF-p run
            mTimestepActive = mTimestep ;
            mCoeffsDirty = true ;

            // if the field is already linked, the old-timestep storage was
            // sized with the previous order and must be recreated, otherwise
            // shift_fields() indexes past the end of mFieldData
            if ( mField != nullptr )
            {
                this->create_old_dof_fields() ;

                // same hazard on the calculator side: the qold tables were
                // sized by allocate() under the previous order ( the coupled
                // main() configures the thermal method after its dof manager
                // exists ), so an order change must rebuild them — otherwise
                // Calculator::qold( aStep ) reads past the table on the
                // first higher-order step
                DofManager * tDofMgr = reinterpret_cast< DofManager * >( mField );
                for ( Block * tBlock : tDofMgr->blocks() )
                {
                    if ( tBlock->calculator() != nullptr )
                    {
                        tBlock->calculator()->init_qold_table();
                    }
                }
                for ( SideSet * tSideSet : tDofMgr->sidesets() )
                {
                    if ( tSideSet->calculator() != nullptr )
                    {
                        tSideSet->calculator()->init_qold_table();
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::history_field_labels(
                Cell< string > & aParents,
                Cell< string > & aHistory,
                const uint       aDepth ) const
        {
            // SINGLE SOURCE for the numbered-history naming convention:
            // must stay in lockstep with create_old_dof_fields() below
            // ( "%s%u" over mDofLabels — the fields shift_fields rotates ).
            // aDepth = 0 means the full mOrder; callers that only touch
            // the levels a resumed step reads pass that depth instead —
            // handing a DEEPER, never-filled level to a consumer that
            // walks entity indices is a bounds fault on EDGE/FACE fields,
            // whose vectors stay empty until something sizes them
            aParents.clear() ;
            aHistory.clear() ;

            const uint tDepth = aDepth == 0 ? mOrder : std::min( aDepth, mOrder ) ;

            for ( const string & tLabel : mDofLabels )
            {
                // mDofLabels repeats a label once per multiplicity
                // ( higher-order edge/face dofs ); the FIELD is one object,
                // so emit each name once
                bool tSeen = false ;
                for ( const string & tKnown : aParents )
                {
                    if ( tKnown == tLabel )
                    {
                        tSeen = true ;
                        break ;
                    }
                }
                if ( tSeen ) continue ;

                aParents.push( tLabel );
                for ( uint k = 0; k < tDepth; ++k )
                {
                    aHistory.push( tLabel + std::to_string( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::create_old_dof_fields()
        {
            string tFormat  = "%s%u" ;
            uint tNumFields = mDofLabels.size() ;

            mFieldData.set_size( mOrder+1, {} );

            // get the mesh
            uint c = 0 ;

            BELFEM_ASSERT( mMesh != nullptr, "Mesh not set" );

            // placeholder for old timesteps
            for( uint k=0; k<mOrder; ++k )
            {

                Cell< Vector< real > * > & tFieldData( mFieldData( ++c  ) );
                tFieldData.set_size( tNumFields, nullptr );

                for( uint f=0; f<tNumFields; ++f )
                {
                    string tOldLabel = sprint( tFormat.c_str(), mDofLabels( f ).c_str(), k );

                    EntityType tType = entity_type( mDofLabels( f ) );

                    // create the new field if it doesn't exist already
                    Vector< real > & tField = mMesh->field_exists( tOldLabel ) ? mMesh->field_data( tOldLabel ) :
                        mMesh->create_field(
                        tOldLabel, tType );

                    // don't write this field
                    mMesh->field( tOldLabel )->set_write_to_file_flag( false );

                    // add pointer to data container
                    tFieldData( f ) = &tField ;
                }
            }

            // placeholder for current timestep
            Cell< Vector< real > * > & tFieldData( mFieldData( 0  ) );
            tFieldData.set_size( tNumFields, nullptr );

            for( uint f=0; f<tNumFields; ++f )
            {

                EntityType tType = entity_type( mDofLabels( f ) );

                // create the new field if it doesn't exist already
                Vector< real > & tField = mMesh->field_exists( mDofLabels( f ) ) ? mMesh->field_data( mDofLabels( f ) ) :
                    mMesh->create_field(
                    mDofLabels( f ), tType );


                // add pointer to data container
                tFieldData( f ) = &tField ;
            }

        }

        void
        IWG_Timestep::shift_fields()
        {
            // initialize the dof manager if it hasn't been done already
            mField->initialize();

            uint tNumFields = mDofLabels.size() ;

            for( int k=mOrder-1; k>=0; k-- )
            {
                Cell< Vector< real > * > & tTargets( mFieldData( k+1 ) );
                Cell< Vector< real > * > & tSources( mFieldData( k ) );

                for( uint f=0; f<tNumFields; ++f )
                {
                    Vector< real > & tTarget = *tTargets( f );
                    Vector< real > & tSource = *tSources( f );

                    // shift values
                    tTarget = tSource ;
                }
            }

            // shift timesteps; remember the value that falls off the end so
            // reset_fields can undo this shift completely ( BDF5 reads mH(3) )
            mHDropped = mH( 3 );
            mH( 3 ) = mH( 2 );
            mH( 2 ) = mH( 1 );
            mH( 1 ) = mH( 0 );
            mH( 0 ) = mDeltaTime ;

            // a new step has started: recompute the BDF coefficients lazily
            // once delta_time has been set for the upcoming step
            ++mStepCount ;
            mCoeffsDirty = true ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_bdf_coefficients()
        {
            // only BDF2-5 use the alpha/beta coefficients and the startup
            // ramp; every other scheme scales M by one and must run exactly
            // the configured method ( gating on mMethod also protects the
            // temporary MassOnly/StiffnessOnly switch in EigenValues from
            // stale mOrder state )
            switch( mMethod )
            {
                case( EulerMethod::BackwardDifference2 ) :
                case( EulerMethod::BackwardDifference3 ) :
                case( EulerMethod::BackwardDifference4 ) :
                case( EulerMethod::BackwardDifference5 ) :
                {
                    break ;
                }
                default :
                {
                    mAlpha = 1.0 ;
                    mOrderActive = 1 ;
                    mTimestepActive = mTimestep ;
                    return ;
                }
            }

            // startup ramp: during the first steps of a BDF-p run only
            // min( p, steps taken ) history states exist, so we run the
            // highest order the history supports ( step 1: BDF1,
            // step 2: BDF2, ... ) and never read an unpopulated mH slot.
            // The rule itself lives in ramped_order(), so the log can report
            // the same number before this lazy recompute has run.
            const uint tOrder = this->ramped_order() ;
            mOrderActive = tOrder ;

            BELFEM_ASSERT( mDeltaTime > 0.0,
                "delta_time must be set before the BDF coefficients are computed" );

#if !defined( NDEBUG ) || defined( DEBUG )
            // the coefficient formulas read tOrder-1 previous step sizes
            for ( uint k=0; k+1<tOrder; ++k )
            {
                BELFEM_ASSERT( mH( k ) > 0.0,
                    "Invalid BDF history step size mH(%u)", ( unsigned int ) k );
            }
#endif

            // update BDF coefficients using current mDeltaTime (upcoming step)
            // and mH values (previous steps, already shifted)
            switch( tOrder )
            {
                case( 2 ) :
                {
                    this->compute_bdf_coeffs_2() ;
                    break ;
                }
                case( 3 ) :
                {
                    this->compute_bdf_coeffs_3() ;
                    break ;
                }
                case( 4 ) :
                {
                    this->compute_bdf_coeffs_4() ;
                    break ;
                }
                case( 5 ) :
                {
                    this->compute_bdf_coeffs_5() ;
                    break ;
                }
                default:
                {
                    // all schemes of order <= 1 scale the mass matrix by one
                    mAlpha = 1.0 ;
                    break ;
                }
            }

            // select the scheme that runs this step
            if ( tOrder == mOrder )
            {
                mTimestepActive = mTimestep ;
            }
            else if ( mHaveStiffness )
            {
                switch( tOrder )
                {
                    case( 2 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf2 ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf3 ;
                        break ;
                    }
                    case( 4 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf4 ;
                        break ;
                    }
                    default :
                    {
                        mTimestepActive = & IWG_Timestep::bdf1 ;
                        break ;
                    }
                }
            }
            else
            {
                switch( tOrder )
                {
                    case( 2 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf2_nok ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf3_nok ;
                        break ;
                    }
                    case( 4 ) :
                    {
                        mTimestepActive = & IWG_Timestep::bdf4_nok ;
                        break ;
                    }
                    default :
                    {
                        mTimestepActive = & IWG_Timestep::bdf1_nok ;
                        break ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::reset_fields()
        {
            uint tNumFields = mDofLabels.size() ;

            // reverse the shift: move all slots back by one
            for( uint k=0; k<mOrder; ++k )
            {
                Cell< Vector< real > * > & tTargets( mFieldData( k ) );
                Cell< Vector< real > * > & tSources( mFieldData( k+1 ) );

                for( uint f=0; f<tNumFields; ++f )
                {
                    Vector< real > & tTarget = *tTargets( f );
                    Vector< real > & tSource = *tSources( f );

                    // un-shift values
                    tTarget = tSource ;
                }
            }

            // restore the step size of the last completed step, so that the
            // shift_fields of the retry pushes it back onto mH( 0 ) instead
            // of the size of the abandoned step
            mDeltaTime = mH( 0 );

            // reverse the timestep shift, including the slot the shift
            // dropped off the end ( read by BDF5 )
            mH( 0 ) = mH( 1 );
            mH( 1 ) = mH( 2 );
            mH( 2 ) = mH( 3 );
            mH( 3 ) = mHDropped ;

            // roll back the step counter for the startup order ramp
            if ( mStepCount > 0 )
            {
                --mStepCount ;
            }
            mCoeffsDirty = true ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::make_savepoint()
        {
            uint tNumFields = mDofLabels.size() ;

            // ( re- ) allocate the snapshot storage; the size can change when
            // set_timestepping_method is called after the field was linked
            if ( mFieldSnapshot.size() != mOrder + 1 )
            {
                mFieldSnapshot.set_size( mOrder+1, {} );

                for( uint k=0; k<=mOrder; ++k )
                {
                    mFieldSnapshot( k ).set_size( tNumFields, {} );
                }

                mHSnapshot.set_size( mH.length(), 0.0 );
            }

            for( uint k=0; k<=mOrder; ++k )
            {
                Cell< Vector< real > * > & tSources( mFieldData( k ) );
                Cell< Vector< real > >   & tTargets( mFieldSnapshot( k ) );

                for( uint f=0; f<tNumFields; ++f )
                {
                    tTargets( f ) = *tSources( f );
                }
            }

            mHSnapshot         = mH ;
            mDeltaTimeSnapshot = mDeltaTime ;
            mHDroppedSnapshot  = mHDropped ;
            mStepCountSnapshot = mStepCount ;
            mHaveSavepoint     = true ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::restore_savepoint()
        {
            BELFEM_ERROR( mHaveSavepoint, "no savepoint to restore" );

            uint tNumFields = mDofLabels.size() ;

            for( uint k=0; k<=mOrder; ++k )
            {
                Cell< Vector< real > >   & tSources( mFieldSnapshot( k ) );
                Cell< Vector< real > * > & tTargets( mFieldData( k ) );

                for( uint f=0; f<tNumFields; ++f )
                {
                    *tTargets( f ) = tSources( f );
                }
            }

            mH          = mHSnapshot ;
            mDeltaTime  = mDeltaTimeSnapshot ;
            mHDropped   = mHDroppedSnapshot ;
            mStepCount  = mStepCountSnapshot ;
            mCoeffsDirty = true ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::save_history_state(
            Vector< real > & aH,
            uint           & aStepCount,
            real           & aLastDeltaTime ) const
        {
            aH = mH ;
            aStepCount = mStepCount ;
            aLastDeltaTime = mDeltaTime ;
        }

//------------------------------------------------------------------------------

        bool
        IWG_Timestep::restore_history_state(
            const Vector< real > & aH,
            const uint             aStepCount,
            const real             aLastDeltaTime )
        {
            // a zero step count behind a populated history would silently
            // re-anchor the ramp at BDF1 — the exact failure this API removes
            if ( aStepCount == 0 )
            {
                return false ;
            }

            if ( aH.length() != mH.length() )
            {
                return false ;
            }

            if ( ! std::isfinite( aLastDeltaTime ) || aLastDeltaTime <= 0.0 )
            {
                return false ;
            }

            for ( uint k = 0 ; k < aH.length() ; ++k )
            {
                if ( ! std::isfinite( aH( k ) ) || aH( k ) < 0.0 )
                {
                    return false ;
                }
            }

            // the slots a resumed full-order step and its possible retry
            // read must be positive; mH( 4 ) is the never-written spare and
            // zero there is normal
            const uint r = std::min( { aStepCount - 1,
                                       static_cast< uint >( mOrder - 1 ),
                                       4u } ) ;
            for ( uint k = 0 ; k < r ; ++k )
            {
                if ( ! ( aH( k ) > 0.0 ) )
                {
                    return false ;
                }
            }

            mH = aH ;

            // the ramp saturates at mOrder; capping avoids carrying an
            // arbitrarily large counter for no information
            mStepCount = std::min( aStepCount, static_cast< uint >( mOrder ) ) ;

            // the LAST COMPLETED step size: the first shift pushes this into
            // mH( 0 ), reconstructing the uninterrupted history exactly
            mDeltaTime = aLastDeltaTime ;

            mCoeffsDirty = true ;

            return true ;
        }

//------------------------------------------------------------------------------

        TimestepMatrices *
        IWG_Timestep::matrices()
        {
            return mTimeStepMatrices ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_mkf( Element * aElement )
        {
            BELFEM_ERROR( false, "compute_mkf not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Timestep::collect_qhist()
        {
            // qhist = beta0*q0 - beta1*q1 + beta2*q2 - beta3*q3 + beta4*q4,
            // truncated to mOrderActive ( which follows the startup ramp, not
            // mOrder ). The alternating signs are applied HERE — mBeta holds
            // positive magnitudes — which is what makes this the single source
            // of truth for the pattern. Full contract, including the lifetime
            // of the returned reference, is on the declaration in the header.
            Vector< real > & aSwap = mCalc->qswap();

            switch ( mOrderActive )
            {
                case( 2 ) :
                {
                    aSwap.vector_data()  = mBeta( 0 ) * mCalc->qold(0).vector_data();
                    aSwap.vector_data() -= mBeta( 1 ) * mCalc->qold(1).vector_data();
                    break ;
                }
                case( 3 ) :
                {
                    aSwap.vector_data()  = mBeta( 0 ) * mCalc->qold(0).vector_data();
                    aSwap.vector_data() -= mBeta( 1 ) * mCalc->qold(1).vector_data();
                    aSwap.vector_data() += mBeta( 2 ) * mCalc->qold(2).vector_data();
                    break ;
                }
                case( 4 ) :
                {
                    aSwap.vector_data()  = mBeta( 0 ) * mCalc->qold(0).vector_data();
                    aSwap.vector_data() -= mBeta( 1 ) * mCalc->qold(1).vector_data();
                    aSwap.vector_data() += mBeta( 2 ) * mCalc->qold(2).vector_data();
                    aSwap.vector_data() -= mBeta( 3 ) * mCalc->qold(3).vector_data();
                    break ;
                }
                case( 5 ) :
                {
                    aSwap.vector_data()  = mBeta( 0 ) * mCalc->qold(0).vector_data();
                    aSwap.vector_data() -= mBeta( 1 ) * mCalc->qold(1).vector_data();
                    aSwap.vector_data() += mBeta( 2 ) * mCalc->qold(2).vector_data();
                    aSwap.vector_data() -= mBeta( 3 ) * mCalc->qold(3).vector_data();
                    aSwap.vector_data() += mBeta( 4 ) * mCalc->qold(4).vector_data();
                    break ;
                }
                default :
                {
                    // order <= 1 : the history is the last step
                    aSwap.vector_data() = mCalc->qold(0).vector_data();
                    break ;
                }
            }
            return aSwap ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_jacobian_and_rhs( Element        * aElement,
                                                Matrix< real > & aJacobian,
                                                Vector< real > & aRHS )
        {

            // lazy update, flagged by shift_fields / reset_fields; by the time
            // the first element of a step is assembled, delta_time is already
            // set for the upcoming step
            if ( mCoeffsDirty )
            {
                this->compute_bdf_coefficients() ;
                mCoeffsDirty = false ;
            }

            aJacobian.fill(0.0) ;
            aRHS.fill(0.0) ;

            this->compute_mkf(aElement);

            //mTimeStepMatrices->assemble_J( mDeltaTime ) ;
            mTimeStepMatrices->assemble_dJdx( mDeltaTime, mAlpha ) ;

            this->compute_timestep( aJacobian, aRHS ) ;
            mCalc->K() = mTimeStepMatrices->K()*mDeltaTime ;

        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_timestep(
            Matrix< real > & aJ,
            Vector< real > & aRHS )
        {
            BELFEM_ASSERT( mTimestepActive != nullptr,
                "Timestepping method not set" );

            (this->*mTimestepActive )( aJ, aRHS );
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::set_field( DofManagerBase * aField )
        {
            IWG::set_field( aField );
            this->create_old_dof_fields() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::explicit_euler( Matrix< real > & aJ,
                                      Vector< real > & aRHS  )
        {
            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();
            aRHS.vector_data() += ( aJ.matrix_data() - mTimeStepMatrices->K().matrix_data()*mDeltaTime ) * mCalc->qold();
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::crank_nicolson( Matrix< real > & aJ,
                                      Vector< real > & aRHS  )
        {

            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;
            aRHS.vector_data() += ( aJ.matrix_data() - mTimeStepMatrices->K().matrix_data()*0.5 * mDeltaTime ) * mCalc->qold().vector_data() ;
            aJ.matrix_data() += mTimeStepMatrices->K().matrix_data()*0.5 * mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::galerkin( Matrix< real > & aJ,
                                  Vector< real > & aRHS )
        {

            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;
            aRHS.vector_data() += ( aJ.matrix_data() - mTimeStepMatrices->K().matrix_data()*mDeltaTime / 3) * mCalc->qold().vector_data();
            aJ.matrix_data() += mTimeStepMatrices->K().matrix_data()*mDeltaTime / 3 ;
            aJ.matrix_data() += mTimeStepMatrices->K().matrix_data()*mDeltaTime / 3 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::mass_matrix_only(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();
            aRHS.fill(0.0) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::stiffness_matrix_only(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->K().matrix_data();
            aRHS.fill(0.0) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf1(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {

            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();
            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // Update b: b = h * f + M * q0
            aRHS.vector_data() += aJ.matrix_data() * mCalc->qold().vector_data() ;

            // Update A: A = M + h * K
            aJ.matrix_data() += mTimeStepMatrices->K().matrix_data()*mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf2(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history ( beta0*q0 - beta1*q1 ), single source
            // of truth shared with the nonlinear-mass producers
            const Vector<real> &tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // M = alpha * M + h * K
            aJ.matrix_data() *= mAlpha ;
            aJ += mTimeStepMatrices->K().matrix_data()*mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf3(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // M = alpha * M + h * K
            aJ.matrix_data() *= mAlpha ;
            aJ += mTimeStepMatrices->K().matrix_data()*mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf4(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 - beta3 * q3 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // M = alpha * M + h * K
            aJ.matrix_data() *= mAlpha ;
            aJ += mTimeStepMatrices->K().matrix_data()*mDeltaTime ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf5(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            // Scale stiffness matrix K and forcing vector f by h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 - beta3 * q3 + beta4 * q4 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // M = alpha * M + h * K
            aJ.matrix_data() *= mAlpha ;
            aJ += mTimeStepMatrices->K().matrix_data()*mDeltaTime ;
        }

 //------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf1_nok(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {

            // Scaling forcing terms by h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // Update b: b = h * f + M * q0
            aRHS.vector_data() += aJ.matrix_data() * mCalc->qold().vector_data() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf2_nok(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            // Scaling forcing terms by * h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history ( beta0*q0 - beta1*q1 ), single source
            // of truth shared with the nonlinear-mass producers
            const Vector<real> &tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // M = alpha * M
            aJ.matrix_data() *= mAlpha ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf3_nok(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            // Scale stiffness matrix K and forcing vector f by h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // A = alpha * M
            aJ.matrix_data() *= mAlpha ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf4_nok(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            // Scale stiffness matrix K and forcing vector f h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 - beta3 * q3 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // A = alpha * M
            aJ.matrix_data() *= mAlpha ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::bdf5_nok(
            Matrix< real > & aJ,
            Vector< real > & aRHS  )
        {
            // Scale stiffness matrix K and forcing vector f by h
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            // beta-weighted history, single source of truth
            const Vector<real> & tSwap = this->collect_qhist();

            // b = h*f + M * ( beta0*q0  - beta1*q1 + beta2 * q2 - beta3 * q3 + beta4 * q4 )
            aRHS.vector_data() += aJ.matrix_data() * tSwap.vector_data();

            // A = alpha * M
            aJ.matrix_data() *= mAlpha ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Timestep::derivative( Matrix< real > & aJ,
                                  Vector< real > & aRHS  )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS -= mTimeStepMatrices->K().matrix_data() * mCalc->qold().vector_data() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::euler_no_stiffness( Matrix< real > & aJ,
                                          Vector< real > & aRHS )
        {
            aJ.matrix_data() += mTimeStepMatrices->M().matrix_data();

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;
            aRHS.vector_data() += aJ.matrix_data() * mCalc->qold().vector_data() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::derivative_no_stiffness( Matrix< real > & aJ,
                                                Vector< real > & aRHS )
        {
            // do nothing
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::static_subfield( Matrix< real > & aJ,
                                       Vector< real > & aRHS )
        {

            aRHS.vector_data() += mTimeStepMatrices->f()*mDeltaTime ;

            aJ.matrix_data() = mTimeStepMatrices->K().matrix_data() * mDeltaTime ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_bdf_coeffs_2()
        {
            // Variable-step BDF2 coefficients from Lagrange derivative formula
            // h = upcoming step, h1 = mH(0) = previous step
            real tH  = mDeltaTime ;
            real tH1 = mH( 0 ) ;

            // Cumulative sums: S2 = h + h1
            real tS2 = tH + tH1 ;

            // alpha = (2h + h1) / (h + h1)
            mAlpha = ( 2.0 * tH + tH1 ) / tS2 ;

            // beta0 = (h + h1) / h1  (coefficient of q_n)
            mBeta( 0 ) = tS2 / tH1 ;

            // beta1 = h^2 / (h1 * (h + h1))  (coefficient of q_{n-1})
            mBeta( 1 ) = tH * tH / ( tH1 * tS2 ) ;
        }
            
//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_bdf_coeffs_3()
        {
            // Variable-step BDF3 coefficients from Lagrange derivative formula
            // h = upcoming step, h1 = mH(0), h2 = mH(1)
            real tH  = mDeltaTime ;
            real tH1 = mH( 0 ) ;
            real tH2 = mH( 1 ) ;

            // Cumulative sums
            real tS2 = tH + tH1 ;
            real tS3 = tS2 + tH2 ;

            // alpha = 1 + h/S2 + h/S3
            mAlpha = 1.0 + tH / tS2 + tH / tS3 ;

            // beta0 = S2 * S3 / (h1 * (h1 + h2))
            mBeta( 0 ) = tS2 * tS3 / ( tH1 * ( tH1 + tH2 ) ) ;

            // beta1 = h^2 * S3 / (S2 * h1 * h2)
            mBeta( 1 ) = tH * tH * tS3 / ( tS2 * tH1 * tH2 ) ;

            // beta2 = h^2 * S2 / (S3 * (h1 + h2) * h2)
            mBeta( 2 ) = tH * tH * tS2 / ( tS3 * ( tH1 + tH2 ) * tH2 ) ;
        }
            
//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_bdf_coeffs_4()
        {
            // Variable-step BDF4 coefficients from Lagrange derivative formula
            // h = upcoming step, h1 = mH(0), h2 = mH(1), h3 = mH(2)
            real tH  = mDeltaTime ;
            real tH1 = mH( 0 ) ;
            real tH2 = mH( 1 ) ;
            real tH3 = mH( 2 ) ;

            // Cumulative sums
            real tS2 = tH + tH1 ;
            real tS3 = tS2 + tH2 ;
            real tS4 = tS3 + tH3 ;

            // alpha = 1 + h/S2 + h/S3 + h/S4
            mAlpha = 1.0 + tH / tS2 + tH / tS3 + tH / tS4 ;

            // beta0 = S2 * S3 * S4 / (h1 * (h1+h2) * (h1+h2+h3))
            mBeta( 0 ) = tS2 * tS3 * tS4
                / ( tH1 * ( tH1 + tH2 ) * ( tH1 + tH2 + tH3 ) ) ;

            // beta1 = h^2 * S3 * S4 / (S2 * h1 * h2 * (h2+h3))
            mBeta( 1 ) = tH * tH * tS3 * tS4
                / ( tS2 * tH1 * tH2 * ( tH2 + tH3 ) ) ;

            // beta2 = h^2 * S2 * S4 / (S3 * (h1+h2) * h2 * h3)
            mBeta( 2 ) = tH * tH * tS2 * tS4
                / ( tS3 * ( tH1 + tH2 ) * tH2 * tH3 ) ;

            // beta3 = h^2 * S2 * S3 / (S4 * (h1+h2+h3) * (h2+h3) * h3)
            mBeta( 3 ) = tH * tH * tS2 * tS3
                / ( tS4 * ( tH1 + tH2 + tH3 ) * ( tH2 + tH3 ) * tH3 ) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_bdf_coeffs_5()
        {
            // Variable-step BDF5 coefficients from Lagrange derivative formula
            // h = upcoming step, h1 = mH(0), h2 = mH(1), h3 = mH(2), h4 = mH(3)
            real tH  = mDeltaTime ;
            real tH1 = mH( 0 ) ;
            real tH2 = mH( 1 ) ;
            real tH3 = mH( 2 ) ;
            real tH4 = mH( 3 ) ;

            // Cumulative sums
            real tS2 = tH + tH1 ;
            real tS3 = tS2 + tH2 ;
            real tS4 = tS3 + tH3 ;
            real tS5 = tS4 + tH4 ;

            // alpha = 1 + h/S2 + h/S3 + h/S4 + h/S5
            mAlpha = 1.0 + tH / tS2 + tH / tS3 + tH / tS4 + tH / tS5 ;

            // beta0 = S2*S3*S4*S5 / (h1*(h1+h2)*(h1+h2+h3)*(h1+h2+h3+h4))
            mBeta( 0 ) = tS2 * tS3 * tS4 * tS5
                / ( tH1 * ( tH1 + tH2 ) * ( tH1 + tH2 + tH3 )
                    * ( tH1 + tH2 + tH3 + tH4 ) ) ;

            // beta1 = h^2*S3*S4*S5 / (S2*h1*h2*(h2+h3)*(h2+h3+h4))
            mBeta( 1 ) = tH * tH * tS3 * tS4 * tS5
                / ( tS2 * tH1 * tH2 * ( tH2 + tH3 ) * ( tH2 + tH3 + tH4 ) ) ;

            // beta2 = h^2*S2*S4*S5 / (S3*(h1+h2)*h2*h3*(h3+h4))
            mBeta( 2 ) = tH * tH * tS2 * tS4 * tS5
                / ( tS3 * ( tH1 + tH2 ) * tH2 * tH3 * ( tH3 + tH4 ) ) ;

            // beta3 = h^2*S2*S3*S5 / (S4*(h1+h2+h3)*(h2+h3)*h3*h4)
            mBeta( 3 ) = tH * tH * tS2 * tS3 * tS5
                / ( tS4 * ( tH1 + tH2 + tH3 ) * ( tH2 + tH3 ) * tH3 * tH4 ) ;

            // beta4 = h^2*S2*S3*S4 / (S5*(h1+h2+h3+h4)*(h2+h3+h4)*(h3+h4)*h4)
            mBeta( 4 ) = tH * tH * tS2 * tS3 * tS4
                / ( tS5 * ( tH1 + tH2 + tH3 + tH4 ) * ( tH2 + tH3 + tH4 )
                    * ( tH3 + tH4 ) * tH4 ) ;
        }

//------------------------------------------------------------------------------
    }
}
