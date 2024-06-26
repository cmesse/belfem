//
// Created by christian on 7/28/21.
//

#include "typedefs.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "fn_dot.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

            IWG_Timestep::IWG_Timestep(
                    const IwgType aType,
                    const IwgMode aMode,
                    const SymmetryMode aSymmetryMode,
                    const DofMode      aDofMode,
                    const SideSetDofLinkMode aSideSetDofLinkMode ) :
                    IWG( aType, aMode, aSymmetryMode, aDofMode, aSideSetDofLinkMode )
            {

            }

//------------------------------------------------------------------------------

            void
            IWG_Timestep::set_timestep( const real aDeltaTime )
            {
                BELFEM_ASSERT( mField != nullptr, "IWG was not linked to a field" );

                if( mField->parent()->master() == mField->rank() )
                {
                    mDeltaTime = aDeltaTime;

                    Vector< real > tDeltaTime(
                            mField->parent()->comm_table().length(),
                            aDeltaTime );

                    send( mField->parent()->comm_table(), tDeltaTime );
                }
                else
                {
                    receive( mField->parent()->master(), mDeltaTime );
                }
            }

//------------------------------------------------------------------------------

            real
            IWG_Timestep::compute_Ttheta( const uint aK )
            {
                return ( 1.0 - mTheta )
                             * dot( mGroup->n( aK ), mGroup->work_psi() )
                        + mTheta
                             * dot( mGroup->n( aK ), mGroup->work_phi() );
            }

//------------------------------------------------------------------------------

            real
            IWG_Timestep::compute_T0( const uint aK )
            {
                return dot( mGroup->n( aK ), mGroup->work_psi() );
            }

//------------------------------------------------------------------------------

            real
            IWG_Timestep::compute_T1( const uint aK )
            {
                return dot( mGroup->n( aK ), mGroup->work_phi() );
            }

//------------------------------------------------------------------------------

            void
            IWG_Timestep::set_euler_method( const real aTheta )
            {
                mTheta = aTheta ;
                if( aTheta == 0.0 )
                {
                    // set the interpolation function
                    mComputeT = & IWG_Timestep::compute_T0 ;
                }
                else if ( aTheta == 1.0 )
                {
                    mComputeT = & IWG_Timestep::compute_T1 ;
                }
                else
                {
                    mComputeT = & IWG_Timestep::compute_Ttheta ;
                }
            }

//------------------------------------------------------------------------------

            void
            IWG_Timestep::set_euler_method( const EulerMethod aEulerMethod )
            {
                switch( aEulerMethod )
                {
                    case( EulerMethod::ForwardExplicit ) :
                    {
                        // set the parameter
                        mTheta = 0.0 ;

                        // set the interpolation function
                        mComputeT
                            = & IWG_Timestep::compute_T0 ;

                        break ;
                    }
                    case( EulerMethod::CrankNicolson ) :
                    {
                        // set the parameter
                        mTheta = 0.5 ;

                        // set the interpolation function
                        mComputeT
                            = & IWG_Timestep::compute_Ttheta ;
                        break ;
                    }
                    case( EulerMethod::Galerkin ) :
                    {
                        // set the parameter
                        mTheta = 2.0 / 3.0 ;

                        // set the interpolation function
                        mComputeT
                            = & IWG_Timestep::compute_Ttheta ;
                        break ;
                    }
                    case( EulerMethod::BackwardImplicit ) :
                    {
                        // set the parameter
                        mTheta = 1.0 ;

                        // set the interpolation function
                        mComputeT
                            = & IWG_Timestep::compute_T1 ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Unknown EulerMethod" );
                    }
                }
            }

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */