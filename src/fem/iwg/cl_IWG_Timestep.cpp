    //
// Created by christian on 7/18/23.
//

#include "typedefs.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "fn_dot.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Calculator.hpp"


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
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::set_timestepping_method( const EulerMethod aMethod, const bool aHaveStiffness )
        {
            // remember the method
            mMethod = aMethod ;

            if( aHaveStiffness )
            {
                switch( aMethod )
                {
                    case( EulerMethod::Static ) :
                    {
                        mTimestep = & IWG_Timestep::static_subfield ;
                        break ;
                    }
                    case( EulerMethod::ForwardExplicit ) :
                    {
                        mTimestep = & IWG_Timestep::explicit_euler ;
                        break;
                    }
                    case( EulerMethod::CrankNicolson ) :
                    {
                        mTimestep = & IWG_Timestep::crank_nicolson ;
                        break ;
                    }
                    case( EulerMethod::Galerkin ) :
                    {
                        mTimestep = & IWG_Timestep::galerkin ;
                        break ;
                    }
                    case( EulerMethod::BackwardImplicit ) :
                    {
                        mTimestep = & IWG_Timestep::backwards_euler ;
                        break ;
                    }
                    case( EulerMethod::Derivative ) :
                    {
                        mTimestep = & IWG_Timestep::derivative ;
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
                    case( EulerMethod::ForwardExplicit ) :
                    case( EulerMethod::CrankNicolson ) :
                    case( EulerMethod::Galerkin ) :
                    case( EulerMethod::BackwardImplicit ) :
                    {
                        mTimestep = & IWG_Timestep::euler_no_stiffness ;
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
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_mkf( Element * aElement,
                     Matrix< real > & aM,
                     Matrix< real > & aK,
                     Vector< real > & aF,
                     Vector< real > & aQ )
        {
            BELFEM_ERROR( false, "compute_mkf not implemented for this IWG" );
        }

//------------------------------------------------------------------------------

        void
        IWG_Timestep::compute_jacobian_and_rhs( Element        * aElement,
                                                Matrix< real > & aJacobian,
                                                Vector< real > & aRHS )
        {

            this->compute_mkf(
                    aElement,
                    aJacobian,
                    mCalc->K(),
                    aRHS,
                    mCalc->q0() );


        }

//------------------------------------------------------------------------------
    }
}
