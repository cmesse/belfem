    //
// Created by christian on 7/18/23.
//

#include "typedefs.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "fn_dot.hpp"
#include "cl_IWG_Timestep.hpp"


namespace belfem
{
    namespace fem
    {

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
    }
}
