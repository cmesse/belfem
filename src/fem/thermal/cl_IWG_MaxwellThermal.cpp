//
// Created by christian on 10/27/25.
//
#include "cl_IWG_MaxwellThermal.hpp"

#include "matrices/mt_thermal_phi.hpp"
#include "matrices/mt_thermal_h.hpp"

namespace belfem
{
    namespace fem
    {
        IWG_MaxwellThermal::IWG_MaxwellThermal(
            ModelDimensionality  aModelDimensionality,
                   const IwgType aType ,
                   const IwgMode aMode ) :
            IWG_TransientHeatConduction( aModelDimensionality, aType, aMode )
        {
            this->initialize( IwgType::MaxwellThermal ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::M ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::K ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::F ) ;

            // Newton tangent blocks for T_h_newton
            mTimeStepMatrices->set_flag( MatrixFlag::dKdX_times_x ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::dMdX_times_x ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::dMdX_times_h ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::dFdX ) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_MaxwellThermal::compute_mkf(
            Element * aElement)
        {
            // link calculator with element
            mGroup->calculator()->link( aElement );

            mTimeStepMatrices->reset() ;

            // hot-path backstop; the always-active check sits in
            // link_to_group()
            BELFEM_ASSERT( mFunMKF != nullptr,
                "thermal kernel function not set for group %lu ( domain type %s )",
                ( long unsigned int ) mGroup->id(),
                to_string( mGroup->domain_type() ).c_str() );

            // call the function
            ( *mFunMKF )( mGroup->calculator(), mTimeStepMatrices );

        }

//------------------------------------------------------------------------------

        void
        IWG_MaxwellThermal::link_to_group( Group  * aGroup )
        {
            IWG::link_to_group(aGroup);

            mGroup = aGroup ;
            mCalc = aGroup->calculator() ;

            // kernel dispatch by domain type. This runs for empty groups
            // too, so the IWG-wide pointer never carries the previous
            // group's kernel; a domain type with no thermal kernel resets
            // it to nullptr and only errors when the group assembles
            switch ( aGroup->domain_type() )
            {
                case( DomainType::Conductor ) :
                case( DomainType::ThinShell ) :
                {
                    // collapsed kernel ( maxwell_kernel_collapse R11/R12 ):
                    // material math is dispatched inside calculator::MaxwellData,
                    // no material lookup at link time; the Newton variant adds
                    // the dcp/dT, dlambda/dT and drho/dT tangent blocks
                    mFunMKF = this->algorithm() == SolverAlgorithm::NewtonRaphson ?
                          & T_h_newton
                        : & T_h_picard ;
                    break ;
                }
                case( DomainType::Ferro ) :
                case( DomainType::Air ) :
                case( DomainType::Buffer ) :
                {
                    mFunMKF = & T_phi ;
                    break ;
                }
                default:
                {
                    mFunMKF = nullptr ;
                    break ;
                }
            }

            if( aGroup->number_of_elements() == 0 )
            {
                return;

            }

            // always-active guard: an assembled group without a kernel
            // must fail here, loudly, and not by jumping through the
            // pointer in compute_mkf()
            BELFEM_ERROR( mFunMKF != nullptr,
                "no thermal kernel function for group %lu ( domain type %s )",
                ( long unsigned int ) aGroup->id(),
                to_string( aGroup->domain_type() ).c_str() );

            uint tNumDofs = mGroup->parent() != nullptr ? mGroup->elements()(0)->number_of_local_dofs() : 0 ;
            mTimeStepMatrices->initialize( tNumDofs ) ;

        }

//------------------------------------------------------------------------------

        void
        IWG_MaxwellThermal::create_custom_vectors_and_matrices( Calculator * aCalc )
        {
            // n x 1 workspace for the conductivity tangent in T_h_newton:
            // Btg = trans( B ) * grad T, outer-multiplied with the shape row
            uint n = aCalc->group()->number_of_nodes_per_element() ;
            aCalc->create_matrix( "Btg", n, 1 );
        }
    }
}
