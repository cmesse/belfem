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

#include "cl_IWG_MaxwellPostproc.hpp"
#include "matrices/mt_maxwell_l2_phi.hpp"
#include "matrices/mt_maxwell_l2_h.hpp"
#include "matrices/mt_maxwell_l2_b.hpp"
#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    namespace fem
    {

        IWG_MaxwellPostproc::IWG_MaxwellPostproc(
            const maxwell::Formulation aFormulation,
            const ModelDimensionality aDimensionality,
                  const bool aHigherOrder ) :
            IWG( IwgType::Maxwell,
                         aDimensionality,
                         IwgMode::Direct, // might change to direct
                         SymmetryMode::Unsymmetric,
                         DofMode::AllBlocksEqual,
                         SideSetDofLinkMode::FacetOnly ),
            mFormulation( aFormulation )
        {
            switch( aFormulation )
            {
                case maxwell::Formulation::L2PhiH :
                {
                    mUseEdges = false ;
                    if( aDimensionality == ModelDimensionality::ThreeD )
                    {
                        mDofFields = { "Hx", "Hy", "Hz" } ;
                        mOtherFields = { "phi" };
                    }
                    else
                    {
                        mDofFields = { "Hx", "Hy" } ;
                        mOtherFields = { "phi", "Hz" };
                    }
                    break ;
                }
                case maxwell::Formulation::L2EdgeH :
                {
                    if( aDimensionality == ModelDimensionality::ThreeD )
                    {
                        mDofFields = { "Hx", "Hy", "Hz" } ;
                    }
                    else
                    {
                        mDofFields = { "Hx", "Hy" } ;
                        mOtherFields = { "Hz" };
                    }

                    mUseEdges = true ;
                    if( aHigherOrder == false )
                    {
                        mEdgeDofMultiplicity = 1.0 ;
                        mFaceDofMultiplicity = 0.0 ;
                        mOtherFields = { "edge_h" };
                    }
                    else
                    {
                        mEdgeDofMultiplicity = 2.0 ;
                        mFaceDofMultiplicity = 2.0 ;
                        mOtherFields = { "edge_h", "face_h" };
                    }
                    break ;
                }
                case maxwell::Formulation::L2PhiB :
                {
                    mUseEdges = false ;
                    if( aDimensionality == ModelDimensionality::ThreeD )
                    {
                        mDofFields = { "Bx", "By", "Bz" } ;
                        mOtherFields = { "phi" };
                    }
                    else
                    {
                        mDofFields = { "Bx", "By" } ;
                        mOtherFields = { "phi", "Bz" };
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "invalid formulation");
                }
            }

            this->initialize( IwgType::Maxwell );
        }

        void
        IWG_MaxwellPostproc::create_custom_vectors_and_matrices( Calculator * aCalc )
        {
            uint tNumNodes = mesh::number_of_nodes( aCalc->group()->element_type() ) ;
            uint tNumDim   = mesh::dimension( aCalc->group()->element_type() ) ;

            aCalc->create_matrix("L2N", tNumDim, tNumNodes * tNumDim  );
            aCalc->create_vector( "nedelec_h", aCalc->num_nedelec_dofs() );
            aCalc->create_vector( "nedelec_a", aCalc->num_nedelec_dofs() );
            aCalc->create_vector( "a_z", tNumNodes, EntityType::NODE );
            aCalc->create_matrix( "HelpA", 2, tNumNodes );
            aCalc->create_vector( "HelpB", 2 );
            aCalc->create_vector( "h", tNumDim );

            switch( mFormulation  )
            {
                case maxwell::Formulation::L2PhiH :
                {
                    mFunKF = & maxwell::l2_phi ;
                    break ;
                }
                case maxwell::Formulation::L2PhiB :
                {
                    mFunKF = & maxwell::l2_phi_ferro ;
                    break ;
                }
                case maxwell::Formulation::L2EdgeH :
                {
                    mFunKF = & maxwell::l2_h ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid formulation");
                }
            }
        }

        void
        IWG_MaxwellPostproc::compute_jacobian_and_rhs(
            Element      * aElement,
            Matrix<real> & aJacobian,
            Vector<real> & aRHS)
        {
            mCalc->link( aElement );
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            ( *mFunKF )( mCalc, aJacobian, aRHS ) ;
        }



    }
}