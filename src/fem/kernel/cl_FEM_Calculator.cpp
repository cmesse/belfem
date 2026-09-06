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

#include <iostream>


#include "assert.hpp"
#include "commtools.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Calculator.hpp"
#include "meshtools.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_entity_type.hpp"
#include "fn_norm.hpp"

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "nedelec/cl_EF_PENTA6TS.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "cl_EdgeFunctionFactory.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Controller.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_SideSet.hpp"
#include "fn_max.hpp"

namespace belfem
{
    namespace fem
    {
        namespace calculator
        {
//------------------------------------------------------------------------------

            VectorData::VectorData( const string & aLabel,
                                    uint aSize,
                                    const EntityType aType )  :
                                    mLabel( aLabel ),
                                    mType( aType )
            {
                mVectorData.set_size( aSize, BELFEM_QUIET_NAN );
            }

//------------------------------------------------------------------------------

            MatrixData::MatrixData( const string & aLabel,
                                    const uint aNumRows,
                                    const uint aNumCols ) :
                                    mLabel( aLabel )
            {
                mMatrixData.set_size( aNumRows, aNumCols, BELFEM_QUIET_NAN );
            }

//------------------------------------------------------------------------------

            MaxwellData::MaxwellData( Calculator * aCalculator,
                                      Kernel * aMaxwellKernel,
                                      Kernel * aThermalKernel ) :
                    mMaxwellCalculator( aMaxwellKernel->dofmgr()->block( aCalculator->group()->id() )->calculator() ),
                    // the block_exists() check matters: block() returns the
                    // EmptyBlock ( with a live calculator ) for unknown ids,
                    // which must never be adopted as thermal peer ( side
                    // connector blocks are not part of the thermal kernel )
                    mThermalCalculator( aThermalKernel != nullptr
                        && aThermalKernel->dofmgr()->block_exists( aCalculator->group()->id() )
                        ? aThermalKernel->dofmgr()->block( aCalculator->group()->id() )->calculator() : nullptr ),
                    mMaterial( aMaxwellKernel->dofmgr()->block( aCalculator->group()->id() )->material() ),
                    mTime( aMaxwellKernel->controller()->time() ),
                    mH(  this->link_vector( "h" ) ),
                    mHn( this->link_vector( "hn" ) ),
                    mHt( this->link_vector( "ht" ) ),
                    mB(  this->link_vector( "b" ) ),
                    mBn( this->link_vector( "bn" ) ),
                    mBt( this->link_vector( "bt" ) ),
                    mJ(  this->link_vector( "j" ) ),
                    mN(  this->link_vector( "normal" ) )
            {
                mLastIndex.set_size( static_cast< uint >( MaxwellDataValue::UNDEFINED), BELFEM_UINT_MAX );

                bool tIsThinShell = aCalculator->group()->domain_type() == DomainType::ThinShell ;

                bool tIsNedelec    = aCalculator->group()->domain_type() == DomainType::Conductor || tIsThinShell ;

                BELFEM_ERROR( mMaterial != nullptr, "No material was assigned for block %lu",
                    ( long unsigned int ) aCalculator->group()->id() );

                bool tIsConstantMu = mMaterial->is_constant( MaterialProperty::mu );

                // short-circuit: constant_property() asserts is_constant(), which
                // fails for any material with a b-h curve ( mu slot holds NaN )
                bool tIsConstantMu0 = tIsConstantMu
                        && mMaterial->constant_property( MaterialProperty::mu ) == constant::mu0 ;

                // flag if we need jc
                bool tIsHTS = mMaterial->have( MaterialProperty::jc ) ;

                // flag if rho depends on the field magnitude and the b-j angle;
                // interrogate the property rather than the material type, so
                // lookup and user-defined materials keep their field dependence
                bool tIsMetal = ! tIsHTS && mMaterial->depends( MaterialProperty::rho, MaterialDependency::normB ) ;

                // lambda field dependence is routed independently of rho: an HTS
                // taper carries a field-independent lambda(T) while a normal
                // conductor may carry a field-dependent lambda(T,normB,beta).
                // Follow the property, matching legacy mt_thermal_h.cpp:81-87.
                bool tLambdaFieldDependent =
                        mMaterial->depends( MaterialProperty::lambda, MaterialDependency::normB ) ;

                // beta conventions are per material family and must never mix:
                // ( b, n ) bn_angle for HTS / REBCO, whose lambda is T-only,
                // ( b, j ) bj_angle for normal metals. A material combining jc
                // with a field-dependent lambda would alias the two conventions
                // in the shared beta slot — reject it loudly at setup
                BELFEM_ERROR( ! ( tIsHTS && tLambdaFieldDependent ),
                    "material %s combines jc with a field-dependent lambda - beta conventions would alias",
                    mMaterial->label().c_str() );

                // upper edge of the physical temperature window for compute_T:
                // transient iterates are clamped into [ gTmin, mTmax ]
                if ( mMaterial->is_constant( MaterialProperty::T_max ) )
                {
                    mTmax = mMaterial->constant_property( MaterialProperty::T_max );
                }

                // density at the undeformed-mesh reference temperature ( see
                // the member note: density( T ) would double-count expansion ).
                //
                // The correction multiplies BOTH branches. It used to sit on
                // the ref_density branch alone, which made it dead for every
                // alloy: Alloy::create_splines calls set_custom( density ),
                // and set_custom raises the have-flag, so an alloy always
                // takes the first branch even though it also defines
                // ref_density. The key was silently ignored by exactly the
                // materials it exists for.
                if ( mMaterial->have( MaterialProperty::density ) )
                {
                    mDensity = mMaterial->density( gTroom );
                }
                else if ( mMaterial->have( MaterialProperty::ref_density ) )
                {
                    mDensity = mMaterial->ref_density();
                }

                // scales the density that enters the thermal mass matrix
                // ( mt_thermal_h: M += N' * density * cp * N * dV ) so a
                // meshed volume can carry the mass of the thinner layer it
                // really contains. Defaults to 1.0, so untouched decks are
                // unaffected. Deliberately NOT applied to any conductivity.
                mDensity *= mMaterial->constant_property(
                        MaterialProperty::density_correction );

                if ( aCalculator->mesh()->number_of_dimensions() == 2 )
                {
                    mFunX = & MaxwellData::compute_x_2d ;
                }
                else
                {
                    mFunX = & MaxwellData::compute_x_3d ;
                }
                if ( mThermalCalculator != nullptr )
                {
                    mFunT = & MaxwellData::compute_T_fem ;
                }
                else
                {
                    mFunT = & MaxwellData::compute_T_const ;
                }

                // field derivatives of rho ( the B/beta Newton tangent
                // channel ): default to zero; the metal branches bind both,
                // the HTS branches bind the |B| channel ( jc/n lookup tables )
                // and deliberately leave beta at zero ( see below )
                mFundRhodB    = & MaxwellData::return_zero ;
                mFundRhodBeta = & MaxwellData::return_zero ;

                // artificial volumetric heat load ( material heating plugin ):
                // bound once, so the thermal kernels pay one indirect call and
                // no branch per integration point when there is none
                mFunHeat = mMaterial->have_heating()
                        ? & MaxwellData::compute_heatload_user
                        : & MaxwellData::return_zero ;

                bool tIsSideConnector =
                       aCalculator->group()->domain_type() == DomainType::LeftCoating
                    || aCalculator->group()->domain_type() == DomainType::RightCoating ;

                if ( tIsSideConnector )
                {
                    // edge-coating wall: recovery element on the tape slit
                    // edge. Field assembly ( h = ht + hb + hn via the master
                    // layer element, seam-node T ) is connector-specific;
                    // the rho family is the plain metal one.

                    // belt + braces: assign_materials already gates on this
                    BELFEM_ERROR( mMaterial->type() == MaterialType::PureMetal,
                        "side connector block %lu: material %s must be a pure metal",
                        ( long unsigned int ) aCalculator->group()->id(),
                        mMaterial->label().c_str() );

                    // the wall kernel has no dMdx tangent blocks, so a
                    // field-dependent mu is not supported here ( decision
                    // 2026-08-09: plated walls are copper, mu = mu0 )
                    BELFEM_ERROR( tIsConstantMu0 || tIsConstantMu,
                        "side connector block %lu: material %s must have a constant mu",
                        ( long unsigned int ) aCalculator->group()->id(),
                        mMaterial->label().c_str() );

                    mFunH = & MaxwellData::compute_h_side_connector ;
                    mFunB = & MaxwellData::compute_b_bulk ;
                    mFunT = & MaxwellData::compute_T_side_connector ;

                    mFunMu = tIsConstantMu0 ?
                          & MaxwellData::compute_mu_0
                        : & MaxwellData::compute_mu_const ;
                    mFundMudH = & MaxwellData::compute_dmu_zero ;

                    // metal rho family, verbatim reuse
                    mFunRho       = & MaxwellData::compute_rho_metal ;
                    mFundRhodT    = & MaxwellData::compute_drhodT_metal ;
                    mFundRhodJ    = & MaxwellData::return_zero ;
                    mFundRhodB    = & MaxwellData::compute_drhodb_metal ;
                    mFundRhodBeta = & MaxwellData::compute_drhodbeta_metal ;

                    // wall temperature source: seam nodes if a thermal kernel
                    // exists, gTbulk otherwise ( see compute_T_side_connector )
                    mHaveSeamT = ( aThermalKernel != nullptr );

                    // scratch for the master in-plane field
                    mWork.set_size( 3, BELFEM_QUIET_NAN );
                }
                else if ( tIsThinShell )
                {
                    // no edge/node test: ThinShell BLOCKS are routed through
                    // the Conductor dof table ( edge_h ) by
                    // FieldList::collect_block_dofs, so a shell block is always
                    // edge-interpolated. Not the FieldList::ThinShell member:
                    // that is the SIDESET table and it does carry phi
                    mFunH = & MaxwellData::compute_h_ts_edge ;

                    if ( tIsConstantMu0 )
                    {
                        mFunMu = & MaxwellData::compute_mu_0 ;
                        mFundMudH = & MaxwellData::compute_dmu_zero ;
                    }
                    else if ( tIsConstantMu )
                    {
                        mFunMu = & MaxwellData::compute_mu_const ;
                        mFundMudH = & MaxwellData::compute_dmu_zero ;
                    }
                    else
                    {
                        mFunMu = & MaxwellData::compute_mu_h ;
                        mFundMudH = & MaxwellData::compute_dmu_material ;
                    }

                    mFunB = & MaxwellData::compute_b_ts ;


                    if ( tIsMetal )
                    {
                        mFunRho       = & MaxwellData::compute_rho_metal ;
                        mFundRhodT    = & MaxwellData::compute_drhodT_metal ;
                        mFundRhodJ    = & MaxwellData::return_zero ;
                        mFundRhodB    = & MaxwellData::compute_drhodb_metal ;
                        mFundRhodBeta = & MaxwellData::compute_drhodbeta_metal ;
                    }
                    else if ( tIsHTS )
                    {
                        // jc(T,|B|,θ) / n(T,|B|,θ) from a lookup table
                        // make rho field-dependent, so the |B| channel of the
                        // Newton tangent must be bound — for constant jc/n
                        // the derivative evaluates to exactly zero and the
                        // kernel early-outs, so those decks are unchanged.
                        // The T channel ( T-leg, 2026-08-13 ) is bound
                        // per law as well: it feeds T_h_newton's
                        // quench-feedback block and, unlike |B|, keeps a
                        // nonzero rho_n term for constant-jc materials.
                        // mFundRhodBeta deliberately STAYS return_zero: the
                        // HTS angle is bn_angle ( field to tape normal ),
                        // while add_rho_field_tangent differentiates
                        // bj_angle ( field to current, metal Kohler ) —
                        // binding it would apply the wrong ∂β/∂q rows
                        // ( 2026-08-13 audit, both voices )
                        if ( mMaterial->have_defect() )
                        {
                            if ( mMaterial->use_piecewise() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_piecewise_ts_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_piecewise_ts_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_piecewise_ts_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_piecewise_ts_defect ;
                            }
                            else if ( mMaterial->use_riva() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_riva_ts_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_riva_ts_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_riva_ts_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_riva_ts_defect ;
                            }
                            else
                            {
                                mFunRho    = & MaxwellData::compute_rho_powerlaw_ts_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_powerlaw_ts_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_powerlaw_ts_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_powerlaw_ts_defect ;
                            }
                        }
                        else
                        {
                            if ( mMaterial->use_piecewise() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_piecewise_ts ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_piecewise_ts ;
                                mFundRhodB = & MaxwellData::compute_drhodb_piecewise_ts ;
                                mFundRhodT = & MaxwellData::compute_drhodT_piecewise_ts ;
                            }
                            else if ( mMaterial->use_riva() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_riva_ts ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_riva_ts ;
                                mFundRhodB = & MaxwellData::compute_drhodb_riva_ts ;
                                mFundRhodT = & MaxwellData::compute_drhodT_riva_ts ;
                            }
                            else
                            {
                                mFunRho    = & MaxwellData::compute_rho_powerlaw_ts ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_powerlaw_ts ;
                                mFundRhodB = & MaxwellData::compute_drhodb_powerlaw_ts ;
                                mFundRhodT = & MaxwellData::compute_drhodT_powerlaw_ts ;
                            }
                        }
                    }
                    else
                    {
                        mFunRho       = & MaxwellData::compute_rho_bulk ;
                        mFundRhodJ    = & MaxwellData::return_zero ;
                        mFundRhodT    = & MaxwellData::compute_drhodT_bulk ;
                    }
                }
                else // bulk
                {
                    // Buffer deliberately takes the node path too: buffer
                    // blocks carry nodal phi ( Air dof table ), so
                    // h = -grad phi; a thin-shell h would run facet
                    // machinery that volume blocks do not have ( F4 )
                    if ( tIsNedelec )
                    {
                        mFunH = & MaxwellData::compute_h_bulk_edge ;
                    }
                    else
                    {
                        mFunH = & MaxwellData::compute_h_bulk_node ;
                    }

                    if ( tIsConstantMu0 )
                    {
                        mFunMu = & MaxwellData::compute_mu_0 ;
                        mFundMudH = & MaxwellData::compute_dmu_zero ;
                    }
                    else if ( tIsConstantMu )
                    {
                        mFunMu = & MaxwellData::compute_mu_const ;
                        mFundMudH = & MaxwellData::compute_dmu_zero ;
                    }
                    else
                    {
                        mFunMu = & MaxwellData::compute_mu_h ;
                        mFundMudH = & MaxwellData::compute_dmu_material ;
                    }


                    mFunB = & MaxwellData::compute_b_bulk ;

                    // Every branch must bind mFundRhodT, exactly as the thin
                    // shell branch above does. It is the thermal Newton tangent
                    // channel ( T_h_newton needs drho/dT for the Joule term ),
                    // so a magnetic-only run never touches it and a missing
                    // assignment stays invisible until the first coupled Newton
                    // iteration calls through a null member pointer.
                    if ( tIsMetal )
                    {
                        mFunRho       = & MaxwellData::compute_rho_metal ;
                        mFundRhodT    = & MaxwellData::compute_drhodT_metal ;
                        mFundRhodJ    = & MaxwellData::return_zero ;
                        mFundRhodB    = & MaxwellData::compute_drhodb_metal ;
                        mFundRhodBeta = & MaxwellData::compute_drhodbeta_metal ;
                    }
                    else if ( tIsHTS )
                    {
                        // |B| + T channels, cf. the thin-shell branch
                        // above. The bulk beta is beta_dummy(), so the beta
                        // channel stays return_zero here for that reason as
                        // well
                        if ( mMaterial->have_defect() )
                        {
                            if ( mMaterial->use_piecewise() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_piecewise_bulk_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_piecewise_bulk_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_piecewise_bulk_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_piecewise_bulk_defect ;
                            }
                            else if ( mMaterial->use_riva() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_riva_bulk_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_riva_bulk_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_riva_bulk_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_riva_bulk_defect ;
                            }
                            else
                            {
                                mFunRho    = & MaxwellData::compute_rho_powerlaw_bulk_defect ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_powerlaw_bulk_defect ;
                                mFundRhodB = & MaxwellData::compute_drhodb_powerlaw_bulk_defect ;
                                mFundRhodT = & MaxwellData::compute_drhodT_powerlaw_bulk_defect ;
                            }
                        }
                        else
                        {
                            if ( mMaterial->use_piecewise() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_piecewise_bulk ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_piecewise_bulk ;
                                mFundRhodB = & MaxwellData::compute_drhodb_piecewise_bulk ;
                                mFundRhodT = & MaxwellData::compute_drhodT_piecewise_bulk ;
                            }
                            else if ( mMaterial->use_riva() )
                            {
                                mFunRho    = & MaxwellData::compute_rho_riva_bulk ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_riva_bulk ;
                                mFundRhodB = & MaxwellData::compute_drhodb_riva_bulk ;
                                mFundRhodT = & MaxwellData::compute_drhodT_riva_bulk ;
                            }
                            else
                            {
                                mFunRho    = & MaxwellData::compute_rho_powerlaw_bulk ;
                                mFundRhodJ = & MaxwellData::compute_drhodj_powerlaw_bulk ;
                                mFundRhodB = & MaxwellData::compute_drhodb_powerlaw_bulk ;
                                mFundRhodT = & MaxwellData::compute_drhodT_powerlaw_bulk ;
                            }
                        }
                    }
                    else
                    {
                        mFunRho       = & MaxwellData::compute_rho_bulk ;
                        mFundRhodJ    = & MaxwellData::return_zero ;
                        mFundRhodT    = & MaxwellData::compute_drhodT_bulk ;
                    }
                }

                // Nothing below may leave a rho-family pointer null: the thermal
                // Newton dereferences mFundRhodT on every element of every
                // conducting block, and a null there is a jump to address 0
                // rather than a diagnosable error.
                BELFEM_ASSERT( mFunRho != nullptr,
                    "MaxwellData: mFunRho was not assigned for material %s",
                    mMaterial->label().c_str() );
                BELFEM_ASSERT( mFundRhodT != nullptr,
                    "MaxwellData: mFundRhodT was not assigned for material %s",
                    mMaterial->label().c_str() );
                BELFEM_ASSERT( mFundRhodJ != nullptr,
                    "MaxwellData: mFundRhodJ was not assigned for material %s",
                    mMaterial->label().c_str() );

                // lambda dispatch is decoupled from the rho branch above so that
                // every material (metal, HTS, bulk) gets a valid lambda pointer.
                if ( tLambdaFieldDependent )
                {
                    mFunLambda    = & MaxwellData::compute_lambda_metal ;
                    mFundLambdadT = & MaxwellData::compute_dlambdadT_metal ;
                }
                else
                {
                    mFunLambda    = & MaxwellData::compute_lambda_bulk ;
                    mFundLambdadT = & MaxwellData::compute_dlambdadT_bulk ;
                }
            }

//------------------------------------------------------------------------------

            void
            MaxwellData::prepare_side_connector_frame()
            {
                Element * tElement = mMaxwellCalculator->element();
                Element * tMaster  = tElement->reference();

                BELFEM_ASSERT( tElement->element()->type() == ElementType::HEX8TB,
                    "expect a HEX8TB wall element" );

                // master layer-block calculator, linked to this wall's master
                mReferenceCalc = mMaxwellCalculator->group()->parent()->block(
                    tMaster->element()->block_id() )->calculator();
                mReferenceCalc->link( tMaster );

                // normal field component, evaluated on the master ( writes the
                // master's "hn" and "n" workspaces; linear elements: constant
                // per element, so k = 0 ). Both results are copied over: the
                // wall's own "hn" enters h in compute_h_side_connector, the
                // normal seeds the frame below
                mHn = compute_hn( mReferenceCalc, 0 );

                Vector< real > & tNormal   = mMaxwellCalculator->vector( "normal" );
                Vector< real > & tTangent  = mMaxwellCalculator->vector( "tangent" );
                Vector< real > & tBinomial = mMaxwellCalculator->vector( "binomial" );

                tNormal = mReferenceCalc->vector( "normal" );

                // frame from the wall corner nodes: tangent along the curve
                // direction, binomial = n x t. The two vectors hold corner
                // COORDINATES until the normalization below
                switch ( mMaxwellCalculator->group()->domain_type() )
                {
                    case DomainType::LeftCoating :
                    {
                        tElement->element()->node( 3 )->get_coords( tBinomial );
                        tElement->element()->node( 2 )->get_coords( tTangent );
                        break ;
                    }
                    case DomainType::RightCoating :
                    {
                        tElement->element()->node( 0 )->get_coords( tBinomial );
                        tElement->element()->node( 1 )->get_coords( tTangent );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "invalid domain type on block %lu",
                            ( long unsigned int ) mMaxwellCalculator->group()->id() );
                    }
                }
                tTangent -= tBinomial ;
                tTangent /= norm( tTangent );
                tBinomial = cross( tNormal, tTangent );

                mBinomialVec = & tBinomial ;
                mTseamVec    = & mMaxwellCalculator->vector( "Tseam" );

                // seam temperatures from the mesh field ( thermal runs only );
                // the wall itself is not part of the thermal problem. Reads
                // are ORIGINAL-NORMALIZED and STATION-ORDERED: the wall nodes
                // are decoupled duplicates, so the tape value lives at
                // node->original(), and the bilinear interpolation in
                // compute_T_side_connector expects the wall station order —
                // the recovery facet's nodes are in master canonical face
                // order and must not be indexed by station. The visualization
                // copy onto the wall nodes' own field entries lives in the
                // MaxwellPostprocessor: assembly is rank-local and the
                // save-time field gather only collects dof-carrying entities,
                // so writes from here would be lost on parallel runs
                if ( mHaveSeamT )
                {
                    const Vector< real > & T = mMaxwellCalculator->group()
                        ->parent()->mesh()->field_data( "T" );

                    // station pairing of the lateral faces: ( 3|0, 2|1, 6|5,
                    // 7|4 ); the tape side is the ( 0,1,5,4 ) face for the
                    // left connector and the ( 3,2,6,7 ) face for the right
                    const uint tLeftFace[ 4 ]  = { 3, 2, 6, 7 };
                    const uint tRightFace[ 4 ] = { 0, 1, 5, 4 };

                    const uint * tTapeFace =
                        mMaxwellCalculator->group()->domain_type()
                            == DomainType::LeftCoating ?
                        tRightFace : tLeftFace ;

                    Vector< real > & tT = * mTseamVec ;

                    for ( uint k=0; k<4; ++k )
                    {
                        tT( k ) = T( tElement->element()->node(
                            tTapeFace[ k ] )->original()->index() );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            MaxwellData::side_connector_xi_eta(
                const real aPsi, real & aXi, real & aEta ) const
            {
                uint s = mMaxwellCalculator->element()->facet()->index_on_master() ;

                switch ( mMaxwellCalculator->group()->domain_type() )
                {
                    case DomainType::LeftCoating :
                    {
                        switch ( s )
                        {
                            case 0 :
                            {
                                aXi = 0.5 * ( aPsi + 1. );
                                aEta = 1. - aXi ;
                                return ;
                            }
                            case 1 :
                            {
                                aXi = 0. ;
                                aEta = 0.5 * ( 1. + aPsi );
                                return ;
                            }
                            case 2 :
                            {
                                aXi = 0.5 * ( 1. + aPsi );
                                aEta = 0. ;
                                return ;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "invalid side index" );
                            }
                        }
                    }
                    case DomainType::RightCoating :
                    {
                        switch ( s )
                        {
                            case 0 :
                            {
                                aEta = 0.5 * ( aPsi + 1. );
                                aXi = 1. - aEta ;
                                return ;
                            }
                            case 1 :
                            {
                                aXi = 0. ;
                                aEta = 0.5 * ( 1. - aPsi );
                                return ;
                            }
                            case 2 :
                            {
                                aXi = 0.5 * ( 1. - aPsi );
                                aEta = 0. ;
                                return ;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "invalid side index" );
                            }
                        }
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "invalid domain type" );
                    }
                }
            }

//------------------------------------------------------------------------------

            const Vector< real > &
            MaxwellData::compute_h_side_connector( const uint aIndex )
            {
                if ( ! mFrameCurrent )
                {
                    this->prepare_side_connector_frame();
                    mFrameCurrent = true ;
                }

                // tangential part from the wall's own edge dofs
                mHt = mMaxwellCalculator->E( aIndex ) * mMaxwellCalculator->q() ;

                // in-plane field of the master tape element at the mapped point
                const Matrix< real > & tPoints =
                    mMaxwellCalculator->integration()->points() ;

                real xi ;
                real eta ;
                this->side_connector_xi_eta( tPoints( 0, aIndex ), xi, eta );

                const Matrix< real > & Em = maxwell::side_connector_edge_function(
                    mReferenceCalc, xi, eta, tPoints( 2, aIndex ) );

                mWork = Em * mReferenceCalc->q() ;

                // h = hb + ht + hn
                const Vector< real > & tBinomial = * mBinomialVec ;
                mH  = dot( mWork, tBinomial ) * tBinomial ;
                mH += mHt ;
                mH += mHn ;

                return mH ;
            }

//------------------------------------------------------------------------------

            real
            MaxwellData::compute_T_side_connector( const uint aIndex )
            {
                if ( ! this->is_current( MaxwellDataValue::T, aIndex ) )
                {
                    real T ;

                    if ( ! mHaveSeamT )
                    {
                        // magnetic-only run: bulk temperature
                        T = gTbulk ;
                    }
                    else
                    {
                        if ( ! mFrameCurrent )
                        {
                            this->prepare_side_connector_frame();
                            mFrameCurrent = true ;
                        }

                        // bilinear interpolation between the seam supports
                        const Matrix< real > & tPoints =
                            mMaxwellCalculator->integration()->points() ;
                        const Vector< real > & tT = * mTseamVec ;

                        real xi   = tPoints( 0, aIndex );
                        real zeta = tPoints( 2, aIndex );

                        real a = 1. - xi ;
                        real b = 1. + xi ;

                        T = 0.25 * ( ( 1. - zeta ) * ( a * tT( 0 ) + b * tT( 1 ) )
                                   + ( 1. + zeta ) * ( b * tT( 2 ) + a * tT( 3 ) ) );
                    }

                    // same contract as compute_T_fem: clamp into the physical
                    // window, zero dT-derivatives while clamped
                    mTClamped = ( T < gTmin ) || ( T > mTmax );
                    mT = mTClamped ? std::clamp( T, gTmin, mTmax ) : T ;

                    this->set( MaxwellDataValue::T, aIndex );
                }
                return mT ;
            }

        }

//------------------------------------------------------------------------------

        namespace maxwell
        {
            // defined here rather than in mt_maxwell_h.cpp: the function is
            // the friend of EF_PENTA6TS that MaxwellData's side connector
            // path evaluates ( see cl_EF_PENTA6TS.hpp )
            const Matrix< real > &
            side_connector_edge_function( Calculator * aCalc, const real xi, const real eta, const real zeta )
            {
                // this is a copy of the precompute function
                real g[ 3 ];
                real h[ 3 ];
                real f[ 2 ];

                g[ 0 ] = xi ;
                g[ 1 ] = eta ;
                g[ 2 ] = 1.-xi-eta ;

                h[ 0 ] = eta ;
                h[ 1 ] = g[ 2 ];
                h[ 2 ] = xi ;

                f[ 0 ] = 0.5 * ( 1.-zeta );
                f[ 1 ] = 0.5 * ( 1.+zeta );

                BELFEM_ASSERT( aCalc->group()->element_type() == ElementType::PENTA6TS, "Expect PENTA6TS element type" );
                auto * tEdgeFunction = reinterpret_cast< EF_PENTA6TS * >( aCalc->edge_function() );

                // have been computed through linking with element
                const real * NablaXi   = tEdgeFunction->mNablaXi ;
                const real * NablaEta  = tEdgeFunction->mNablaEta ;
                const real * NablaZeta = tEdgeFunction->mNablaZeta ;
                const real * s         = tEdgeFunction->mS ;
                Matrix< real > & E     = tEdgeFunction->mE ;

                E( 0, 0 ) = g[ 0 ] * NablaEta[ 0 ] - h[ 0 ] * NablaXi[ 0 ];
                E( 1, 0 ) = g[ 0 ] * NablaEta[ 1 ] - h[ 0 ] * NablaXi[ 1 ];
                E( 2, 0 ) = g[ 0 ] * NablaEta[ 2 ] - h[ 0 ] * NablaXi[ 2 ];

                E( 0, 1 ) = g[ 1 ] * NablaZeta[ 0 ] - h[ 1 ] * NablaEta[ 0 ];
                E( 1, 1 ) = g[ 1 ] * NablaZeta[ 1 ] - h[ 1 ] * NablaEta[ 1 ];
                E( 2, 1 ) = g[ 1 ] * NablaZeta[ 2 ] - h[ 1 ] * NablaEta[ 2 ];

                E( 0, 2 ) = g[ 2 ] * NablaXi[ 0 ] - h[ 2 ] * NablaZeta[ 0 ];
                E( 1, 2 ) = g[ 2 ] * NablaXi[ 1 ] - h[ 2 ] * NablaZeta[ 1 ];
                E( 2, 2 ) = g[ 2 ] * NablaXi[ 2 ] - h[ 2 ] * NablaZeta[ 2 ];

                // duplicate and finalize upper face
                E( 0, 3 ) = s[ 3 ] * E( 0, 0 ) * f[ 1 ];
                E( 1, 3 ) = s[ 3 ] * E( 1, 0 ) * f[ 1 ];
                E( 2, 3 ) = s[ 3 ] * E( 2, 0 ) * f[ 1 ];

                E( 0, 4 ) = s[ 4 ] * E( 0, 1 ) * f[ 1 ];
                E( 1, 4 ) = s[ 4 ] * E( 1, 1 ) * f[ 1 ];
                E( 2, 4 ) = s[ 4 ] * E( 2, 1 ) * f[ 1 ];

                E( 0, 5 ) = s[ 5 ] * E( 0, 2 ) * f[ 1 ];
                E( 1, 5 ) = s[ 5 ] * E( 1, 2 ) * f[ 1 ];
                E( 2, 5 ) = s[ 5 ] * E( 2, 2 ) * f[ 1 ];

                E( 0, 0 ) *= s[ 0 ] * f[ 0 ];
                E( 1, 0 ) *= s[ 0 ] * f[ 0 ];
                E( 2, 0 ) *= s[ 0 ] * f[ 0 ];

                E( 0, 1 ) *= s[ 1 ] * f[ 0 ];
                E( 1, 1 ) *= s[ 1 ] * f[ 0 ];
                E( 2, 1 ) *= s[ 1 ] * f[ 0 ];

                E( 0, 2 ) *= s[ 2 ] * f[ 0 ];
                E( 1, 2 ) *= s[ 2 ] * f[ 0 ];
                E( 2, 2 ) *= s[ 2 ] * f[ 0 ];

                return E ;
            }
        }



//------------------------------------------------------------------------------

        Calculator::Calculator( Group * aGroup, const ModelDimensionality aDimensionality ) :
                mGroup( aGroup ),
                mMesh( aGroup->parent()->mesh() ),
                mDimensionality( aDimensionality ),
                mTimestep( aGroup->parent()->iwg()->delta_time() )
        {
            mFunBJAngle = aDimensionality == ModelDimensionality::ThreeD ?
                & Calculator::bj_angle_3d : &Calculator::bj_angle_2d;
        }

//------------------------------------------------------------------------------

        Calculator::Calculator( Group * aGroup, Mesh * aMesh ) :
                mGroup( aGroup ),
                mMesh( aMesh ),
                mDimensionality( ModelDimensionality::ThreeD ),
                mTimestep( aMesh->time_stamp() )
        {
            mFunBJAngle = & Calculator::bj_angle_3d ;
        }

//------------------------------------------------------------------------------

        Calculator::~Calculator()
        {
            if ( mMaxwellData != nullptr )
            {
                delete mMaxwellData ;
            }

            // delete calculator matrices
            for( calculator::MatrixData * tMatrix : mMatrices )
            {
                delete tMatrix ;
            }

            // delete calculator vectors
            for( calculator::VectorData * tVector : mVectors )
            {
                delete tVector ;
            }

            if( mEdgeFunction != nullptr )
            {
                delete mEdgeFunction ;
            }
            for( EdgeFunction * tEdgeFuction : mEdgeFunctionsMaster )
            {
                delete tEdgeFuction ;
            }
            for( EdgeFunction * tEdgeFuction : mEdgeFunctionsSlave )
            {
                delete tEdgeFuction ;
            }
            if( mDomainIntegration != nullptr )
            {
                delete mDomainIntegration ;
            }
            if ( mLinearIntegration != nullptr )
            {
                delete mLinearIntegration ;
            }
        }

//------------------------------------------------------------------------------


        void
        Calculator::initialize_integration(
                const ElementType       aElementType,
                const InterpolationType aInterpolationType )
        {
            if( mDomainIntegration != nullptr )
            {
                delete mDomainIntegration ;
            }
            mDomainIntegration = new IntegrationData( aElementType,
                                                aInterpolationType );

            if ( mLinearIntegration != nullptr )
            {
                delete mLinearIntegration ;
            }

            mLinearIntegration = new IntegrationData( mesh::linear_element_type( aElementType ),
                                                      aInterpolationType );
        }

//------------------------------------------------------------------------------

        void
        Calculator::set_integration_order( const uint aOrder )
        {
            mIntegrationOrder = aOrder ;
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );
            mNumberOfCornerNodes = mesh::number_of_corner_nodes( mGroup->element_type() );

            if( mDomainIntegration  != nullptr )
            {
                mDomainIntegration->populate( aOrder,
                    mGroup->parent() != nullptr ? mGroup->parent()->integration_scheme() : IntegrationScheme::GAUSS );
            }
            if( mLinearIntegration != nullptr )
            {
                mLinearIntegration->populate( aOrder,
                     mGroup->parent() != nullptr ? mGroup->parent()->integration_scheme() : IntegrationScheme::GAUSS );
            }
            this->allocate_memory() ;

            mNumberOfIntegrationPoints = mDomainIntegration->weights().length() ;

            if( mEdgeFunction != nullptr )
            {
                mEdgeFunction->precompute( mDomainIntegration->points() );
            }

            if( mGroup != nullptr )
            {
                mGroup->initialize_lookup_tables( aOrder );
            }
        }

//------------------------------------------------------------------------------

        void
        Calculator::allocate_memory()
        {
            if( mGroup->parent() == nullptr )
            {
                return  ;
            }

            if( mGroup->parent()->iwg() == nullptr )
            {
                return  ;
            }

            if(     mGroup->type() == GroupType::BLOCK ||
                    mGroup->type() == GroupType::SIDESET )
            {
                // get pointer to equation object
                IWG * tEquation = mGroup->parent()->iwg();

                // reset the list
                mVectors.clear() ;

                // reset the map
                mVectorMap.clear() ;

                // flags for special cases
                bool tHaveH = false ;

                for ( const string & tLabel: tEquation->all_fields())
                {
                    // the size of the vector
                    uint tSize = 0;

                    EntityType tType = entity_type( tLabel );

                    if( tLabel == "edge_h" )
                    {
                        tHaveH = true ;
                    }

                    // determine size
                    switch ( tType )
                    {
                        case ( EntityType::EDGE ) :
                        {
                            tSize = tEquation->edge_multiplicity() * mesh::number_of_edges( mGroup->element_type());
                            break ;
                        }
                        case ( EntityType::FACE ) :
                        {
                            tSize = tEquation->face_multiplicity() * mesh::number_of_faces( mGroup->element_type());
                            break ;
                        }
                        case ( EntityType::CELL ) :
                        {
                            tSize = tEquation->cell_multiplicity();
                            break ;
                        }
                        case ( EntityType::FACET ) :
                        {
                            tSize = tEquation->lambda_multiplicity();
                            break ;
                        }
                        default:
                        {
                            tSize = mesh::number_of_nodes( mGroup->element_type());
                        }
                    }

                    // allocate size
                    this->create_vector( tLabel, tSize, tType );
                }

                if( tHaveH )
                {
                    // clear memory
                    if( mEdgeFunction != nullptr )
                    {
                        delete mEdgeFunction ;
                        mEdgeFunction = nullptr ;
                    }
                    for( EdgeFunction * tEdgeFuction : mEdgeFunctionsMaster )
                    {
                        delete tEdgeFuction ;
                    }
                    mEdgeFunctionsMaster.clear() ;
                    for( EdgeFunction * tEdgeFuction : mEdgeFunctionsSlave )
                    {
                        delete tEdgeFuction ;
                    }
                    mEdgeFunctionsSlave.clear() ;

                    EdgeFunctionFactory tFactory;

                    switch ( mGroup->type() )
                    {
                        case ( GroupType::BLOCK ) :
                        {
                            if( mMesh->block( mGroup->id() )->has_edges() )
                            {
                                mEdgeFunction = tFactory.create_edge_function( mGroup->element_type() );
                                mEdgeFunction->precompute( mDomainIntegration->points() );
                            }

                            break ;
                        }
                        default:
                        {
                            uint tNumDimensions = mesh::dimension( mGroup->master_type() );

                            // master edge functions
                            {
                                uint n = mesh::number_of_facets( mGroup->master_type() );
                                mEdgeFunctionsMaster.set_size( n, nullptr );
                                for( uint s=0; s<n; ++s )
                                {
                                    mEdgeFunctionsMaster( s ) = tFactory.create_edge_function( mGroup->master_type() );
                                    mEdgeFunctionsMaster( s )->precompute(mGroup->master_integration(s)->points()) ;
                                }
                                this->create_vector("ntBs",  mesh::number_of_nodes( mGroup->master_type() ) );
                                this->create_vector("ntEm",  mEdgeFunctionsMaster( 0 )->ndofs() );
                                this->create_matrix("ntEm", 1, mEdgeFunctionsMaster( 0 )->ndofs() );
                                this->create_matrix("nxEm", tNumDimensions == 2 ? 1 : 3, mEdgeFunctionsMaster( 0 )->ndofs() );
                            }

                            // slave edge functions (needed for ghost and any H-H interface).
                            // Sized by the total number of slave integration
                            // permutations (cumulative over facet × orientation)
                            // so each entry is precomputed at the correct
                            // (facet, orientation) integration points and
                            // indexed the same way as mSlaveIntegration.
                            if( mGroup->slave_type() != ElementType::EMPTY )
                            {
                                uint tNumFacets = mesh::number_of_facets( mGroup->slave_type() );
                                uint tNumPermutations = 0 ;
                                for ( uint f = 0; f < tNumFacets; ++f )
                                {
                                    tNumPermutations += mesh::number_of_orientations(
                                        mGroup->slave_type(), f );
                                }
                                mEdgeFunctionsSlave.set_size( tNumPermutations, nullptr );
                                for( uint s = 0; s < tNumPermutations; ++s )
                                {
                                    mEdgeFunctionsSlave( s ) = tFactory.create_edge_function( mGroup->slave_type() );
                                    mEdgeFunctionsSlave( s )->precompute( mGroup->slave_integration( s )->points() );
                                }
                                this->create_vector("ntEs", mEdgeFunctionsSlave( 0 )->ndofs() );
                                this->create_matrix("ntEs", 1, mEdgeFunctionsSlave( 0 )->ndofs() );
                                this->create_matrix("nxEs", tNumDimensions == 2 ? 1 : 3, mEdgeFunctionsSlave( 0 )->ndofs() );
                            }
                        }
                    }
                }

                mGroup->parent()->iwg()->create_custom_vectors_and_matrices( this );

                ElementType tType = mGroup->type() == GroupType::BLOCK ? mGroup->element_type() : mGroup->master_type() ;

                switch( tType )
                {
                    case ElementType::TRI3 :
                    {
                        mFunNedelecDataH = & Calculator::nedelec_data_linear_h ;
                        mFunNedelecDataA = nullptr ;
                        break ;
                    }
                    case ElementType::TRI6 :
                    {
                        mFunNedelecDataH = & Calculator::nedelec_data_quadratic_h_2d ;
                        mFunNedelecDataA = nullptr ;
                        break ;
                    }
                    case ElementType::TET4 :
                    {
                        mFunNedelecDataH = & Calculator::nedelec_data_linear_h ;
                        mFunNedelecDataA = & Calculator::nedelec_data_linear_a ;
                        break  ;
                    }
                    case ElementType::TET10 :
                    {
                        mFunNedelecDataH = & Calculator::nedelec_data_quadratic_h_3d ;
                        mFunNedelecDataA = & Calculator::nedelec_data_quadratic_a_3d ;
                        break  ;
                    }
                    case ElementType::PENTA6TS :
                    case ElementType::QUAD4TS :
                    case ElementType::HEX8TS :
                    case ElementType::HEX8TB :
                    case ElementType::HEX8 :
                    case ElementType::QUAD4 :
                    {
                        mFunNedelecDataH = & Calculator::nedelec_data_linear_h ;
                        mFunNedelecDataA = nullptr ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Unsupported Element Type");
                    }
                }
            }
        }

//------------------------------------------------------------------------------


        uint
        Calculator::slave_integration_index_2d( const mesh::Facet * aFacet )
        {
            uint tNumOrient = mesh::number_of_orientations(
                mGroup->slave_type(), aFacet->index_on_slave() );

            if ( tNumOrient <= 1 )
            {
                // regular TRI/QUAD: one table per facet, no orientation
                return aFacet->index_on_slave() ;
            }
            else
            {
                // thin-shell QUAD: multiple orientations per facet
                return aFacet->index_on_slave() * tNumOrient
                     + aFacet->orientation_on_slave() - 1 ;
            }
        }

        uint
        Calculator::slave_integration_index_tet( const mesh::Facet * aFacet )
        {
            return aFacet->index_on_slave() * 3
                 + aFacet->orientation_on_slave() - 1 ;
        }

        uint
        Calculator::slave_integration_index_hex( const mesh::Facet * aFacet )
        {
            return aFacet->index_on_slave() * 4
                 + aFacet->orientation_on_slave() - 1 ;
        }

        uint
        Calculator::slave_integration_index_penta( const mesh::Facet * aFacet )
        {
            // PENTA (shell and volume): 3 QUAD side facets (4 orientations
            // each) followed by bottom/top TRI facets (3 orientations each).
            // Cumulative offset table = {0, 4, 8, 12, 15}. PENTA*TS now shares
            // this layout with volume PENTA (see Option-B canonicalization).
            static const uint tOffsets[] = { 0, 4, 8, 12, 15 };

            return tOffsets[ aFacet->index_on_slave() ]
                 + aFacet->orientation_on_slave() - 1 ;
        }

//------------------------------------------------------------------------------

        void
        Calculator::allocate()
        {
            // a group with only aura elements must allocate as well: the
            // postprocessors link its calculator during patch recovery
            if(  mGroup->parent() != nullptr && (
                    ( mGroup->number_of_elements() == 0
                      && mGroup->aura_elements().size() == 0 )
                    || this->integration() == nullptr ) )
            {
                return;
            }


            // make sure that allocation funciton is only called once
            BELFEM_ERROR( ! mIsAllocated, "calculator is already allocated" );
            mIsAllocated = true ;

            // build the timestep history table here rather than in
            // link( Group * ): allocate() is reached by BOTH block and
            // sideset calculators via DofManager::init_work(), after the
            // IWG timestepping method has been configured
            this->init_qold_table();

            // build the maxwell data helper if the kernels have been
            // registered via link_maxwell(). This must happen BEFORE the
            // "done if this is a block" early return below — blocks are
            // exactly the groups that need the helper. Work vectors are no
            // precondition: MaxwellData::link_vector() creates missing ones
            // on demand.
            if ( mMaxwellKernel != nullptr
                 && mGroup->type() == GroupType::BLOCK
                 && mGroup->domain_type() != DomainType::Air
                 // an aura-only block may be absent from the peer dof manager,
                 // whose aura is expanded from a different owned-element set
                 && mMaxwellKernel->dofmgr()->block_exists( mGroup->id() )
                 // side connector walls are never part of the thermal kernel
                 // but still need the helper ( the ctor nulls the thermal
                 // peer via its own block_exists check )
                 && ( mThermalKernel == nullptr
                      || mThermalKernel->dofmgr()->block_exists( mGroup->id() )
                      || mGroup->domain_type() == DomainType::LeftCoating
                      || mGroup->domain_type() == DomainType::RightCoating ) )
            {
                mMaxwellData = new calculator::MaxwellData( this, mMaxwellKernel, mThermalKernel );

                this->select_link_element_dispatcher();
            }
            else
            {
                mFunLinkElement = & Calculator::link_element_default;
            }

            // number of nodes per element
            mNumberOfNodes = mesh::number_of_nodes( mGroup->element_type() );

            // corner nodes per element
            mNumberOfCornerNodes = mesh::number_of_corner_nodes( mGroup->element_type() );

            // number of edges per element
            //uint tNumEdges = mesh::number_of_edges( mGroup->element_type() );

            // number of faces per element
            //uint tNumFaces = mesh::number_of_faces( mGroup->element_type() );

            // uint tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->element_type() );

            // get the number of dimensions
            uint tNumDimensions = mMesh->number_of_dimensions() ;

            // get the number of dofs from first element. a group with only
            // aura elements has no dofs on this proc; the recovery paths
            // never touch the dof-sized work arrays
            uint tNumDofs = mGroup->parent() != nullptr && mGroup->number_of_elements() > 0 ?
                mGroup->elements()(0)->number_of_local_dofs() : 0 ;

            // matrices for X-Coordinates
            mX.set_size( mNumberOfNodes, tNumDimensions, BELFEM_QUIET_NAN );
            mXc.set_size( mNumberOfCornerNodes, tNumDimensions, BELFEM_QUIET_NAN );

            // matrix for Jacobian and its inverse
            mJ = this->create_matrix( "J", tNumDimensions, tNumDimensions );
            mInvJ = this->create_matrix( "InvJ", tNumDimensions, tNumDimensions );

            // special node and matrix interpolators
            switch( mGroup->parent() == nullptr ? IwgType::UNDEFINED : mGroup->parent()->iwg()->type() )
            {
                case( IwgType::Gradient2D ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 2, tNumDofs );
                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );

                    // link functions
                    mFunN = & Calculator::N2D ;
                    mFunB = & Calculator::Bscalar ;


                    break;
                }
                case( IwgType::Gradient3D ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 3, tNumDofs );

                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );

                    // link functions
                    mFunN = & Calculator::N3D ;
                    mFunB = & Calculator::Bscalar ;

                    break;
                }
                case( IwgType::PlaneStress ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 2, tNumDofs );

                    // matrices for gradient operator
                    mdN = this->create_matrix( "dN", 2, 2*mNumberOfNodes );
                    mB = this->create_matrix( "B", 3, tNumDofs );

                    // link functions
                    mFunN = & Calculator::N2D ;
                    mFunB = & Calculator::Bplanestress ;
                    break ;
                }
                case( IwgType::LinearElasticity ) :
                {
                    // matrix for node interpolation
                    mN = this->create_matrix( "N", 3, tNumDofs );

                    // matrices for gradient operator
                    mdN = this->create_matrix( "dN", 3, 3*mNumberOfNodes );
                    mB  = this->create_matrix( "B", 6, tNumDofs );

                    // link functions
                    mFunN = & Calculator::N3D ;
                    mFunB = & Calculator::Bvoigt ;
                    break ;
                }
                default:
                {
                    mN = this->create_matrix( "N", 1, mNumberOfNodes );
                    mB = this->create_matrix( "B", tNumDimensions, mNumberOfNodes );

                    mN->matrix().fill( 0.0 );
                    mB->matrix().fill( 0.0 );

                    mFunN = & Calculator::Nscalar ;
                    mFunB = & Calculator::Bscalar ;

                    if ( mGroup->type() == GroupType::BLOCK )
                    {
                        break ;
                    }
                    // set functions for master element
                    if( mGroup->master_type() != ElementType::EMPTY )
                    {
                        uint tNumMasterNodes = mesh::number_of_nodes( mGroup->master_type() ) ;

                        mNm = this->create_matrix( "Nm", 1, tNumMasterNodes );
                        mBm = this->create_matrix( "Bm", tNumDimensions, tNumMasterNodes );
                        mInvJm = this->create_matrix( "InvJm", tNumDimensions, tNumDimensions );

                        mFunNm = & Calculator::Nscalar_master ;
                        mFunBm = & Calculator::Bscalar_master ;

                        mntBm  = this->create_matrix("ntBm", 1, tNumMasterNodes );
                        mnxBm  = this->create_matrix( "nxBm", tNumDimensions == 3 ? 3 : 1, tNumMasterNodes );
                    }

                    // set functions for slave element
                    if( mGroup->slave_type() != ElementType::EMPTY  )
                    {
                        uint tNumSlaveNodes = mesh::number_of_nodes( mGroup->slave_type() );

                        mNs = this->create_matrix( "Ns", 1, tNumSlaveNodes );
                        mBs = this->create_matrix( "Bs", tNumDimensions, tNumSlaveNodes );
                        mInvJs = this->create_matrix( "InvJs", tNumDimensions, tNumDimensions );

                        mFunNs = & Calculator::Nscalar_slave ;
                        mFunBs = & Calculator::Bscalar_slave ;

                        mntBs  = this->create_matrix( "ntBs", 1, tNumSlaveNodes );
                        mnxBs  = this->create_matrix( "nxBs", tNumDimensions == 3 ? 3 : 1, tNumSlaveNodes );
                    }

                    break ;
                }
            }

            // stiffness matrix
            mK.set_size( tNumDofs, tNumDofs, BELFEM_QUIET_NAN );

            // Newton correction matrix
            mJN.set_size( tNumDofs, tNumDofs, BELFEM_QUIET_NAN );

            // load vector
            mf.set_size( tNumDofs, BELFEM_QUIET_NAN );

            // dof vectors
            mq0.set_size( tNumDofs, BELFEM_QUIET_NAN );
            mq.set_size( tNumDofs, BELFEM_QUIET_NAN );
            mqswap.set_size( tNumDofs, BELFEM_QUIET_NAN );

            // link function to invert J
            switch( tNumDimensions )
            {
                case( 2 ) :
                {
                    mFunInvertJ = & inv2 ;
                    break ;
                }
                case( 3 ) :
                {
                    mFunInvertJ = & inv3 ;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid number of dimensions") ;
                }
            }

            // done if this is a block
            if ( mGroup->type() == GroupType::BLOCK )
            {
                switch( mGroup->parent()->iwg()->model_dimensionality() )
                {
                    case( ModelDimensionality::TwoD ) :
                    {
                        if( mesh::geometry_type( mGroup->element_type() ) == GeometryType::TRI )
                        {
                            mFundV = & Calculator::dV_tri6_tet10 ;
                            break ;
                        }
                        else
                        {
                            if ( mGroup->element_type() == ElementType::QUAD4TS )
                            {
                                // | det J | : these layer elements are wound
                                // clockwise by construction, see dV_quad4ts.
                                // The magnetic solve reaches the wrapper too,
                                // via its own Calculator, but takes the
                                // edge-function branch, where the value is
                                // positive for valid input and abs is the
                                // identity
                                mFundV = & Calculator::dV_quad4ts ;
                            }
                            else if ( mGroup->element_type() == ElementType::QUAD4 )
                            {
                                // bilinear quad: det J varies over the
                                // element, and dV_hex caches it per
                                // integration point ( key aIndex, not 0 ).
                                // Its edge-function branch takes the null
                                // path here: the edge function factory
                                // cannot build a plain QUAD4
                                mFundV = & Calculator::dV_hex ;
                            }
                            else
                            {
                                BELFEM_ERROR( false, "Higher order thin shells are not implemented!");
                            }
                        }
                        break ;
                    }
                    case( ModelDimensionality::AxSymmX ) :
                    {
                        mFundV = & Calculator::dV_axsymmx ;
                        break ;
                    }
                    case( ModelDimensionality::AxSymmY ) :
                    {
                        mFundV = & Calculator::dV_axsymmy ;
                        break ;
                    }
                    case( ModelDimensionality::ThreeD ) :
                    {
                        switch( mesh::geometry_type( mGroup->element_type() ) )
                        {
                            case GeometryType::TET :
                            {
                                mFundV = & Calculator::dV_tri6_tet10 ;
                                break ;
                            }
                            case GeometryType::PENTA :
                            {
                                if ( mGroup->element_type() == ElementType::PENTA6TS )
                                {
                                    mFundV = & Calculator::dV_ts ;
                                }
                                else
                                {
                                    BELFEM_ERROR( false, "Higher order thin shells are not implemented!");
                                }
                                break ;
                            }
                            case GeometryType::HEX :
                            {
                                if ( mGroup->element_type() == ElementType::HEX8TS )
                                {
                                    mFundV = & Calculator::dV_ts ;
                                }
                                else
                                {
                                    mFundV = & Calculator::dV_hex ;
                                }
                                break ;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "You haven't implemented the quads yet!");
                            }
                        }
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid Model Dimensionality");
                    }
                }
                switch( mGroup->parent() != nullptr ? mGroup->parent()->iwg()->model_dimensionality() : ModelDimensionality::ThreeD )
                {
                    case( ModelDimensionality::TwoD ) :
                    {
                        mFunInvJ = & Calculator::invJ2D3D ;
                        break ;
                    }
                    case( ModelDimensionality::AxSymmX ) :
                    {
                        mFunInvJ = & Calculator::invJaxsym ;
                        break ;
                    }
                    case( ModelDimensionality::AxSymmY ) :
                    {
                        mFunInvJ = & Calculator::invJaxsym ;
                        break ;
                    }
                    case( ModelDimensionality::ThreeD ) :
                    {
                        mFunInvJ = & Calculator::invJ2D3D ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid Model Dimensionality");
                    }
                }

                return ;
            }

            switch( mGroup->parent() != nullptr ? mGroup->parent()->iwg()->model_dimensionality() : ModelDimensionality::ThreeD )
            {
                case( ModelDimensionality::TwoD ) :
                case( ModelDimensionality::ThreeD ) :
                {
                    mFundS = & Calculator::dS_cartesian ;
                    break ;
                }
                case( ModelDimensionality::AxSymmX ) :
                {
                    mFundS = & Calculator::dS_axsymmx ;
                    break ;
                }
                case( ModelDimensionality::AxSymmY ) :
                {
                    mFundS = & Calculator::dS_axsymmy ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Model Dimensionality");
                }
            }

            // check if we are allocating master and slave elements
            if( mGroup->master_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->master_type() );

                mXm.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                // allocate normal vector
                mNormal.set_size( tNumDimensions, BELFEM_QUIET_NAN );

                mJm = this->create_matrix( "Jm", tNumDimensions, tNumDimensions );

                // get the interpolation order of the master block
                InterpolationOrder tOrder = mesh::interpolation_order( mGroup->master_type() ) ;

                switch( mesh::geometry_type( mGroup->master_type() ) )
                {
                    case GeometryType::TRI :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_tri_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_tri_curved ;
                        }
                        break ;
                    }
                    case GeometryType::QUAD :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_quad_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_quad_curved ;
                        }
                        break ;
                    }
                    case GeometryType::TET :
                    {
                        if( tOrder == InterpolationOrder::LINEAR )
                        {
                            mFunNormal = & Calculator::normal_tet_straight ;
                        }
                        else
                        {
                            mFunNormal = & Calculator::normal_tet_curved ;
                        }
                        break ;
                    }
                    case GeometryType::HEX :
                    {
                        mFunNormal = & Calculator::normal_hex ;
                        break ;
                    }
                    case GeometryType::PENTA :
                    {
                        // Thin-shell PENTA*TS shares facet topology with
                        // volume PENTA (5 facets, CCW outward normals), so
                        // both dispatch through normal_penta.
                        mFunNormal = & Calculator::normal_penta ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "No normal function assigned");
                    }
                }
            }

            if( mGroup->slave_type() != ElementType::EMPTY )
            {
                uint tNumNodes = mesh::number_of_nodes( mGroup->slave_type() );

                mXs.set_size( tNumNodes, tNumDimensions, BELFEM_QUIET_NAN );

                //tNumNedelecDofs = mesh::number_of_nedelec_dofs( mGroup->slave_type() );

                // mIndexXs = this->create_matrix( "Xs", tNumNodes, tNumDimensions );
                // mIndexBs = this->create_matrix( "Bs", tNumNodes, tNumDimensions );

                mJs = this->create_matrix( "Js", tNumDimensions, tNumDimensions );

                switch( mesh::geometry_type( mGroup->slave_type() ) )
                {
                    case( GeometryType::TRI ) :
                    case( GeometryType::QUAD ) :
                    {
                        mFunSlaveIntegrationIndex = & Calculator::slave_integration_index_2d  ;
                        break ;
                    }
                    case( GeometryType::TET ) :
                    {
                        mFunSlaveIntegrationIndex = & Calculator::slave_integration_index_tet  ;
                        break ;
                    }
                    case( GeometryType::HEX ) :
                    {
                        mFunSlaveIntegrationIndex = & Calculator::slave_integration_index_hex  ;
                        break ;
                    }
                    case( GeometryType::PENTA ) :
                    {
                        mFunSlaveIntegrationIndex = & Calculator::slave_integration_index_penta ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "unsupported type for slave side of facets");
                    }
                }
            }
            else
            {
                mFunSlaveIntegrationIndex = nullptr ;
            }

        }

//-----------------------------------------------------------------------------

        void
        Calculator::select_link_element_dispatcher()
        {
            if ( mThermalKernel == nullptr
                 // side connector walls have no thermal peer block — their
                 // temperature comes from the seam nodes, so they must link
                 // like a plain maxwell group even on coupled runs
                 || mGroup->domain_type() == DomainType::LeftCoating
                 || mGroup->domain_type() == DomainType::RightCoating )
            {
                mFunLinkElement = & Calculator::link_element_maxwell;
            }
            else if ( mGroup->parent()->iwg()->type() == IwgType::Maxwell )
            {
                mFunLinkElement = & Calculator::link_element_maxwell_thermal;
            }
            else if ( mGroup->parent()->iwg()->type() == IwgType::MaxwellThermal )
            {
                mFunLinkElement = & Calculator::link_element_thermal_maxwell;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid IwgType" );
            }
        }

//-----------------------------------------------------------------------------

        void
        Calculator::link_maxwell( Kernel * aMaxwellKernel, Kernel * aThermalKernel )
        {
            mMaxwellKernel = aMaxwellKernel ;
            mThermalKernel = aThermalKernel ;

            // the helper only exists on block groups; sidesets keep a nullptr
            // and the maxwell_data() assert catches misuse
            if ( mGroup->type() != GroupType::BLOCK || mGroup->domain_type() == DomainType::Air )
            {
                return ;
            }

            // if the calculator is already allocated, the work vectors exist
            // and we can (re)build immediately, e.g. when the thermal kernel
            // is attached after initialization
            if ( mIsAllocated )
            {
                if ( mMaxwellData != nullptr )
                {
                    delete mMaxwellData ;
                }
                mMaxwellData = new calculator::MaxwellData( this, mMaxwellKernel, mThermalKernel );

                // the dispatcher choice depends on mThermalKernel, which may
                // have changed since allocate() ran (the thermal kernel
                // initializes before Controller::set_thermal_kernel relinks us)
                this->select_link_element_dispatcher();
            }
        }

//-----------------------------------------------------------------------------

        void
        Calculator::link( Group * aGroup )
        {
            mGroup = aGroup ;

            if( mGroup->master_type() != ElementType::EMPTY )
            {
                uint n = mEdgeFunctionsMaster.size() ;
                for( uint s=0; s<n; ++s )
                {
                    mEdgeFunctionsMaster( s )->precompute( mGroup->master_integration( s )->points() );
                }
            }

            if( mGroup->slave_type() != ElementType::EMPTY )
            {
                uint n = mEdgeFunctionsSlave.size() ;
                for( uint s=0; s<n; ++s )
                {
                    mEdgeFunctionsSlave( s )->precompute( mGroup->slave_integration( s )->points() );
                }
            }

            mIsLinear = mesh::interpolation_order( mGroup->element_type() ) == InterpolationOrder::LINEAR ;
        }

//-----------------------------------------------------------------------------

        void
        Calculator::link( Element * aElement )
        {
            BELFEM_ERROR( mIsAllocated, "calculator is not allocated");

            // reset vectors
            for( calculator::VectorData * tVector : mVectors )
            {
                tVector->set_index( BELFEM_UINT_MAX );
            }

            // reset matrices
            for( calculator::MatrixData * tMatrix : mMatrices )
            {
                tMatrix->set_index( BELFEM_UINT_MAX );
            }

            // reset determinant
            mDetJ     = BELFEM_QUIET_NAN ;
            mDetJIndex = BELFEM_UINT_MAX ;
            mSurfaceIncrement = BELFEM_QUIET_NAN ;

            mNormalIndex = BELFEM_UINT_MAX ;

            aElement->get_node_coors( mX );

            mIsCurved = aElement->element()->is_curved() ;
            mIsLinear = mesh::interpolation_order( aElement->element()->type() ) == InterpolationOrder::LINEAR ;

            // for curved elements, we can copy the corner nodes
            if( ! mIsCurved )
            {
                for( uint j=0; j<mX.n_cols(); ++j )
                {
                    for( uint i=0; i<mNumberOfCornerNodes; ++i )
                    {
                        mXc( i, j ) = mX( i, j );
                    }
                }
            }

            mB->set_index( BELFEM_UINT_MAX );
            mJ->set_index( BELFEM_UINT_MAX );
            if( aElement->master() != nullptr )
            {
                aElement->master()->get_node_coors( mXm );
                mMasterIndex = aElement->facet()->index_on_master() ;
                mMasterIntegration = mGroup->master_integration( mMasterIndex );
                if ( mEdgeFunctionsMaster.size() > 0 )
                {
                    mEdgeFunctionMaster = mEdgeFunctionsMaster( mMasterIndex );
                }

                mJm->set_index( BELFEM_UINT_MAX );
                mBm->set_index( BELFEM_UINT_MAX );
                mNumberOfNodesOnMaster = aElement->master()->element()->number_of_nodes() ;
                mNumberOfEdgesOnMaster = aElement->master()->element()->number_of_edges() ;
                mNumberOfFacesOnMaster = aElement->master()->element()->number_of_faces() ;

                if( mGroup->parent()->iwg()->enrich_sidesets() )
                {
                    Block * tMasterBlock = mGroup->parent()->block( aElement->master()->element()->block_id() ) ;

                    mMasterVolumeIntegration = tMasterBlock->integration();
                    mVolumeEnrichment        = tMasterBlock->enrichment_data( mMasterIndex );
                    mSideSetEnrichment       = mGroup->enrichment_data( mMasterIndex );
                }
            }

            uint tSlaveIntIndex = BELFEM_UINT_MAX ;
            if( aElement->slave() != nullptr )
            {
                aElement->slave()->get_node_coors( mXs );
                tSlaveIntIndex = ( this->*mFunSlaveIntegrationIndex )( aElement->facet() );
                mSlaveIntegration = mGroup->slave_integration( tSlaveIntIndex );
                mJs->set_index( BELFEM_UINT_MAX );
                mBs->set_index( BELFEM_UINT_MAX );
                mNumberOfNodesOnSlave  = aElement->slave()->element()->number_of_nodes() ;
                mNumberOfEdgesOnSlave  = aElement->slave()->element()->number_of_edges() ;
                mNumberOfFacesOnSlave  = aElement->slave()->element()->number_of_faces() ;
            }

            if( mEdgeFunction != nullptr )
            {
                mEdgeFunction->link( aElement );
            }

            if( mEdgeFunctionsMaster.size() > 0 && aElement->facet () != nullptr )
            {
                mEdgeFunctionMaster = mEdgeFunctionsMaster( aElement->facet()->index_on_master() );

                BELFEM_ASSERT( aElement->master() != nullptr,
                    "master element of %lu (block %u, %s ), owned by %u is null.",
                    ( long unsigned int ) aElement->id(),
                    ( long unsigned int ) aElement->element()->block_id(),
                    to_string( aElement->element()->type() ).c_str(),
                    ( unsigned int ) aElement->element()->owner() );

                mEdgeFunctionMaster->link( aElement->master() );

                if( mEdgeFunctionsSlave.size() > 0 && aElement->facet()->has_slave() )
                {
                    // Use the flattened (facet, orientation) index that
                    // already resolved to mSlaveIntegration above, so the
                    // edge function is the one precomputed at the matching
                    // integration points.
                    BELFEM_ASSERT( tSlaveIntIndex != BELFEM_UINT_MAX,
                        "slave integration index was not resolved" );
                    mEdgeFunctionSlave = mEdgeFunctionsSlave( tSlaveIntIndex );
                    mEdgeFunctionSlave->link( aElement->slave() );
                }
            }

            ( this->*mFunLinkElement )( aElement );
        }

        void
        Calculator::link( mesh::Facet * aFacet )
        {
            BELFEM_ASSERT( aFacet->master() != nullptr, "Expect a master element");
            BELFEM_ASSERT( aFacet->slave() != nullptr, "Expect a slave element");
            // grab coordinates of master
            for ( uint j=0; j<3; ++j )
            {
                for( uint i=0; i<aFacet->master()->number_of_nodes(); ++i )
                {
                    mXm( i, j ) = aFacet->master()->node( i )->x( j );
                }
            }

            // grab coordinates of slave
            for ( uint j=0; j<3; ++j )
            {
                for( uint i=0; i<aFacet->slave()->number_of_nodes(); ++i )
                {
                    mXs( i, j ) = aFacet->slave()->node( i )->x( j );
                }
            }

            // grab coordinates of facet
            // grab coordinates of slave
            for ( uint j=0; j<3; ++j )
            {
                for( uint i=0; i<aFacet->number_of_nodes(); ++i )
                {
                    mX( i, j ) = aFacet->node( i )->x( j );
                }
            }

            // reset vectors
            for( calculator::VectorData * tVector : mVectors )
            {
                tVector->set_index( BELFEM_UINT_MAX );
            }

            // reset matrices
            for( calculator::MatrixData * tMatrix : mMatrices )
            {
                tMatrix->set_index( BELFEM_UINT_MAX );
            }

            // reset determinant
            mDetJ     = BELFEM_QUIET_NAN ;
            mDetJIndex = BELFEM_UINT_MAX ;
            mSurfaceIncrement = BELFEM_QUIET_NAN ;

            mNormalIndex = BELFEM_UINT_MAX ;
            mB->set_index( BELFEM_UINT_MAX );
            mJ->set_index( BELFEM_UINT_MAX );
            mInvJ->set_index( BELFEM_UINT_MAX );
            mBm->set_index( BELFEM_UINT_MAX );
            mJm->set_index( BELFEM_UINT_MAX );
            mInvJm->set_index( BELFEM_UINT_MAX );
            mBs->set_index( BELFEM_UINT_MAX );
            mJs->set_index( BELFEM_UINT_MAX );
            mInvJs->set_index( BELFEM_UINT_MAX );
            mMasterIndex = aFacet->index_on_master() ;
            mMasterIntegration = mGroup->master_integration( aFacet->index_on_master() );

            BELFEM_ASSERT( mFunSlaveIntegrationIndex != nullptr,
                "slave integration index function not set" );
            mSlaveIntegration = mGroup->slave_integration(
                ( this->*mFunSlaveIntegrationIndex )( aFacet ) );

            mIsCurved = aFacet->master()->is_curved() || aFacet->slave()->is_curved() ;
            mIsLinear = mesh::interpolation_order( aFacet->element()->type() ) == InterpolationOrder::LINEAR ;
        }


//-----------------------------------------------------------------------------

        Calculator *
        Calculator::get_normal_calculator(
                Vector< real > & aPhiM,
                Vector< real > & aPhiS,
                bool           & aMasterIsConductor,
                bool           & aSlaveIsConductor )
        {
            // get the mesh facet
            mesh::Facet * tMeshFacet = this->element()->facet();

            BELFEM_ASSERT( tMeshFacet != nullptr, "Expect a mesh facet for element %lu",
                ( long unsigned int ) this->element()->id() );

            // get the fem sideset: the group's parent is the dof manager
            // that owns this calculator ( the controller detour used before
            // hard-errors on a controller-less kernel, e.g. the test fixtures )
            SideSet * tSideSet = this->group()->parent()->sideset( tMeshFacet->sideset_id() );

            // get the fem facet
            Element * tFacet = tSideSet->element( tMeshFacet->id() );

            // get the other calculator
            Calculator * aNormalCalc = tSideSet->calculator() ;

            BELFEM_ASSERT(  aNormalCalc != nullptr, "Expect a normal calculator for sideset %u", ( uint )tMeshFacet->element()->block_id() );

            BELFEM_ASSERT(  tFacet != nullptr, "Expect a facet for sideset %u", ( uint )tMeshFacet->element()->block_id() );

            BELFEM_ASSERT(  tMeshFacet->master() != nullptr, "Expect a master element for facet %lu", ( long unsigned int ) tFacet->id() );
            BELFEM_ASSERT(  tMeshFacet->slave()  != nullptr, "Expect a slave element for facet %lu", ( long unsigned int ) tFacet->id() );

            // link calculator and facet
            aNormalCalc->link(  tFacet );

            // get master element
            mesh::Element * tMaster = tMeshFacet->master();

            // get the slave element
            mesh::Element * tSlave = tMeshFacet->slave() ;

            // an h-conductor side has no potential: its trace comes from
            // the edge field ( compute_h_trace ), the nodal phi there is
            // bookkeeping and must not be read
            aMasterIsConductor = this->volume_is_conductor( tMaster );
            aSlaveIsConductor  = this->volume_is_conductor( tSlave );

            const Vector< real > & tPhi = this->group()->parent()->mesh()->field_data( "phi" );

            if ( ! aMasterIsConductor )
            {
                aPhiM.set_size(  tMaster->number_of_nodes() );
                for ( uint k=0; k<tMaster->number_of_nodes(); ++k )
                {
                    aPhiM( k ) = tPhi( tMaster->node( k )->index() );
                }
            }

            if ( ! aSlaveIsConductor )
            {
                aPhiS.set_size(  tSlave->number_of_nodes() );
                for ( uint k=0; k<tSlave->number_of_nodes(); ++k )
                {
                    aPhiS( k ) = tPhi( tSlave->node( k )->index() );
                }
            }
            return aNormalCalc ;
        }

//------------------------------------------------------------------------------

        bool
        Calculator::volume_is_conductor( const mesh::Element * aVolume ) const
        {
            return mMesh->block( aVolume->block_id() )->domain_type() == DomainType::Conductor ;
        }

//------------------------------------------------------------------------------

        calculator::VectorData *
        Calculator::create_vector( const string & aLabel, const uint aSize, const EntityType aType )
        {

            calculator::VectorData * aVector = nullptr ;

            if( mVectorMap.key_exists( aLabel ) )
            {
                aVector = mVectorMap( aLabel );
                aVector->vector().set_size( aSize );
                BELFEM_ASSERT( aVector->entity_type() == aType,
                    "Vector %s has already been created but is of wrong entity type",
                    aLabel.c_str() );

            }
            else
            {
                aVector = new calculator::VectorData( aLabel, aSize, aType );

                // add vector to user vectors
                mVectors.push( aVector );

                // add vector to map
                mVectorMap[ aLabel ] = aVector ;
            }

            return aVector ;
        }

//------------------------------------------------------------------------------

        calculator::MatrixData *
        Calculator::create_matrix( const string & aLabel,
                             const uint aNumRows,
                             const uint aNumCols )
        {

            calculator::MatrixData * aMatrix = nullptr ;

            if( mMatrixMap.key_exists( aLabel ) )
            {
                aMatrix = mMatrixMap( aLabel );
                aMatrix->matrix().set_size( aNumRows, aNumCols );
            }
            else
            {
                aMatrix = new calculator::MatrixData( aLabel,
                                                               aNumRows,
                                                               aNumCols );

                // add vector to user vectors
                mMatrices.push( aMatrix );

                // add vector to map
                mMatrixMap[ aLabel ] = aMatrix ;
            }

            // fill matrix with zeros
            aMatrix->matrix().fill( 0.0 );

            return aMatrix ;
        }

//------------------------------------------------------------------------------

        void
        Calculator::print_dofs()
        {
            mGroup->parent()->iwg()->print_dofs( mElement, false );
        }

//------------------------------------------------------------------------------

        void
        Calculator::print_local_dofs()
        {
            mGroup->parent()->iwg()->print_dofs( mElement, true );
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_linear( const string & aEdgeField )
        {
            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str() );

            // grab data object
            calculator::VectorData * tVectorData = mVectorMap( aEdgeField );

            // get ref to field on mesh
            Vector< real > & tField = mMesh->field_data( aEdgeField );

            // grab the vector object
            Vector< real > & aData = tVectorData->vector() ;

            // loop over all edges
            for( uint e=0; e< mElement->element()->number_of_edges(); ++e )
            {
                aData( e ) = tField( mElement->element()->edge( e )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_master_h()
        {
            BELFEM_ASSERT( mElement->master() != nullptr,
                "Element %lu has no master", ( long unsigned int ) mElement->id() );

            // the edge field is the live dof storage on every rank ( q()
            // reads it too ), and the only source for aura masters
            const Vector< real > & tField = mMesh->field_data( "edge_h" );

            // sized for the master type by create_custom_vectors_and_matrices
            Vector< real > & aData = mVectorMap( "nedelec_h" )->vector() ;

            mesh::Element * tMaster = mElement->master()->element() ;

            BELFEM_ASSERT( aData.length() == tMaster->number_of_edges(),
                "nedelec_data_master_h expects one dof per edge ( %u dofs vs %u edges )",
                ( unsigned int ) aData.length(),
                ( unsigned int ) tMaster->number_of_edges() );

            for( uint e=0; e< tMaster->number_of_edges(); ++e )
            {
                aData( e ) = tField( tMaster->edge( e )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_slave_h()
        {
            BELFEM_ASSERT( mElement->slave() != nullptr,
                "Element %lu has no slave", ( long unsigned int ) mElement->id() );

            const Vector< real > & tField = mMesh->field_data( "edge_h" );

            // sized for the slave type by create_custom_vectors_and_matrices
            Vector< real > & aData = mVectorMap( "nedelec_h_s" )->vector() ;

            mesh::Element * tSlave = mElement->slave()->element() ;

            BELFEM_ASSERT( aData.length() == tSlave->number_of_edges(),
                "nedelec_data_slave_h expects one dof per edge ( %u dofs vs %u edges )",
                ( unsigned int ) aData.length(),
                ( unsigned int ) tSlave->number_of_edges() );

            for( uint e=0; e< tSlave->number_of_edges(); ++e )
            {
                aData( e ) = tField( tSlave->edge( e )->index() );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_quadratic_2d(
                const string & aEdgeField,
                const string & aFaceField,
                const string & aVectorLabel )
        {

            BELFEM_ASSERT(
                    mMesh->field( aEdgeField )->entity_type() == EntityType::EDGE,
                    "Field '%s' is not an edge field", aEdgeField.c_str());


            BELFEM_ASSERT(
                    mMesh->field( aFaceField )->entity_type() == EntityType::FACE,
                    "Field '%s' is not an edge field", aFaceField.c_str());


            Vector< real > & aData = this->vector(aVectorLabel );

            Vector< real > & tEdgeData = mMesh->field_data( aEdgeField );
            Vector< real > & tFaceData = mMesh->field_data( aFaceField );

            uint tCount = 0;

            for ( uint e = 0; e < mElement->element()->number_of_edges(); ++e )
            {

                // get index of edge
                index_t tIndex = mElement->element()->edge( e )->index();

                // check direction of edge
                if ( mElement->edge_direction( e ))
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                }
                else
                {
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex + 1 );
                    aData( tCount++ ) = tEdgeData( tIndex + tIndex );
                }
            }

            // write data into container
            index_t tIndex = 2 * mElement->element()->index();
            aData( tCount++ ) = tFaceData( tIndex );
            aData( tCount++ ) = tFaceData( tIndex + 1 );

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::nedelec_data_quadratic_3d(
            const string & aEdgeField,
            const string & aFaceField,
            const string & aVectorLabel )
        {
            // get ref to edge field on mesh
            Vector< real > & tEdgeField = mMesh->field_data( aEdgeField );

            // get ref to face field on mesh
            Vector< real > & tFaceField = mMesh->field_data( aFaceField );

            // the result vector
            Vector< real > & aData = this->vector( aVectorLabel );

            // initialize counter
            uint tCount = 0 ;

            uint tNumEdges = mElement->element()->number_of_edges() ;
            uint tNumFaces = mElement->element()->number_of_faces() ;

            // check length of memory container
            BELFEM_ASSERT(
                    aData.length() >= 2 * ( tNumEdges + tNumFaces ),
                    "Length of vector does not fit ( is %u, but need at least %u )",
                    ( unsigned int ) aData.length(),
                    ( unsigned int ) 2 * ( tNumEdges + tNumFaces ) );

            for( uint e=0; e< tNumEdges; ++e )
            {

                // get index of edge
                index_t tIndex = mElement->element()->edge( e )->index() ;

                // check direction of edge
                if( mElement->edge_direction( e ) )
                {
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex );
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex + 1 );
                }
                else
                {
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex + 1 );
                    aData( tCount ++ ) = tEdgeField( tIndex + tIndex  );
                }
            }

            // face dofs orientation are handeled by modifying the shape functions
            // therefore, we just populate here
            for( uint f=0; f<tNumFaces; ++f )
            {
                // get index of face
                index_t tIndex = mElement->element()->face( f )->index() ;

                // write data into container
                aData( tCount++ ) = tFaceField( tIndex + tIndex );
                aData( tCount++ ) = tFaceField( tIndex + tIndex + 1 );
            }

            return aData ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tri_straight( const uint aIndex  )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                // remember the index
                mNormalIndex = aIndex;

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) = mXm( 1, 1 ) - mXm( 0, 1 );
                        mNormal( 1 ) = mXm( 0, 0 ) - mXm( 1, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {

                        mNormal( 0 ) = mXm( 2, 1 ) - mXm( 1, 1 );
                        mNormal( 1 ) = mXm( 1, 0 ) - mXm( 2, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = mXm( 0, 1 ) - mXm( 2, 1 );
                        mNormal( 1 ) = mXm( 2, 0 ) - mXm( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tri_curved( const uint aIndex  )
        {
            if( ! mIsCurved )
            {
                return this->normal_tri_straight( aIndex );
            }
            else if( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix
                const Matrix< real > & J = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) =
                                J( 1, 1 )
                                - J( 0, 1 );

                        mNormal( 1 ) =
                                J( 0, 0 )
                                - J( 1, 0 );

                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) = -J( 1, 1 );
                        mNormal( 1 ) = J( 1, 0 );

                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = J( 0, 1 );
                        mNormal( 1 ) = -J( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // if this edge was straight, this would be the length of this side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_quad_straight( const uint aIndex )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                // remember the index
                mNormalIndex = aIndex;

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) = mXm( 1, 1 ) - mXm( 0, 1 );
                        mNormal( 1 ) = mXm( 0, 0 ) - mXm( 1, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) = mXm( 2, 1 ) - mXm( 1, 1 );
                        mNormal( 1 ) = mXm( 1, 0 ) - mXm( 2, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = mXm( 3, 1 ) - mXm( 2, 1 );
                        mNormal( 1 ) = mXm( 2, 0 ) - mXm( 3, 0 );
                        break;
                    }
                    case ( 3 ) :
                    {
                        mNormal( 0 ) = mXm( 0, 1 ) - mXm( 3, 1 );
                        mNormal( 1 ) = mXm( 3, 0 ) - mXm( 0, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // finally, we must adapt this value, since along the edge
                // we integrate from -1 to 1 rather than 0 to 1
                mSurfaceIncrement *= 0.5;
            }

            // now we can return the vector
            return mNormal;
        }
//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_quad_curved( const uint aIndex )
        {
            if( ! mIsCurved )
            {
                return this->normal_quad_straight( aIndex );
            }
            else if( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix
                const Matrix< real > & J = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case ( 0 ) :
                    {
                        mNormal( 0 ) =  J( 0, 1 );
                        mNormal( 1 ) = -J( 0, 0 );
                        break;
                    }
                    case ( 1 ) :
                    {
                        mNormal( 0 ) =  J( 1, 1 );
                        mNormal( 1 ) = -J( 1, 0 );
                        break;
                    }
                    case ( 2 ) :
                    {
                        mNormal( 0 ) = -J( 0, 1 );
                        mNormal( 1 ) =  J( 0, 0 );
                        break;
                    }
                    case ( 3 ) :
                    {
                        mNormal( 0 ) = -J( 1, 1 );
                        mNormal( 1 ) =  J( 1, 0 );
                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this value will contain the length of the side
                mSurfaceIncrement = norm( mNormal );

                // now let's norm the vector
                mNormal /= mSurfaceIncrement;

                // no multiplication of mSurfaceIncrement with 0.5 here!
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tet_straight( const uint aIndex )
        {
            if( mNormalIndex == BELFEM_UINT_MAX )
            {
                mNormalIndex = aIndex ;
                switch ( mMasterIndex )
                {
                    case 0 :
                    {
                        mNormal( 0 ) = (mXm(0,1)-mXm(1,1))*(mXm(0,2)-mXm(3,2))-(mXm(0,2)-mXm(1,2))*(mXm(0,1)-mXm(3,1));
                        mNormal( 1 ) = (mXm(0,2)-mXm(1,2))*(mXm(0,0)-mXm(3,0))-(mXm(0,0)-mXm(1,0))*(mXm(0,2)-mXm(3,2));
                        mNormal( 2 ) = (mXm(0,0)-mXm(1,0))*(mXm(0,1)-mXm(3,1))-(mXm(0,1)-mXm(1,1))*(mXm(0,0)-mXm(3,0));

                        break ;
                    }
                    case 1 :
                    {
                        mNormal( 0 ) = (mXm(1,1)-mXm(2,1))*(mXm(1,2)-mXm(3,2))-(mXm(1,2)-mXm(2,2))*(mXm(1,1)-mXm(3,1));
                        mNormal( 1 ) = (mXm(1,2)-mXm(2,2))*(mXm(1,0)-mXm(3,0))-(mXm(1,0)-mXm(2,0))*(mXm(1,2)-mXm(3,2));
                        mNormal( 2 ) = (mXm(1,0)-mXm(2,0))*(mXm(1,1)-mXm(3,1))-(mXm(1,1)-mXm(2,1))*(mXm(1,0)-mXm(3,0));
                        break ;
                    }
                    case 2 :
                    {
                        mNormal( 0 ) = (mXm(0,2)-mXm(2,2))*(mXm(2,1)-mXm(3,1))-(mXm(0,1)-mXm(2,1))*(mXm(2,2)-mXm(3,2));
                        mNormal( 1 ) = (mXm(0,0)-mXm(2,0))*(mXm(2,2)-mXm(3,2))-(mXm(0,2)-mXm(2,2))*(mXm(2,0)-mXm(3,0));
                        mNormal( 2 ) = (mXm(0,1)-mXm(2,1))*(mXm(2,0)-mXm(3,0))-(mXm(0,0)-mXm(2,0))*(mXm(2,1)-mXm(3,1));
                        break ;
                    }
                    case 3 :
                    {
                        mNormal( 0 ) = (mXm(0,2)-mXm(1,2))*(mXm(0,1)-mXm(2,1))-(mXm(0,1)-mXm(1,1))*(mXm(0,2)-mXm(2,2));
                        mNormal( 1 ) = (mXm(0,0)-mXm(1,0))*(mXm(0,2)-mXm(2,2))-(mXm(0,2)-mXm(1,2))*(mXm(0,0)-mXm(2,0));
                        mNormal( 2 ) = (mXm(0,1)-mXm(1,1))*(mXm(0,0)-mXm(2,0))-(mXm(0,0)-mXm(1,0))*(mXm(0,1)-mXm(2,1));

                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "invalid face index ");
                    }
                }

                // This is a bit confusing: the sum of the weights for a triangle is 0.5,
                // so this should be multiplied by 2. On the other hand, however,
                // this is twice the surface of the triangle since we span a parallelogram.
                // Eventually, we multiply by 0.5 * 2 = 1
                mSurfaceIncrement = norm( mNormal );
                mNormal /= mSurfaceIncrement ;
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_tet_curved( const uint aIndex )
        {
            if( ! mIsCurved )
            {
                return this->normal_tet_straight( aIndex );
            }
            else if ( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix from the master element
                const Matrix< real > & J = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case 0 :
                    {
                        mNormal( 0 ) = J( 0, 1 ) * J( 1, 2 )-J( 0, 2 ) * J( 1, 1 );
                        mNormal( 1 ) = J( 0, 2 ) * J( 1, 0 )-J( 0, 0 ) * J( 1, 2 );
                        mNormal( 2 ) = J( 0, 0 ) * J( 1, 1 )-J( 0, 1 ) * J( 1, 0 );
                        break ;
                    }
                    case 1 :
                    {
                        mNormal( 0 ) = J( 1, 1 ) * J( 2, 2 )-J( 1, 2 ) * J( 2, 1 );
                        mNormal( 1 ) = J( 1, 2 ) * J( 2, 0 )-J( 1, 0 ) * J( 2, 2 );
                        mNormal( 2 ) = J( 1, 0 ) * J( 2, 1 )-J( 1, 1 ) * J( 2, 0 );
                        break ;
                    }
                    case 2 :
                    {
                        mNormal( 0 ) = J( 0, 2 ) * J( 2, 1 )-J( 0, 1 ) * J( 2, 2 );
                        mNormal( 1 ) = J( 0, 0 ) * J( 2, 2 )-J( 0, 2 ) * J( 2, 0 );
                        mNormal( 2 ) = J( 0, 1 ) * J( 2, 0 )-J( 0, 0 ) * J( 2, 1 );
                        break ;
                    }
                    case 3 :
                    {
                        mNormal( 0 ) = ( J( 2, 1 )-J( 0, 1 ) ) * ( J( 1, 2 )-J( 2, 2 ) )+( J( 0, 2 )-J( 2, 2 ) ) * ( J( 1, 1 )-J( 2, 1 ) );
                        mNormal( 1 ) = ( J( 0, 0 )-J( 2, 0 ) ) * ( J( 1, 2 )-J( 2, 2 ) )-( J( 0, 2 )-J( 2, 2 ) ) * ( J( 1, 0 )-J( 2, 0 ) );
                        mNormal( 2 ) = ( J( 2, 0 )-J( 0, 0 ) ) * ( J( 1, 1 )-J( 2, 1 ) )+( J( 0, 1 )-J( 2, 1 ) ) * ( J( 1, 0 )-J( 2, 0 ) );
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "invalid face index ");
                    }
                }

                // this gives us 2x the surface as with normal_tet_straight
                // again, remember that the sum of all weights on a triangle is 0.5, so we are good!
                mSurfaceIncrement = norm( mNormal );

                mNormal /= mSurfaceIncrement ;

            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_penta( const uint aIndex )
        {
            if ( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix for the master element
                const Matrix< real > & J = this->Jm( aIndex );

                Vector< real > & n = mNormal;

                switch ( mMasterIndex )
                {
                    case 0  :
                    {
                        n( 0 ) = J( 2, 1 ) * ( J( 0, 2 )-J( 1, 2 ) )-J( 2, 2 ) * ( J( 0, 1 )-J( 1, 1 ) );
                        n( 1 ) = J( 2, 2 ) * ( J( 0, 0 )-J( 1, 0 ) )-J( 2, 0 ) * ( J( 0, 2 )-J( 1, 2 ) );
                        n( 2 ) = J( 2, 0 ) * ( J( 0, 1 )-J( 1, 1 ) )-J( 2, 1 ) * ( J( 0, 0 )-J( 1, 0 ) );

                        // norm corresponds to surface
                        mSurfaceIncrement = norm( mNormal );
                        mNormal /= mSurfaceIncrement ;

                        // using quadrilateral weights, Σ w = 4, hence divide by four
                        mSurfaceIncrement *= 0.25;

                        break ;
                    }
                    case 1  :
                    {
                        n( 0 ) = J( 1, 2 ) * J( 2, 1 )-J( 1, 1 ) * J( 2, 2 );
                        n( 1 ) = J( 1, 0 ) * J( 2, 2 )-J( 1, 2 ) * J( 2, 0 );
                        n( 2 ) = J( 1, 1 ) * J( 2, 0 )-J( 1, 0 ) * J( 2, 1 );

                        // norm corresponds to surface
                        mSurfaceIncrement = norm( mNormal );
                        mNormal /= mSurfaceIncrement ;

                        // using quadrilateral weights, Σ w = 4, hence divide by four
                        mSurfaceIncrement *= 0.25;

                        break ;
                    }
                    case 2 :
                    {
                        n( 0 ) = J( 0, 1 ) * J( 2, 2 )-J( 0, 2 ) * J( 2, 1 );
                        n( 1 ) = J( 0, 2 ) * J( 2, 0 )-J( 0, 0 ) * J( 2, 2 );
                        n( 2 ) = J( 0, 0 ) * J( 2, 1 )-J( 0, 1 ) * J( 2, 0 );

                        // norm corresponds to surface
                        mSurfaceIncrement = norm( mNormal );
                        mNormal /= mSurfaceIncrement ;

                        // using quadrilateral weights, Σ w = 4, hence divide by four
                        mSurfaceIncrement *= 0.25;

                        break ;
                    }
                    case 3 :
                    {
                        n( 0 ) = J( 0, 2 ) * J( 1, 1 )-J( 0, 1 ) * J( 1, 2 );
                        n( 1 ) = J( 0, 0 ) * J( 1, 2 )-J( 0, 2 ) * J( 1, 0 );
                        n( 2 ) = J( 0, 1 ) * J( 1, 0 )-J( 0, 0 ) * J( 1, 1 );

                        // norm corresponds to two times the surface
                        mSurfaceIncrement = norm( mNormal );
                        mNormal /= mSurfaceIncrement ;

                        // using triangular weights, Σ w = 0.5, hence no further multiplication needed


                        break ;
                    }
                    case 4:
                    {
                        n( 0 ) = J( 0, 1 ) * J( 1, 2 )-J( 0, 2 ) * J( 1, 1 );
                        n( 1 ) = J( 0, 2 ) * J( 1, 0 )-J( 0, 0 ) * J( 1, 2 );
                        n( 2 ) = J( 0, 0 ) * J( 1, 1 )-J( 0, 1 ) * J( 1, 0 );


                        // norm corresponds to two times the surface
                        mSurfaceIncrement = norm( mNormal );
                        mNormal /= mSurfaceIncrement ;

                        // using triangular weights, Σ w = 0.5, hence no further multiplication needed

                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }
            }

            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::normal_hex( const uint aIndex )
        {
            if ( mNormalIndex != aIndex )
            {
                // remember the index
                mNormalIndex = aIndex;

                // compute the Jacobian matrix for the master element
                const Matrix< real > & J = this->Jm( aIndex );

                switch ( mMasterIndex )
                {
                    case( 0 ) :
                    {
                        // face eta=-1: n = cross(dX/dξ, dX/dζ)
                        mNormal( 0 ) = J( 0, 1 ) * J( 2, 2 ) - J( 0, 2 ) * J( 2, 1 );
                        mNormal( 1 ) = J( 0, 2 ) * J( 2, 0 ) - J( 0, 0 ) * J( 2, 2 );
                        mNormal( 2 ) = J( 0, 0 ) * J( 2, 1 ) - J( 0, 1 ) * J( 2, 0 );
                        break ;
                    }
                    case( 1 ) :
                    {
                        // face xi=+1: n = cross(dX/dη, dX/dζ)
                        mNormal( 0 ) = J( 1, 1 ) * J( 2, 2 ) - J( 1, 2 ) * J( 2, 1 );
                        mNormal( 1 ) = J( 1, 2 ) * J( 2, 0 ) - J( 1, 0 ) * J( 2, 2 );
                        mNormal( 2 ) = J( 1, 0 ) * J( 2, 1 ) - J( 1, 1 ) * J( 2, 0 );
                        break ;
                    }
                    case( 2 ) :
                    {
                        // face eta=+1: n = cross(dX/dζ, dX/dξ)
                        mNormal( 0 ) = J( 0, 2 ) * J( 2, 1 ) - J( 0, 1 ) * J( 2, 2 );
                        mNormal( 1 ) = J( 0, 0 ) * J( 2, 2 ) - J( 0, 2 ) * J( 2, 0 );
                        mNormal( 2 ) = J( 0, 1 ) * J( 2, 0 ) - J( 0, 0 ) * J( 2, 1 );
                        break ;
                    }
                    case( 3 ) :
                    {
                        // face xi=-1: n = cross(dX/dζ, dX/dη)
                        mNormal( 0 ) = J( 2, 1 ) * J( 1, 2 ) - J( 2, 2 ) * J( 1, 1 );
                        mNormal( 1 ) = J( 2, 2 ) * J( 1, 0 ) - J( 2, 0 ) * J( 1, 2 );
                        mNormal( 2 ) = J( 2, 0 ) * J( 1, 1 ) - J( 2, 1 ) * J( 1, 0 );
                        break ;
                    }
                    case( 4 ) :
                    {
                        // face zeta=-1: n = cross(dX/dη, dX/dξ)
                        mNormal( 0 ) = J( 1, 1 ) * J( 0, 2 ) - J( 1, 2 ) * J( 0, 1 );
                        mNormal( 1 ) = J( 1, 2 ) * J( 0, 0 ) - J( 1, 0 ) * J( 0, 2 );
                        mNormal( 2 ) = J( 1, 0 ) * J( 0, 1 ) - J( 1, 1 ) * J( 0, 0 );
                        break ;
                    }
                    case( 5 ) :
                    {
                        // face zeta=+1: n = cross(dX/dξ, dX/dη)
                        mNormal( 0 ) = J( 0, 1 ) * J( 1, 2 ) - J( 0, 2 ) * J( 1, 1 );
                        mNormal( 1 ) = J( 0, 2 ) * J( 1, 0 ) - J( 0, 0 ) * J( 1, 2 );
                        mNormal( 2 ) = J( 0, 0 ) * J( 1, 1 ) - J( 0, 1 ) * J( 1, 0 );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid master index for facet" );
                    }
                }

                // this gives us 0.25 * the surface, but the sum of all weights on a quad is 4, so we are good here!
                mSurfaceIncrement = norm( mNormal ) ;

                mNormal /= mSurfaceIncrement ;
            }

            // now we can return the vector
            return mNormal;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::node_data( const string & aNodeField )
        {
            Vector< real > & aVector = mVectorMap( aNodeField )->vector();
            mGroup->parent()->iwg()->collect_node_data( mElement, aNodeField, aVector );
            return aVector ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::qold( const uint aStep )
        {
            for( uint k=0; k<mElement->number_of_local_dofs(); ++k )
            {
                Dof * tDof = mElement->local_dof( k ) ;

                // a too-small table means the timestepping order changed after
                // allocate() without an init_qold_table() refresh
                BELFEM_ASSERT( aStep * ( mMaxDofFieldIndex + 1 )
                        + tDof->field_index() < mQold.size(),
                    "qold table too small for step %u ( stale timestepping order? )",
                    ( unsigned int ) aStep );

                Vector< real > * tField
                    = mQold( aStep * ( mMaxDofFieldIndex + 1 ) + tDof->field_index() );

                BELFEM_ASSERT( tField != nullptr,
                    "no qold entry for step %u, field %lu",
                    ( unsigned int ) aStep,
                    ( long unsigned int ) tDof->field_index() );

                mq0( k ) = ( *tField )( tDof->dof_index_on_field() );
            }

            return mq0 ;
        }

        void
        Calculator::init_qold_table()
        {
            mMaxDofFieldIndex = 0 ;
            mQold.clear();

            // nothing to do if this group has no kernel or no timestepping scheme
            if ( mGroup->parent() == nullptr )
            {
                return ;
            }

            uint tOrder = mGroup->parent()->iwg()->timestepping_order() ;

            if ( tOrder == 0 )
            {
                return ;
            }

            // get the list of dof strings
            const Cell< string > & tFields = mGroup->parent()->iwg()->dof_fields();

            // loop over all fields
            for ( const string & tLabel : tFields )
            {
                mesh::Field * tField = mMesh->field( tLabel );

                // make sure that the fields exist
                for ( uint s=0; s<tOrder; ++s )
                {
                    string tOldLabel = tLabel + std::to_string( s );

                    if ( ! mMesh->field_exists( tOldLabel ) )
                    {
                        mMesh->create_field( tOldLabel, tField->entity_type() );
                    }
                }
                mMaxDofFieldIndex = std::max( mMaxDofFieldIndex, tField->index() ) ;
            }

            // stride is the inclusive maximum plus one, otherwise
            // ( s, max ) and ( s+1, 0 ) would alias
            mQold.set_size( tOrder * ( mMaxDofFieldIndex + 1 ), nullptr );

            for ( uint s=0; s<tOrder; ++s )
            {
                // loop over all fields
                for ( const string & tLabel : tFields )
                {
                    mesh::Field * tField = mMesh->field( tLabel );

                    string tOldLabel = tLabel + std::to_string( s );

                    mesh::Field * tOldField = mMesh->field( tOldLabel );
                    tOldField->set_write_to_file_flag( false );
                    mQold( s * ( mMaxDofFieldIndex + 1 ) + tField->index() ) = & tOldField->data();
                }
            }
        }

//------------------------------------------------------------------------------

        /**
         * return the swap vector for the dofs
         */
        Vector< real > &
        Calculator::qswap()
        {
            return mqswap ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        Calculator::q()
        {
            for( uint k=0; k<mElement->number_of_local_dofs(); ++k )
            {
                Dof * tDof = mElement->local_dof( k );
                mq( k ) = mMesh->field( tDof->field_index() )->data()( tDof->dof_index_on_field() );
            }
            return mq ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::E( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunction != nullptr,
                "Edge function has not been assigned for group %lu",
                ( long unsigned int ) mGroup->id() );

                return mEdgeFunction->E( aIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::C( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunction != nullptr,
                "Edge function has not been assigned for group %lu",
                ( long unsigned int ) mGroup->id() );

            return mEdgeFunction->C( aIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::G( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunction != nullptr,
                "Edge function has not been assigned for group %lu",
                ( long unsigned int ) mGroup->id() );

            return mEdgeFunction->G( aIndex );
        }

//------------------------------------------------------------------------------

        uint
        Calculator::num_nedelec_dofs() const
        {
            if ( mEdgeFunctionMaster != nullptr )
            {
                return mEdgeFunctionMaster->ndofs() ;
            }
            else if( mEdgeFunction != nullptr )
            {
                return  mEdgeFunction->ndofs() ;
            }
            else
            {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::Em( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunctionMaster != nullptr,
                "Master Edge function has not been assigned" );

            return mEdgeFunctionMaster->E( aIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::Cm( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunctionMaster != nullptr,
                "Master Edge function has not been assigned" );

            return mEdgeFunctionMaster->C( aIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::Es( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunctionSlave != nullptr,
                "Slave Edge function has not been assigned" );

            return mEdgeFunctionSlave->E( aIndex );
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Calculator::Cs( const uint aIndex )
        {
            BELFEM_ASSERT( mEdgeFunctionSlave != nullptr,
                "Slave Edge function has not been assigned" );

            return mEdgeFunctionSlave->C( aIndex );
        }

//------------------------------------------------------------------------------

        const Material *
        Calculator::material() const
        {
            BELFEM_ASSERT( mGroup->material() != nullptr,
                "no material assigned for %s %lu",
                mGroup->type() == GroupType::BLOCK ? "block" : "sideset",
                ( long unsigned int ) mGroup->id() );

            return mGroup->material();
        }

//------------------------------------------------------------------------------

        void
        Calculator::link_element_maxwell( Element * aElement )
        {
            mMaxwellData->reset();

            mElement = aElement;

        }

        void
        Calculator::link_element_maxwell_thermal( Element * aElement )
        {
            mMaxwellData->reset();

            mElement = aElement;

            Calculator * tPeer = mMaxwellData->thermal() ;

            // relink the thermal peer unless it already sits on this element
            // ( element() is null before the peer's very first link ). an aura
            // element may be missing from the thermal group, whose aura is
            // expanded from a smaller owned-element set; the peer then keeps
            // its previous element ( postprocessing only, never assembly )
            //
            // KEEPING THE PEER IS DELIBERATE. Do not turn the miss into an
            // assert or an error: this is the sanctioned fallback for a maxwell
            // element with no thermal counterpart, and coupled models added
            // later may need it. The obligation sits on the CONSUMER, which
            // must not read a peer it did not get -- the superconductor postproc
            // compares the peer's element id against its own and uses gTbulk
            // when they differ. Measured on 2 ranks over a bulk
            // superconductor deck: the miss did not occur in 1.48 M peer reads
            // per pass, so the branch is latent rather than dead
            if ( ( tPeer->element() == nullptr || tPeer->element()->id() != aElement->id() )
                 && tPeer->group()->element_exists( aElement->id() ) )
            {
                tPeer->link( tPeer->group()->element( aElement->id() ) );
            }

            BELFEM_ASSERT( this->num_intpoints() == tPeer->num_intpoints(),
                "Number of integration points does not match" );
        }

        void
        Calculator::link_element_thermal_maxwell( Element * aElement )
        {
            mMaxwellData->reset();

            mElement = aElement;

            Calculator * tPeer = mMaxwellData->maxwell() ;

            // relink the maxwell peer unless it already sits on this element
            // ( element() is null before the peer's very first link )
            if ( tPeer->element() == nullptr || tPeer->element()->id() != aElement->id() )
            {
                tPeer->link( tPeer->group()->element( aElement->id() ) );
            }

            BELFEM_ASSERT( this->num_intpoints() == tPeer->num_intpoints(),
                "Number of integration points does not match" );
        }
    }
}
