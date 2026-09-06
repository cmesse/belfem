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
#ifndef BELFEM_CL_FEM_CALCULATOR_HPP
#define BELFEM_CL_FEM_CALCULATOR_HPP
#include <algorithm>
#include <cmath>

#include "cl_Material.hpp"
#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"

#include "cl_Mesh.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "en_IWGs.hpp"
#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        class Element ;
        class EdgeFunction ;
        class Group ;
        class Calculator ;
        class Kernel ;

//------------------------------------------------------------------------------

        namespace calculator
        {

            class VectorData
            {
                const string     mLabel ;
                const EntityType mType  ;
                uint       mIndex = BELFEM_UINT_MAX ;

                Vector< real > mVectorData ;
            public:

                VectorData( const string & aLabel, const uint aSize, const EntityType aType );

                ~VectorData() = default ;

                const string &
                label() const ;

                void
                set_index( const uint aIndex );

                uint
                index() const ;

                Vector< real > &
                vector();

                EntityType
                entity_type() const ;

            };

            class MatrixData
            {
                const string     mLabel ;
                uint       mIndex = BELFEM_UINT_MAX ;
                Matrix< real > mMatrixData ;

            public:

                MatrixData( const string & aLabel, const uint aNumRows, const uint aNumCols );


                ~MatrixData() = default ;

                const string &
                label() const ;

                void
                set_index( const uint aIndex );

                uint
                index() const ;

                Matrix< real > &
                matrix();

            };

            enum class MaxwellDataValue
            {
                T         =  0,
                H         =  1,
                B         =  2,
                j         =  3,
                rho       =  4,
                cp        =  5,
                lambda    =  6,
                beta      =  7,
                normH     =  8,
                normB     =  9,
                normJ     = 10,
                n         = 11,  // <-- normal
                drhodj    = 12,
                drhodT    = 13,
                dcpdT     = 14,
                dlambdadT = 15,
                x         = 16,
                mu        = 17,
                dmudh     = 18,
                drhodb    = 19,
                drhodbeta = 20,
                UNDEFINED = 21
            };

            class MaxwellData
            {
                Calculator * mMaxwellCalculator  = nullptr ;
                Calculator * mThermalCalculator = nullptr ;

                Material   * mMaterial = nullptr ;

                Matrix< real >   mCoords ;
                real mX = 0.0 ;
                real mY = 0.0 ;
                real mZ = 0.0 ;

                const real & mTime ;

                Vector< real > & mH ;
                Vector< real > & mHn ;
                Vector< real > & mHt ;

                Vector< real > & mB ;
                Vector< real > & mBn ;
                Vector< real > & mBt ;

                Vector< real > & mJ ;
                Vector< real > & mN ;

                real mRho    = BELFEM_QUIET_NAN ;
                real mdRhodJ = BELFEM_QUIET_NAN ;
                real mdRhodT = BELFEM_QUIET_NAN ;
                real mdRhodB = BELFEM_QUIET_NAN ;
                real mdRhodBeta = BELFEM_QUIET_NAN ;
                real mBeta   = 0.5 * constant::pi ;

                real mT = BELFEM_QUIET_NAN ;

                real mCp = BELFEM_QUIET_NAN ;
                real mdCpdT = BELFEM_QUIET_NAN ;
                real mLambda = BELFEM_QUIET_NAN ;
                real mdLambdadT = BELFEM_QUIET_NAN ;

                real mNormH = BELFEM_QUIET_NAN ;
                real mNormB = BELFEM_QUIET_NAN ;
                real mNormJ = BELFEM_QUIET_NAN ;

                real mMu       = constant::mu0 ;
                real mdMudH    = 0 ;

                //! true while compute_T had to clamp the FEM temperature into
                //! the physical table window [ gTmin, mTmax ]. While clamped,
                //! all dT-derivatives are zero ( consistent tangent )
                bool mTClamped   = false ;

                //! true while compute_rho had to clamp into [ gRhoMin, gRhoMax ].
                //! While clamped, drho/dJ, drho/dT, drho/dB and drho/dbeta are zero
                bool mRhoClamped = false ;

                //! upper temperature bound for the clamp: the material's T_max
                //! if defined, otherwise unbounded
                real mTmax = BELFEM_REAL_MAX ;

                //! density at the undeformed-mesh reference temperature gTroom.
                //! PHYSICS TRAP: never evaluate density( T ) per point — the
                //! computation runs on the undeformed mesh and the transport
                //! properties are already corrected for thermal expansion, so a
                //! temperature-dependent density would double-count the expansion
                real mDensity = BELFEM_QUIET_NAN ;

                real ( MaxwellData::*mFunT )( const uint aIndex ) = nullptr ;


                real ( MaxwellData::*mFunRho )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundRhodT )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundRhodJ )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundRhodB )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundRhodBeta )( const uint aIndex ) = nullptr ;

                //! artificial volumetric heat load: compute_heatload_user when
                //! the material carries a heating plugin, return_zero otherwise
                real ( MaxwellData::*mFunHeat )( const uint aIndex ) = nullptr ;

                real ( MaxwellData::*mFunMu )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundMudH )( const uint aIndex ) = nullptr ;

                real ( MaxwellData::*mFunLambda )( const uint aIndex ) = nullptr ;
                real ( MaxwellData::*mFundLambdadT )( const uint aIndex ) = nullptr ;

                void ( MaxwellData::*mFunX )( const uint aIndex ) = nullptr ;
                const Vector< real > & ( MaxwellData::*mFunH )( const uint aIndex ) = nullptr ;
                const Vector< real > & ( MaxwellData::*mFunB )( const uint aIndex ) = nullptr ;

                Vector< uint > mLastIndex ;

                //! side connector state ( see compute_h_side_connector ):
                //! per-element master link and frame, rebuilt lazily after
                //! reset() clears mFrameCurrent
                Calculator * mReferenceCalc = nullptr ;
                bool mFrameCurrent = false ;
                //! true while a thermal kernel exists: compute_T_side_connector
                //! interpolates the seam temperatures, otherwise gTbulk
                bool mHaveSeamT = false ;
                //! set by prepare_side_connector_frame(), point into the
                //! "binomial" and "Tseam" workspaces of the wall calculator
                Vector< real > * mBinomialVec = nullptr ;
                Vector< real > * mTseamVec = nullptr ;
                //! scratch for the master in-plane field ( Em * q_master )
                Vector< real > mWork ;

            public:

                MaxwellData( Calculator * aCalculator,
                             Kernel * aMaxwellKernel,
                             Kernel * aThermalKernel ) ;

                ~MaxwellData() = default ;

                const Vector< real > &
                compute_h( const uint aIndex ) ;

                const Vector< real > &
                compute_b( const uint aIndex ) ;

                const Vector< real > &
                compute_j( const uint aIndex ) ;

                void
                compute_x( const uint aIndex ) ;

                real
                compute_rho( const uint aIndex ) ;

                real
                compute_drhodT( const uint aIndex ) ;

                real
                compute_drhodj( const uint aIndex ) ;

                real
                compute_drhodb( const uint aIndex ) ;

                real
                compute_drhodbeta( const uint aIndex ) ;

                real
                compute_T( const uint aIndex ) ;

                real
                compute_cp( const uint aIndex ) ;

                real
                compute_dcpdT( const uint aIndex ) ;

                real
                compute_lambda( const uint aIndex ) ;

                real
                compute_dlambdadT( const uint aIndex ) ;

                void
                reset();

                real
                norm_b( const uint aIndex );

                real
                norm_j( const uint aIndex );

                real
                density() const ;

                bool
                T_clamped() const ;

                bool
                rho_clamped() const ;

                Calculator *
                maxwell() ;

                Calculator *
                thermal() ;

                real
                compute_mu( const uint aIndex ) ;

                real
                compute_dmudh( const uint aIndex ) ;

                /**
                 * artificial volumetric heat load [ W/m³ ] at integration point
                 * aIndex, from the material's heating plugin; zero without one
                 */
                real
                compute_volumetric_heatload( const uint aIndex ) ;

            private:

                Vector< real > &
                link_vector( const string & aLabel );

                bool
                is_current( MaxwellDataValue aValue, const uint aIndex ) const ;

                void
                compute_x_2d( const uint aIndex );

                void
                compute_x_3d( const uint aIndex );

                void
                set( MaxwellDataValue aValue, const uint aIndex );

                real
                compute_T_fem( const uint aIndex ) ;

                real
                compute_T_const( const uint aIndex ) ;

                //! seam-node temperature of the edge-coating wall ( gTbulk
                //! fallback on magnetic-only runs ), compute_T_fem clamp contract
                real
                compute_T_side_connector( const uint aIndex ) ;

                const Vector< real > &
                compute_h_bulk_edge( const uint aIndex ) ;

                const Vector< real > &
                compute_h_ts_edge( const uint aIndex ) ;

                const Vector< real > &
                compute_h_bulk_node( const uint aIndex ) ;

                //! h = ht + hb + hn recovery of the edge-coating wall element
                const Vector< real > &
                compute_h_side_connector( const uint aIndex ) ;

                //! once per element: link the master layer-block calculator,
                //! copy hn / normal, build tangent + binomial, fetch Tseam
                void
                prepare_side_connector_frame() ;

                //! map the wall's psi coordinate onto the master's ( xi, eta )
                void
                side_connector_xi_eta( const real aPsi, real & aXi, real & aEta ) const ;

                real
                compute_mu_0( const uint aIndex );

                real
                compute_mu_const( const uint aIndex );

                real
                compute_mu_h( const uint aIndex );

                const Vector< real > &
                compute_b_bulk( const uint aIndex ) ;

                const Vector< real > &
                compute_b_ts( const uint aIndex ) ;

                real
                compute_lambda_bulk( const uint aIndex );

                real
                compute_dlambdadT_bulk( const uint aIndex );

                real
                compute_lambda_metal( const uint aIndex );

                real
                compute_dlambdadT_metal( const uint aIndex );


                real
                compute_rho_bulk( const uint aIndex );

                real
                compute_drhodT_bulk( const uint aIndex );

                real
                compute_rho_metal( const uint aIndex );

                real
                compute_drhodT_metal( const uint aIndex );

                real
                compute_drhodb_metal( const uint aIndex );

                real
                compute_drhodbeta_metal( const uint aIndex );

                //! HTS drho/d|B| family: jc(T,|B|,θ) / n(T,|B|,θ)
                //! from the lookup table. NOTE deliberately no drhodbeta
                //! twin: the HTS angle is bn_angle ( field to tape normal )
                //! while add_rho_field_tangent differentiates bj_angle
                //! ( field to current, the metal Kohler variable ) — binding
                //! the beta channel would apply the wrong ∂β/∂q rows
                //! ( 2026-08-13 audit, both voices, C1 refuted for β )
                real
                compute_drhodb_powerlaw_ts( const uint aIndex );

                real
                compute_drhodb_powerlaw_ts_defect( const uint aIndex );

                real
                compute_drhodb_piecewise_ts( const uint aIndex );

                real
                compute_drhodb_piecewise_ts_defect( const uint aIndex );

                real
                compute_drhodb_powerlaw_bulk( const uint aIndex );

                real
                compute_drhodb_powerlaw_bulk_defect( const uint aIndex );

                real
                compute_drhodb_piecewise_bulk( const uint aIndex );

                real
                compute_drhodb_piecewise_bulk_defect( const uint aIndex );

                //! HTS drho/dT family ( T-leg, 2026-08-13 ): the
                //! quench-feedback tangent, consumed by T_h_newton only.
                //! Replaces the retired compute_drhodT_hts, whose
                //! b = a·c/(c−a) reconstruction was sign-flipped and
                //! diverged at the flux-flow crossover — the closed forms
                //! live in Material::drho_{powerlaw,piecewise}_dT
                real
                compute_drhodT_powerlaw_ts( const uint aIndex );

                real
                compute_drhodT_powerlaw_ts_defect( const uint aIndex );

                real
                compute_drhodT_piecewise_ts( const uint aIndex );

                real
                compute_drhodT_piecewise_ts_defect( const uint aIndex );

                real
                compute_drhodT_powerlaw_bulk( const uint aIndex );

                real
                compute_drhodT_powerlaw_bulk_defect( const uint aIndex );

                real
                compute_drhodT_piecewise_bulk( const uint aIndex );

                real
                compute_drhodT_piecewise_bulk_defect( const uint aIndex );

                //! riva law ( 2026-08-27 ): the Duron parallel model made
                //! total over the full jc/n table range, see powerlaws.hpp
                real
                compute_rho_riva_ts( const uint aIndex );

                real
                compute_rho_riva_ts_defect( const uint aIndex );

                real
                compute_rho_riva_bulk( const uint aIndex );

                real
                compute_rho_riva_bulk_defect( const uint aIndex );

                real
                compute_drhodj_riva_ts( const uint aIndex );

                real
                compute_drhodj_riva_ts_defect( const uint aIndex );

                real
                compute_drhodj_riva_bulk( const uint aIndex );

                real
                compute_drhodj_riva_bulk_defect( const uint aIndex );

                real
                compute_drhodb_riva_ts( const uint aIndex );

                real
                compute_drhodb_riva_ts_defect( const uint aIndex );

                real
                compute_drhodb_riva_bulk( const uint aIndex );

                real
                compute_drhodb_riva_bulk_defect( const uint aIndex );

                real
                compute_drhodT_riva_ts( const uint aIndex );

                real
                compute_drhodT_riva_ts_defect( const uint aIndex );

                real
                compute_drhodT_riva_bulk( const uint aIndex );

                real
                compute_drhodT_riva_bulk_defect( const uint aIndex );

                real
                compute_rho_powerlaw_bulk( const uint aIndex );

                real
                compute_drhodj_powerlaw_bulk( const uint aIndex );

                real
                compute_rho_powerlaw_ts( const uint aIndex );

                real
                compute_drhodj_powerlaw_ts( const uint aIndex );

                real
                compute_rho_piecewise_bulk( const uint aIndex );

                real
                compute_drhodj_piecewise_bulk( const uint aIndex );

                real
                compute_rho_piecewise_ts( const uint aIndex );

                real
                compute_drhodj_piecewise_ts( const uint aIndex );

                real
                compute_rho_powerlaw_bulk_defect( const uint aIndex );

                real
                compute_drhodj_powerlaw_bulk_defect( const uint aIndex );

                real
                compute_rho_powerlaw_ts_defect( const uint aIndex );

                real
                compute_drhodj_powerlaw_ts_defect( const uint aIndex );

                real
                compute_rho_piecewise_bulk_defect( const uint aIndex );

                real
                compute_drhodj_piecewise_bulk_defect( const uint aIndex );

                real
                compute_rho_piecewise_ts_defect( const uint aIndex );

                real
                compute_drhodj_piecewise_ts_defect( const uint aIndex );

                real
                return_zero( const uint aIndex );

                real
                beta_dummy() const ;

                real
                compute_dmu_zero( const uint aIndex  ) ;

                real
                compute_dmu_material( const uint aIndex  ) ;

                real
                compute_heatload_user( const uint aIndex ) ;

            };
        } // end namespace calculator

//------------------------------------------------------------------------------

        class Calculator
        {
            // link to group
            Group * mGroup = nullptr ;

            // link to mesh
            Mesh  * mMesh = nullptr ;

            const ModelDimensionality mDimensionality ;

            //! timsetep from equation class
            const real & mTimestep ;

            // link to current element
            Element * mElement = nullptr ;

            EdgeFunction * mEdgeFunction       = nullptr ;
            EdgeFunction * mEdgeFunctionMaster = nullptr ;
            EdgeFunction * mEdgeFunctionSlave  = nullptr ;

            Cell< EdgeFunction * > mEdgeFunctionsMaster ;
            Cell< EdgeFunction * > mEdgeFunctionsSlave ;

            // switch telling if we are allocated
            bool mIsAllocated = false ;

            uint mNumberOfNodes = BELFEM_UINT_MAX ;
            uint mNumberOfCornerNodes = BELFEM_UINT_MAX ;
            uint mNumberOfIntegrationPoints = 0 ;

            uint mNumberOfNodesOnMaster = BELFEM_UINT_MAX ;
            uint mNumberOfNodesOnSlave = BELFEM_UINT_MAX ;

            uint mNumberOfEdgesOnMaster = BELFEM_UINT_MAX ;
            uint mNumberOfEdgesOnSlave = BELFEM_UINT_MAX ;

            uint mNumberOfFacesOnMaster = BELFEM_UINT_MAX ;
            uint mNumberOfFacesOnSlave  = BELFEM_UINT_MAX ;
            uint mMasterIndex = BELFEM_UINT_MAX ;

            IntegrationData * mDomainIntegration    = nullptr ;
            IntegrationData * mLinearIntegration    = nullptr ;
            IntegrationData * mThinShellIntegration = nullptr ;
            const IntegrationData * mMasterIntegration    = nullptr ;
            const IntegrationData * mSlaveIntegration     = nullptr ;

            // for enrichment
            const IntegrationData * mMasterVolumeIntegration = nullptr ;
            const IntegrationData * mVolumeEnrichment        = nullptr ;
            const IntegrationData * mSideSetEnrichment       = nullptr ;

            //! flag telling if element is linear
            bool mIsLinear = false ;

            //! flag telling if element is straight or curved
            bool mIsCurved = false ;

            //! stiffness matrix
            Matrix< real > mK ;

            //! mass matrix
            Matrix< real > mM ;

            //! Newton correction matrix
            Matrix< real > mJN ;

            //! load vector
            Vector< real > mf ;

            //! dof vector for last timestep
            Vector< real > mq0 ;

            //! dof vector for next timestep
            Vector< real > mq ;

            //! dof vector for swapping
            Vector< real > mqswap ;

            //! normal vector
            Vector< real > mNormal ;
            uint mNormalIndex      = BELFEM_UINT_MAX ;
            real mSurfaceIncrement = BELFEM_QUIET_NAN ;

            //! dof labels for next timestep
            Cell< string >     mDofLabels ;

            //! custom matrices
            Cell< calculator::VectorData * > mVectors ;
            Cell< calculator::MatrixData * > mMatrices ;

            //! maps
            Map< string , calculator::VectorData * > mVectorMap ;
            Map< string , calculator::MatrixData * > mMatrixMap ;

            //! model parameters, multi purpose
            Vector< real > mModelParameters ;

            // pointers for faster access
            Matrix< real > mX  ;  // node coordinates
            Matrix< real > mXc ;  // node coordinates at corners
            calculator::MatrixData * mJ    = nullptr ;  // jacobian matrix
            calculator::MatrixData * mInvJ = nullptr ;  // inverse of the jacobian matrix
            calculator::MatrixData * mN    = nullptr ;  // node interpolatoin operator
            calculator::MatrixData * mdN   = nullptr ;  // derivative function
            calculator::MatrixData * mB    = nullptr ;  // gradient operator

            // for faces
            Matrix< real > mXm ;
            calculator::MatrixData * mNm = nullptr ;  // node interpolation operator
            calculator::MatrixData * mJm = nullptr ;
            calculator::MatrixData * mInvJm = nullptr ;
            calculator::MatrixData * mBm = nullptr ;
            calculator::MatrixData * mntBm = nullptr ; // trans( normal ) * B
            calculator::MatrixData * mnxBm = nullptr ; // cross( normal, B )

            Matrix< real > mXs ;
            calculator::MatrixData * mNs = nullptr ;  // node interpolation operator
            calculator::MatrixData * mJs = nullptr ;
            calculator::MatrixData * mInvJs = nullptr ;
            calculator::MatrixData * mBs = nullptr ;
            calculator::MatrixData * mntBs = nullptr ; // trans( normal ) * B
            calculator::MatrixData * mnxBs = nullptr ; // cross( normal, B )

            calculator::MaxwellData * mMaxwellData = nullptr ;

            // kernels for the maxwell data helper, set via link_maxwell();
            // construction of the helper is deferred to allocate()
            Kernel * mMaxwellKernel = nullptr ;
            Kernel * mThermalKernel = nullptr ;

            real mDetJ      = BELFEM_QUIET_NAN ;
            uint mDetJIndex = BELFEM_UINT_MAX ;
            real mRadius    = BELFEM_QUIET_NAN ; // only needed if axisymmetric

            // function to compute the flattened slave integration index
            // (cumulative facet/orientation offset)
            uint ( Calculator::*mFunSlaveIntegrationIndex )( const mesh::Facet * aFacet ) ;

            real ( * mFunInvertJ )( const Matrix< real > & aA, Matrix< real > & aB );

            const Vector< real > & ( Calculator::*mFunNormal )( const uint aIndex );

            // function for node interpolator
            const Matrix< real > & (  Calculator::*mFunN )( const uint aIndex );

            // function for gradient operator
            const Matrix< real > & (  Calculator::*mFunB )( const uint aIndex );

            // function for node interpolator master
            const Matrix< real > & (  Calculator::*mFunNm )( const uint aIndex );

            // function for gradient operator master
            const Matrix< real > & (  Calculator::*mFunBm )( const uint aIndex );

            // function for node interpolator slave
            const Matrix< real > & (  Calculator::*mFunNs )( const uint aIndex );

            // function for gradient operator slave
            const Matrix< real > & (  Calculator::*mFunBs )( const uint aIndex );

            // volume increment
            real ( Calculator::*mFundV )( const uint aIndex );

            // surface increment
            real ( Calculator::*mFundS )( const uint aIndex );

            // inverse the jacobian
            const Matrix< real > & ( Calculator::*mFunInvJ )( const uint aIndex );

            const Vector< real > &
            ( Calculator::*mFunCollectNodeData )( const string & aLabel );

            real
            ( Calculator::*mFunBJAngle )(
                const Vector< real > & b,
                const Vector< real > & j,
                real & norm_b,
                real & norm_j ) const ;

            // for theta method
            real mTheta = 1.0 ;
            real mOneMinusTheta = 0.0 ;

            const Vector< real > & ( Calculator::*mFunNedelecDataH )();
            const Vector< real > & ( Calculator::*mFunNedelecDataA )();

            // defaulted so that no allocate() path can ever leave it null
            void ( Calculator::*mFunLinkElement )( Element * aElement )
                = & Calculator::link_element_default ;

            uint mIntegrationOrder = 0 ;

            Cell< Vector< real > * > mQold ;
            index_t mMaxDofFieldIndex = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            //! default constructor
            Calculator( Group * aGroup, const ModelDimensionality aDimensionality );

            // special constrictor used for TET10 T-Matrices
            Calculator( Group * aGroup, Mesh * aMesh );

//------------------------------------------------------------------------------

            ~Calculator() ;

//------------------------------------------------------------------------------

            void
            init_qold_table();

//------------------------------------------------------------------------------

            /**
             * called by dof manager
             */
            void
            allocate();

//------------------------------------------------------------------------------

            void
            link( Group * aGroup );

//------------------------------------------------------------------------------

            void
            link( Element * aElement );

//------------------------------------------------------------------------------

            void
            link( mesh::Facet * aFacet );

//------------------------------------------------------------------------------

            Group *
            group() ;

//------------------------------------------------------------------------------

            Element *
            element() ;

//------------------------------------------------------------------------------

            const Material *
            material() const;

//------------------------------------------------------------------------------

            /**
             * return the stiffness matrix
             */
            Matrix< real > &
            K() ;

//------------------------------------------------------------------------------

            /**
             * return the mass matrix
             */
            Matrix< real > &
            M() ;

//------------------------------------------------------------------------------

            /**
             * return the Newton correction matrix
             */
            Matrix< real > &
            JN() ;

//------------------------------------------------------------------------------

            /**
             * return the load vector
             */
            Vector< real > &
            f() ;

//------------------------------------------------------------------------------

            /**
             * return the dof vector at current timestep
             */
            const Vector< real > &
            q() ;

//------------------------------------------------------------------------------

            /**
             * return the dof vector at an old timestep
             */
            const Vector< real > &
            qold( const uint aStep=0 ) ;

//------------------------------------------------------------------------------

            /**
             * return a swap vector for the dofs
             */
            Vector< real > &
            qswap() ;

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            const Vector< real > &
            node_data( const string & aNodeField );

//------------------------------------------------------------------------------

            bool
            vector_exists( const string aLabel ) const ;

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            Vector< real > &
            vector( const string & aLabel );

//------------------------------------------------------------------------------

            bool
            matrix_exists( const string aLabel ) const ;

//------------------------------------------------------------------------------

            /**
             * return a matrix object
             */
            Matrix< real > &
            matrix( const string & aLabel );

//------------------------------------------------------------------------------

            /**
             * node coordinates on element
             */
            const Matrix< real > &
            X() const ;

//------------------------------------------------------------------------------

            /**
             * jacobian matrix
             */
            const Matrix< real > &
            J( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * inverse of jacobian matrix
             */
            const Matrix< real > &
            invJ( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * jacobian matrix for master
             */
            const Matrix< real > &
            Jm( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * jacobian matrix for slave
             */
            const Matrix< real > &
            Js( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * node interpolator
             */
            const Matrix< real > &
            N( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * edge interpolator
             */
            const Matrix< real > &
            E( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * curl interpolator
             */
            const Matrix< real > &
            C( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * gradient interpolator for nedelec elements
             */
            const Matrix< real > &
            G( const uint aIndex ) ;

//------------------------------------------------------------------------------

            uint
            num_nedelec_dofs() const ;

//------------------------------------------------------------------------------

            /**
             * edge interpolator on master
             */
            const Matrix< real > &
            Em( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * curl interpolator on master
             */
            const Matrix< real > &
            Cm( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * edge interpolator on slave
             */
            const Matrix< real > &
            Es( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * curl interpolator on slave
             */
            const Matrix< real > &
            Cs( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * node interpolator master
             */
            const Matrix< real > &
            Nm( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * node interpolator slave
             */
            const Matrix< real > &
            Ns( const uint aIndex ) ;

//------------------------------------------------------------------------------

            /**
             * node interpolator, but as vector
             */
            const Vector< real > &
            Nvec( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            /**
            * node interpolation function for local values
            */
            real
            node_interp( const uint aIndex, const Vector< real > & aNodeValues ) const;

//------------------------------------------------------------------------------

            /**
             * gradient operator
             */
            const Matrix< real > &
            B( const uint aIndex=0 ) ;

//------------------------------------------------------------------------------

            /**
             * gradient operator master
             */
            const Matrix< real > &
            Bm( const uint aIndex=0 ) ;

//------------------------------------------------------------------------------

            /**
             * gradient operator slave
             */
            const Matrix< real > &
            Bs( const uint aIndex=0 ) ;

//------------------------------------------------------------------------------

            /**
             * node coordinates on master element
             */
            const Matrix< real > &
            Xm() const ;

//------------------------------------------------------------------------------

            /**
             * node coordinates on slave element
             */
            const Matrix< real > &
            Xs() const ;

//------------------------------------------------------------------------------

            /**
             * surface increment
             */
            real
            dS ( const uint aIndex=0 ) ;

//------------------------------------------------------------------------------

            /**
             * volume increment
             */
            real
            dV ( const uint aIndex=0 ) ;

//------------------------------------------------------------------------------

            /**
             * returns the normal of a surface
             * @param aIndex
             * @return
             */
            const Vector< real > &
            normal( const uint aIndex=0 );

//------------------------------------------------------------------------------

            void
            initialize_integration( const ElementType    aElementType,
                                    const InterpolationType aInterpolationType );

//------------------------------------------------------------------------------

            void
            set_integration_order( const uint aOrder );

            uint
            integration_order() const ;

//------------------------------------------------------------------------------

            void
            allocate_memory();

//------------------------------------------------------------------------------

            const IntegrationData *
            integration() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            master_integration() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            slave_integration() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            volume_integration() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            volume_enrichment() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            sideset_enrichment() const ;

//------------------------------------------------------------------------------

            uint
            num_intpoints() const ;

//------------------------------------------------------------------------------

            real
            timestep() const ;

            void
            set_model_parameters( const Vector< real > & aParams );

//------------------------------------------------------------------------------

            const Vector< real > &
            model_parameters() const ;

//------------------------------------------------------------------------------

            calculator::VectorData *
            create_vector( const string & aLabel,
                           const uint aSize,
                           const EntityType aType = EntityType::UNDEFINED );

//------------------------------------------------------------------------------

            calculator::MatrixData *
            create_matrix( const string & aLabel,
                           const uint aNumRows,
                           const uint aNumCols );

//------------------------------------------------------------------------------

            void
            print_dofs();

//------------------------------------------------------------------------------

            void
            print_local_dofs();

//------------------------------------------------------------------------------

            const Vector< real > &
            nedelec_data_h();

//------------------------------------------------------------------------------

            const Vector< real > &
            nedelec_data_a();

//------------------------------------------------------------------------------

            /**
             * links the tape-sideset calculator to the facet of the linked
             * layer element and returns it. Gathers the nodal phi of each
             * volume side that is a phi-region into aPhiM / aPhiS and reports
             * per side whether the volume is an h-conductor instead ( its
             * trace then comes from nedelec_data_master_h / _slave_h )
             */
            Calculator *
            get_normal_calculator(
                Vector< real > & aPhiM,
                Vector< real > & aPhiS,
                bool           & aMasterIsConductor,
                bool           & aSlaveIsConductor ) ;

//------------------------------------------------------------------------------

            /**
             * edge dofs of the master volume of the linked facet, read from
             * the edge_h mesh field ( the live dof storage on every rank,
             * see q() ) into vector( "nedelec_h" ), which is sized for the
             * master type. Linear elements only ( one dof per edge )
             */
            const Vector< real > &
            nedelec_data_master_h();

//------------------------------------------------------------------------------

            /**
             * slave twin of nedelec_data_master_h(), into
             * vector( "nedelec_h_s" )
             */
            const Vector< real > &
            nedelec_data_slave_h();

//------------------------------------------------------------------------------

            /**
             * true if the volume element sits on an h-conductor block
             * ( DomainType::Conductor ), i.e. its field is edge-interpolated
             * and has no nodal potential
             */
            bool
            volume_is_conductor( const mesh::Element * aVolume ) const ;

//------------------------------------------------------------------------------

            Mesh *
            mesh() ;

//------------------------------------------------------------------------------

            bool
            element_is_linear() const ;

//------------------------------------------------------------------------------

            /**
             * register the kernels for the maxwell data helper. Called by the
             * MaxwellFactory (magnetic only) and by
             * Controller::set_thermal_kernel (both kernels) once everything is
             * in place. The helper itself is created at the end of allocate(),
             * when the work vectors exist; calling this afterwards rebuilds it.
             */
            void
            link_maxwell( Kernel * aMaxwellKernel,
                          Kernel * aThermalKernel = nullptr ) ;

//------------------------------------------------------------------------------

            // helper function for material functions
            // computes the angle between magnetic field and current density
            real
            bj_angle( const Vector< real > & b, const Vector< real > & j,
                real & norm_b, real & norm_j ) const ;

            // helper function for material functions
            // computes the unfolded angle [ 0, pi ] between magnetic field and tape normal
            real
            bn_angle( const Vector< real > & b, const Vector< real > & n,
                real & norm_b ) const ;

            calculator::MaxwellData *
            maxwell() ;

//------------------------------------------------------------------------------

            // getter needed by sideset connector
            EdgeFunction *
            edge_function()
            {
                return mEdgeFunction ;
            }

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // picks the mFunLinkElement dispatcher based on the registered
            // kernels; called from allocate() and again from link_maxwell()
            // when the thermal kernel is attached after allocation
            void
            select_link_element_dispatcher();

            void
            link_element_default( Element * aElement );

            void
            link_element_maxwell( Element * aElement );

            void
            link_element_maxwell_thermal( Element * aElement );

            void
            link_element_thermal_maxwell( Element * aElement );

//------------------------------------------------------------------------------

            // Dispatchers for the flattened slave integration index
            // (cumulative over facet, orientation). The corresponding
            // IntegrationData and EdgeFunction are looked up by this index
            // at Calculator::link() time.

            uint
            slave_integration_index_2d( const mesh::Facet * aFacet ) ;

            uint
            slave_integration_index_tet( const mesh::Facet * aFacet ) ;

            uint
            slave_integration_index_hex( const mesh::Facet * aFacet ) ;

            uint
            slave_integration_index_penta( const mesh::Facet * aFacet ) ;

//------------------------------------------------------------------------------

            // node interpolator for scalar fields
            const Matrix< real > &
            Nscalar( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // node interpolator for scalar fields, master
            const Matrix< real > &
            Nscalar_master( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // node interpolator for scalar fields, slave
            const Matrix< real > &
            Nscalar_slave( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // node interpolator for 2D vector fields
            const Matrix< real > &
            N2D( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // node interpolator for 3D vector fields
            const Matrix< real > &
            N3D( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // gradient operator for scalar fields
            const Matrix< real > &
            Bscalar( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // gradient operator for scalar fields master
            const Matrix< real > &
            Bscalar_master( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // gradient operator for scalar fields slave
            const Matrix< real > &
            Bscalar_slave( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // gradient operator for plane stress
            const Matrix< real > &
            Bplanestress( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // gradient operator for 3d mech
            const Matrix< real > &
            Bvoigt( const uint aIndex ) ;

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_tri_straight( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_tri_curved( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_quad_straight( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_quad_curved( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_tet_straight( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_tet_curved( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_penta( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_hex( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_tri6_tet10( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_ts( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_quad4ts( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_hex( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_axsymmx( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_axsymmy( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_cartesian( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_axsymmx( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_axsymmy( const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            nedelec_data_linear_h();

            const Vector< real > &
            nedelec_data_quadratic_h_2d();

            const Vector< real > &
            nedelec_data_quadratic_h_3d();

            const Vector< real > &
            nedelec_data_linear_a();

            const Vector< real > &
            nedelec_data_quadratic_a_3d();

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            const Vector< real > &
            nedelec_data_linear( const string & aEdgeField );

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            const Vector< real > &
            nedelec_data_quadratic_2d( const string & aEdgeField,
                                       const string & aFaceField,
                                       const string & aVectorLabel );

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            const Vector< real > &
            nedelec_data_quadratic_3d( const string & aEdgeField,
                                       const string & aFaceField,
                                       const string & aVectorLabel );

//------------------------------------------------------------------------------

            const Matrix< real > &
            invJ2D3D( const uint aIndex ) ;

//------------------------------------------------------------------------------


            const Matrix< real > &
            invJaxsym( const uint aIndex ) ;

            real
            bj_angle_2d( const Vector< real > & b, const Vector< real > & j, real & norm_b, real & norm_j ) const ;

            real
            bj_angle_3d( const Vector< real > & b, const Vector< real > & j, real & norm_b, real & norm_j ) const ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const string &
        calculator::VectorData::label() const
        {
            return mLabel ;
        }

        inline void
        calculator::VectorData::set_index( const belfem::uint aIndex )
        {
            mIndex = aIndex ;
        }

        inline uint
        calculator::VectorData::index() const
        {
            return mIndex ;
        }

        inline Vector< real > &
        calculator::VectorData::vector()
        {
            return mVectorData ;
        }

        inline EntityType
        calculator::VectorData::entity_type() const
        {
            return mType ;
        }

        inline const string &
        calculator::MatrixData::label() const
        {
            return mLabel ;
        }

        inline void
        calculator::MatrixData::set_index( const uint aIndex )
        {
            mIndex = aIndex ;
        }

        inline uint
        calculator::MatrixData::index() const
        {
            return mIndex ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        calculator::MatrixData::matrix()
        {
            return mMatrixData ;
        }

//------------------------------------------------------------------------------

        inline bool Calculator::vector_exists( const string aLabel ) const
        {
            return  mVectorMap.key_exists( aLabel );
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Calculator::vector( const string & aLabel )
        {
            return mVectorMap( aLabel )->vector() ;
        }

//------------------------------------------------------------------------------

        inline bool Calculator::matrix_exists( const string aLabel ) const
        {
            return  mMatrixMap.key_exists( aLabel );
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Calculator::matrix( const string & aLabel )
        {
            return mMatrixMap( aLabel )->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::X() const
        {
            return mX ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Xm() const
        {
            return mXm ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Xs() const
        {
            return mXs ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::J( const uint aIndex )
        {
            if ( aIndex != mJ->index() )
            {
                mJ->set_index( aIndex );
                mJ->matrix().matrix_data() = mIsCurved ?
                        mDomainIntegration->dNdXi( aIndex ) * mX :
                        mLinearIntegration->dNdXi( aIndex ) * mXc ;
            }

            return mJ->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Jm( const uint aIndex )
        {
            if ( aIndex != mJm->index() )
            {
                mJm->set_index( aIndex );
                mJm->matrix().matrix_data() =
                        mMasterIntegration->dNdXi( aIndex ) * mXm ;
            }

            return mJm->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Js( const uint aIndex )
        {
            if ( aIndex != mJs->index() )
            {
                mJs->set_index( aIndex );
                mJs->matrix().matrix_data() =
                        mSlaveIntegration->dNdXi( aIndex ) * mXs ;
            }

            return mJs->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::invJ( const uint aIndex )
        {
            return ( this->*mFunInvJ )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::invJ2D3D( const uint aIndex )
        {
            if ( aIndex != mInvJ->index() )
            {
                // remember index
                mInvJ->set_index( aIndex );

                // in 2D and 3D, we can directly use this value for dV
                mDetJIndex = aIndex ;

                // compute inverse and remember determinant
                mDetJ = ( * mFunInvertJ )( this->J( aIndex ), mInvJ->matrix() );
            }

            return mInvJ->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::invJaxsym( const uint aIndex )
        {
            if ( aIndex != mInvJ->index() )
            {
                // remember index
                mInvJ->set_index( aIndex );

                // just invert the matrix but do not store the determinant
                ( * mFunInvertJ )( this->J( aIndex ), mInvJ->matrix() );
            }

            return mInvJ->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Nscalar( const uint aIndex )
        {
            return mDomainIntegration->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Nscalar_master( const uint aIndex )
        {
            return mMasterIntegration->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Nscalar_slave( const uint aIndex )
        {
            return mSlaveIntegration->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::N2D( const uint aIndex )
        {
            if( mN->index() != aIndex )
            {
                // precomputed data
                const Vector< real > & tPhi = mDomainIntegration->phi( aIndex );

                // remember the index
                mN->set_index( aIndex );

                // get link to matrix
                Matrix< real > & tN = mN->matrix() ;

                // initialize counter
                uint tCount = 0 ;

                // populate matrix
                for( uint k=0; k<mNumberOfNodes; ++k )
                {
                    tN( 0, tCount++ ) = tPhi( k );
                    tN( 1, tCount++ ) = tPhi( k );
                }
            }
            return mN->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::N3D( const uint aIndex )
        {
            if( mN->index() != aIndex )
            {
                // precomputed data
                const Vector< real > & tPhi = mDomainIntegration->phi( aIndex );

                // remember the index
                mN->set_index( aIndex );

                // get link to matrix
                Matrix< real > & tN = mN->matrix() ;

                // initialize counter
                uint tCount = 0 ;

                // populate matrix
                for( uint k=0; k<mNumberOfNodes; ++k )
                {
                    tN( 0, tCount++ ) = tPhi( k );
                    tN( 1, tCount++ ) = tPhi( k );
                    tN( 2, tCount++ ) = tPhi( k );
                }
            }
            return mN->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bscalar( const uint aIndex )
        {
            if ( aIndex != mB->index() )
            {
                // remember the index
                mB->set_index( aIndex );

                mB->matrix() =
                        this->invJ( aIndex ) * mDomainIntegration->dNdXi( aIndex );
            }

            return mB->matrix() ;
        }


//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bscalar_master( const uint aIndex )
        {

            if ( aIndex != mBm->index() )
            {
                // remember the index
                mBm->set_index( aIndex );

                ( * mFunInvertJ )(this->Jm( aIndex ), mInvJm->matrix() ) ;

                mBm->matrix() =mInvJm->matrix() * mMasterIntegration->dNdXi( aIndex );
            }

            return mBm->matrix() ;
        }

        inline const Matrix< real > &
        Calculator::Bscalar_slave( const uint aIndex )
        {
            if ( aIndex != mBs->index() )
            {
                // remember the index
                mBs->set_index( aIndex );

                ( * mFunInvertJ )(this->Js( aIndex ), mInvJs->matrix() ) ;

                mBs->matrix() =  mInvJs->matrix() * mSlaveIntegration->dNdXi( aIndex );
            }

            return mBs->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bplanestress( const uint aIndex )
        {
            if ( aIndex != mB->index() )
            {
                // remember the index
                mB->set_index( aIndex );

                // compute derivative for scalar field
                Matrix< real > & tdN = mdN->matrix() ;

                // compute derivatives
                tdN = this->invJ( aIndex ) * mDomainIntegration->dNdXi( aIndex );

                // get matrix object
                Matrix< real > & tB = mB->matrix() ;

                // initialize counter
                uint tCount = 0 ;

                // populate data
                for( uint k=0; k<mNumberOfNodes; ++k )
                {
                    tB( 0, tCount )   = tdN( 0, k );
                    tB( 2, tCount++ ) = tdN( 1, k );
                    tB( 1, tCount )   = tdN( 1, k );
                    tB( 2, tCount++ ) = tdN( 0, k );
                }
            }

            return mB->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bvoigt( const uint aIndex )
        {
            if ( aIndex != mB->index() )
            {
                // remember the index
                mB->set_index( aIndex );

                // compute derivative for scalar field
                Matrix< real > & tdN = mdN->matrix() ;

                // compute derivatives
                tdN = this->invJ( aIndex ) * mDomainIntegration->dNdXi( aIndex );

                // get matrix object
                Matrix< real > & tB = mB->matrix() ;

                // initialize counter
                uint tCount = 0 ;

                // populate data
                for( uint k=0; k<mNumberOfNodes; ++k )
                {
                    tB( 0, tCount )   = tdN( 0, k );
                    tB( 4, tCount )   = tdN( 2, k );
                    tB( 5, tCount )   = tdN( 1, k );
                    ++tCount ;

                    tB( 1, tCount )   = tdN( 1, k );
                    tB( 3, tCount )   = tdN( 2, k );
                    tB( 5, tCount )   = tdN( 0, k );
                    ++tCount ;

                    tB( 2, tCount )   = tdN( 2, k );
                    tB( 3, tCount )   = tdN( 1, k );
                    tB( 4, tCount )   = tdN( 0, k );
                    ++tCount ;
                }
            }

            return mB->matrix() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::N( const uint aIndex )
        {
            return ( this->*mFunN )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Nm( const uint aIndex )
        {
            return ( this->*mFunNm )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Ns( const uint aIndex )
        {
            return ( this->*mFunNs )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bm( const uint aIndex )
        {
            return ( this->*mFunBm )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::Bs( const uint aIndex )
        {
            return ( this->*mFunBs )( aIndex );
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::node_interp( const uint aIndex, const Vector< real > & aNodeValues ) const
        {
            return dot( mDomainIntegration->phi( aIndex ),  aNodeValues );
        }

//------------------------------------------------------------------------------
        inline const Vector< real > &
        Calculator::Nvec( const uint aIndex ) const
        {
            return mDomainIntegration->phi( aIndex ) ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::B( const uint aIndex )
        {
            return ( this->*mFunB )( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::normal( const uint aIndex )
        {
            return  ( this->*mFunNormal) ( aIndex );
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dS( const uint aIndex )
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            real adS = ( this->*mFundS )( aIndex );
            BELFEM_ASSERT( adS >= 0.0, "Negative Jacobian determinant" );
            return adS;
#else
            return ( this->*mFundS )( aIndex );
#endif
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV( const uint aIndex )
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            real adV = ( this->*mFundV)( aIndex );
            BELFEM_ERROR( adV >= 0.0, "Negative Jacobian determinant" );
            return adV;
#else
            return ( this->*mFundV)( aIndex );
#endif
        }

//------------------------------------------------------------------------------

        // tri and tet only!
        inline real
        Calculator::dV_tri6_tet10( const uint aIndex )
        {
            // we don't need to recompute this at every point for linear elements
            uint tIndex = mIsCurved ? aIndex : 0 ;

            if( mDetJIndex != tIndex )
            {
                mDetJIndex = tIndex ;
                mDetJ = det( this->J( tIndex ) );
            }
            return mDetJ ;
        }

        inline real
        Calculator::dV_ts( const uint aIndex )
        {
            // Formulation discriminator (not element-specific):
            //   - Magnetic solve: an edge function is attached, and its
            //     most recent C(aIndex) call has already cached the
            //     Jacobian determinant at this integration point in
            //     mDetJ. Reuse it.
            //   - Thermal solve (or any Nedelec-free path): no edge
            //     function; fall back to computing det(J) from the
            //     scalar Lagrange shape on demand.
            // Used by PENTA6TS (where the edge function caches
            // thickness*surface as the element-level volume) and by
            // HEX8TS (where the edge function caches the per-IP
            // determinant of the 3D Jacobian); in both cases the
            // edge-function value is the correct dV weight.
            if (mEdgeFunction != nullptr)
            {
                return mEdgeFunction->det_J() ;
            }
            else
            {
                uint tIndex = mIsCurved ? aIndex : 0 ;

                if( mDetJIndex != tIndex )
                {
                    mDetJIndex = tIndex ;
                    mDetJ = det( this->J( tIndex ) );
                }
                return mDetJ ;
            }
        }

        inline real
        Calculator::dV_quad4ts( const uint aIndex )
        {
            // For any physically meaningful ( positive ) layer thickness,
            // QUAD4TS layer elements are wound CLOCKWISE by construction, so
            // their scalar Lagrange determinant is negative. A negative deck
            // thickness is currently accepted and inverts the stack, which
            // flips that sign; this function is correct either way, but the
            // reasoning below assumes the valid case. That winding is
            // forced, not accidental: the layer element's local edge node
            // order must follow the facet edge direction, because the edge
            // signs are read from the node order ( Element::compute_edge_
            // directions ) while the Nedelec tangent is read from the facet
            // ( EF_QUAD4TS::link ). Given the clockwise extrusion normal of
            // ThinShellFactory::process_nodes_line2, ( bottom0, bottom1,
            // top1, top0 ) is the only ordering that satisfies both, and it
            // is left-handed.
            //
            // The map is still a diffeomorphism, so the thermal assembly is
            // correct up to the volume weight: N is reference-space and
            // B = invJ * dNdXi returns true Cartesian gradients whichever
            // orientation the map has. Change of variables then asks for
            // | det J |, which is what this returns -- an exact weight, not
            // a clamp.
            //
            // The MAGNETIC solve reaches this function too, through its own
            // Calculator for the same block -- the dispatch is by element
            // type, so do NOT assume magnetic never enters here. It takes
            // dV_ts's edge-function branch, where det_J is
            // 0.25 * thickness * length: positive for any thickness that
            // reaches assembly, so the abs is the identity and magnetic
            // results are unchanged. A negative deck thickness would flip
            // that, but it cannot get this far -- Kernel::compute_element_
            // volumes forms facet_area * thickness per layer element and
            // hard-errors on a negative total before any assembly runs.
            //
            // Wraps dV_ts rather than duplicating its body ON PURPOSE. Both
            // functions share the mDetJ / mDetJIndex cache with invJ2D3D,
            // which the thermal kernels populate by calling B( k ) BEFORE
            // dV( k ). At k = 0 that leaves a valid cache entry, so a copy
            // of the body that only wrapped its own  mDetJ = det( ... )
            // assignment would return the SIGNED cached value and never run
            // the abs. Taking it on the returned value is immune to that.
            return std::abs( this->dV_ts( aIndex ) );
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV_hex( const uint aIndex )
        {
            if ( mEdgeFunction != nullptr )
            {
                mEdgeFunction->update_nabla( aIndex );
                return mEdgeFunction->det_J() ;
            }
            else
            {

                if( mDetJIndex != aIndex )
                {
                    mDetJIndex = aIndex ;
                    mDetJ = det( this->J( aIndex ) );
                }
                return mDetJ ;
            }
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV_axsymmx( const uint aIndex )
        {
            if( mDetJIndex != aIndex )
            {
                mDetJIndex = aIndex ;

                // radius contribution
                mDetJ = mIsCurved ? dot(
                        mDomainIntegration->phi( aIndex ).vector_data(), mX.col( 1 ) ) :
                        dot( mLinearIntegration->phi( aIndex ).vector_data(), mXc.col( 1 ) );

                mDetJ *= det( this->J( aIndex ) ) * 2.0 * constant::pi ;
            }
            return mDetJ ;
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV_axsymmy( const uint aIndex )
        {

            if( mDetJIndex != aIndex )
            {
                mDetJIndex = aIndex ;

                // radius contribution
                mDetJ = mIsCurved ? dot(
                        mDomainIntegration->phi( aIndex ).vector_data(), mX.col( 0 ) ) :
                        dot( mLinearIntegration->phi( aIndex ).vector_data(), mXc.col( 0 ) );

                mDetJ *= det( this->J( aIndex ) ) * 2.0 * constant::pi ;
            }
            return mDetJ ;
        }

//------------------------------------------------------------------------------


        inline real
        Calculator::dS_cartesian ( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            // note that the surface increment is scaled reciprocally to the
            // integration weights.
            return mSurfaceIncrement ;
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dS_axsymmx( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            return dot( mMasterIntegration->phi( aIndex ).vector_data(),
                        mXm.col( 1 ) )
                    * mSurfaceIncrement  * 2.0 * constant::pi ;

        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dS_axsymmy( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            return dot( mMasterIntegration->phi( aIndex ).vector_data(),
                        mXm.col( 0 ) )
                   * mSurfaceIncrement * 2.0 * constant::pi ;

        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::integration() const
        {
            return mDomainIntegration ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::master_integration() const
        {
            return mMasterIntegration ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::slave_integration() const
        {
            return mSlaveIntegration ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::volume_integration() const
        {
            return mMasterVolumeIntegration ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::volume_enrichment() const
        {
            return mVolumeEnrichment ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Calculator::sideset_enrichment() const
        {
            return mSideSetEnrichment ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Calculator::K()
        {
            return mK ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Calculator::M()
        {
            return mM ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Calculator::JN()
        {
            return mJN ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Calculator::f()
        {
            return mf ;
        }

//------------------------------------------------------------------------------

        inline uint
        Calculator::num_intpoints() const
        {
            return mNumberOfIntegrationPoints ;
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::timestep() const
        {
            return mTimestep ;
        }

//------------------------------------------------------------------------------
        inline void
        Calculator::set_model_parameters( const Vector< real > & aParams )
        {
            mModelParameters = aParams ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::model_parameters() const
        {
            return mModelParameters ;
        }

//------------------------------------------------------------------------------

        inline Element *
        Calculator::element()
        {
            return mElement ;
        }

//------------------------------------------------------------------------------

        inline Group *
        Calculator::group()
        {
            return mGroup ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_h()
        {
            return (this->*mFunNedelecDataH )();
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_a()
        {
            return (this->*mFunNedelecDataA )();
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_linear_h()
        {
            return this->nedelec_data_linear( "edge_h");
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_linear_a()
        {
            return this->nedelec_data_linear( "edge_a");
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_quadratic_h_2d()
        {
            return this->nedelec_data_quadratic_2d(
                "edge_h",
                "face_h",
                "nedelec_h");
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_quadratic_h_3d()
        {
            return this->nedelec_data_quadratic_3d(
            "edge_h",
            "face_h",
            "nedelec_h");
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Calculator::nedelec_data_quadratic_a_3d()
        {
            return this->nedelec_data_quadratic_3d(
            "edge_a",
            "face_a",
            "nedelec_a");
        }

//------------------------------------------------------------------------------

        inline uint
        Calculator::integration_order() const
        {
            return mIntegrationOrder ;
        }

//------------------------------------------------------------------------------

        inline Mesh * Calculator::mesh()
        {
            return mMesh ;
        }

//------------------------------------------------------------------------------

        inline bool Calculator::element_is_linear() const
        {
            return mIsLinear ;
        }

        inline real
        Calculator::bj_angle( const Vector< real > & b, const Vector< real > & j, real & norm_b, real & norm_j ) const
        {
            return ( this->*mFunBJAngle )( b, j, norm_b, norm_j );
        }

        inline real
        Calculator::bn_angle( const Vector< real > & b, const Vector< real > & n, real & norm_b ) const
        {
            norm_b = norm( b ) ;
            BELFEM_ASSERT( std::abs( norm( n ) - 1.0 ) < 1e-6, "tape normal must be normalized" ) ;

            // UNFOLDED angle between field and tape normal, [ 0, pi ]
            // ( 2026-08-16 ). theta < pi/2 iff b has a component
            // along +n; n is the outward normal of the mid-surface facet's
            // MASTER volume element ( master -> slave ), which is also the
            // layer-stack direction — a minus prefix on the thinshell
            // sidesets key flips both together. Measured jc(theta) tables
            // are asymmetric about pi/2 and consume this angle as is;
            // analytic laws that are even in theta ( ModifiedKim ) fold
            // internally by construction
            return norm_b < 1e-6 ? constant::pi*0.5 :
                std::acos( std::clamp( dot( n, b ) / norm_b, -1.0, 1.0 ) ) ;
        }

        inline real
        Calculator::bj_angle_3d( const Vector< real > & b, const Vector< real > & j, real & norm_b, real & norm_j ) const
        {
            norm_b = norm( b ) ;
            norm_j = norm( j ) ;

            // Angle between magnetic field and current density (beta)
            return norm_b < 1e-6 || norm_j < 1e-6 ?
                        constant::pi*0.5 : std::acos( std::min(std::abs( dot( b, j ) / ( norm_b * norm_j ) ) , 1.0 )) ;
        }

        inline real
        Calculator::bj_angle_2d( const Vector< real > & b, const Vector< real > & j, real & norm_b, real & norm_j ) const
        {
            norm_b = norm( b ) ;
            norm_j = norm( j ) ;

            // Angle between magnetic field and current density (beta)
            // in 2D: current is always perpendicular to model plane
            return constant::pi*0.5 ;
        }

        inline Vector< real > &
        calculator::MaxwellData::link_vector( const string & aLabel )
        {
            if ( ! mMaxwellCalculator->vector_exists( aLabel ) )
            {
                uint tNumDim = mMaxwellCalculator->mesh()->number_of_dimensions() ;

                // for 2D calculations, j is always perpendicular to model plane
                // this is the only special case that we have to catch here
                if ( aLabel == "j" && tNumDim == 2 )
                {
                    tNumDim = 1 ;
                }

                VectorData * tData = mMaxwellCalculator->create_vector( aLabel, tNumDim );
                return tData->vector();
            }
            return mMaxwellCalculator->vector( aLabel );
        }

        inline void
        calculator::MaxwellData::compute_x( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::x, aIndex ) )
            {
                ( this->*mFunX )( aIndex );
                this->set( MaxwellDataValue::x, aIndex );
            }
        }

        inline void
        calculator::MaxwellData::compute_x_2d( const uint aIndex )
        {
            mCoords = mMaxwellCalculator->N( aIndex ) * mMaxwellCalculator->X();
            mX = mCoords( 0, 0 );
            mY = mCoords( 0, 1 );
        }

        inline void
        calculator::MaxwellData::compute_x_3d( const uint aIndex )
        {
            mCoords = mMaxwellCalculator->N( aIndex ) * mMaxwellCalculator->X();
            mX = mCoords( 0, 0 );
            mY = mCoords( 0, 1 );
            mZ = mCoords( 0, 2 );
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_j( const uint aIndex )
        {
            if ( aIndex == mLastIndex( static_cast< uint > ( MaxwellDataValue::j ) ) )
            {
                return mJ ;
            }

            mLastIndex( static_cast< uint > ( MaxwellDataValue::j ) ) = aIndex ;
            mJ = mMaxwellCalculator->C( aIndex ) * mMaxwellCalculator->q();

            return mJ ;
        }


        /**
         * Full field trace of one volume side of a thin-shell facet, for
         * linear elements ( the k = 0 contract of compute_hn ):
         *
         *  - phi-side ( node-interpolated volume: air, ferro, buffer ):
         *    h = -B_sigma( 0 ) * phi_sigma, constant on a linear element
         *    ( Alves et al. 2022b, h0 = -grad phi-, hN = -grad phi+ );
         *
         *  - h-side ( DomainType::Conductor, edge-interpolated, no phi ):
         *    h = integration-weighted mean over the facet rule of
         *    E_sigma( k ) * q_sigma, the conductor's own Nedelec trace.
         *    A lowest-order Nedelec field is a + b x r, so its normal
         *    component is linear over the facet and the weighted mean is
         *    the exact centroid value in 2D and 3D, independent of the
         *    rule's point ordering ( in 3D point 0 of the default 7-point
         *    rule already is the centroid; in 2D it is not ). The rule
         *    rests on [ n . B ] = 0 across the sheet with mu = mu0 on the
         *    conductor side, so n . h is continuous — the factory refuses
         *    a magnetic conductor next to a shell for that reason
         *    ( MaxwellFactory::assign_materials ).
         *
         * aNormalCalc : the tape-sideset calculator from get_normal_calculator()
         * aMaster     : master ( true ) or slave ( false ) side
         * aIsConductor: true if that side is an h-conductor ( edge trace ), false for a phi-region
         * aPhi        : nodal phi of that side ( phi-side only )
         * aScratch    : a d-vector on the caller's calculator
         * aH          : result, a d-vector on the caller's calculator
         */
        inline void
        compute_h_trace(
                Calculator           * aNormalCalc,
                const bool             aMaster,
                const bool             aIsConductor,
                const Vector< real > & aPhi,
                Vector< real >       & aScratch,
                Vector< real >       & aH )
        {
            if ( aIsConductor )
            {
                const Vector< real > & q = aMaster ?
                        aNormalCalc->nedelec_data_master_h() :
                        aNormalCalc->nedelec_data_slave_h() ;

                const Vector< real > & w = aNormalCalc->integration()->weights() ;

                const uint tNumPoints = aNormalCalc->num_intpoints() ;

                real tSum = 0.0 ;
                aH.fill( 0.0 );

                for ( uint k = 0; k < tNumPoints; ++k )
                {
                    aScratch = ( aMaster ? aNormalCalc->Em( k ) : aNormalCalc->Es( k ) ) * q ;
                    aH += w( k ) * aScratch ;
                    tSum += w( k );
                }
                aH /= tSum ;
            }
            else
            {
                aH = -1. * ( aMaster ? aNormalCalc->Bm( 0 ) : aNormalCalc->Bs( 0 ) ) * aPhi ;
            }
        }

        /**
         * Computes the purely-normal magnetic field hn at the thin-shell
         * facet: the average of the master and slave volume traces
         * ( compute_h_trace, per-side dispatch phi-region / h-conductor ),
         * projected onto the facet normal. Result is written into
         * aCalc->vector("hn") and also returned. For flat linear shells
         * call once with k=0 outside the loop.
         */
        inline const Vector< real > &
        compute_hn( Calculator * aCalc , const uint k )
        {
            // normal field
            Vector< real > & hn = aCalc->vector("hn");

            // warning:
            // for linear elements, we only need to do this once
            BELFEM_ASSERT( aCalc->element_is_linear(), "compute_hn doesn't work for non-linear elements because it expects surface integration points but we integrate over the volume" );

            if ( k == 0 )
            {
                // dofs from master element
                Vector< real > & phi_m = aCalc->vector("phi_m");

                // dofs from slave element
                Vector< real > & phi_s = aCalc->vector("phi_s");

                // get calculator for facet and the kind of each volume side
                bool tMasterIsConductor ;
                bool tSlaveIsConductor ;

                Calculator * tCalc = aCalc->get_normal_calculator(
                        phi_m, phi_s, tMasterIsConductor, tSlaveIsConductor );

                // scratch for the per-point products
                Vector< real > & hk = aCalc->vector("hk");

                // h-field on master side
                Vector< real > & hm = aCalc->vector("hm");
                compute_h_trace( tCalc, true, tMasterIsConductor, phi_m, hk, hm );

                // h-field on slave side
                Vector< real > & hs = aCalc->vector("hs");
                compute_h_trace( tCalc, false, tSlaveIsConductor, phi_s, hk, hs );

                // temporarily writing average into hn container
                hn = 0.5 * ( hm + hs ) ;

                // normal vector
                Vector< real > & n = aCalc->vector("normal");
                n = tCalc->normal( k );

                // normal component of h
                hn = dot( hn, n ) * n ;
            }
            return hn;
        }

        inline void
        Calculator::link_element_default( Element * aElement )
        {
            mElement = aElement;
        }



        inline void
        calculator::MaxwellData::reset()
        {
            mLastIndex.fill( BELFEM_UINT_MAX );
            mBeta   = 0.5 * constant::pi ;
            mTClamped   = false ;
            mRhoClamped = false ;
            mFrameCurrent = false ;
        }

        inline bool
        calculator::MaxwellData::is_current( MaxwellDataValue aValue, const uint aIndex ) const
        {
            return mLastIndex( static_cast< uint >( aValue ) ) == aIndex ;
        }

        inline void
        calculator::MaxwellData::set( MaxwellDataValue aValue, const uint aIndex )
        {
            mLastIndex( static_cast< uint >( aValue ) ) = aIndex ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_h( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::H, aIndex ) )
            {
                this->set( MaxwellDataValue::H, aIndex );
                return ( this->*mFunH )( aIndex );

            }
            return mH ;
        }

        inline real
        calculator::MaxwellData::compute_mu( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::mu, aIndex ) )
            {
                this->set( MaxwellDataValue::mu, aIndex );
                mMu = ( this->*mFunMu )( aIndex );
            }
            return mMu ;
        }

        inline real
        calculator::MaxwellData::compute_mu_0( const uint aIndex )
        {
            return constant::mu0 ;
        }

        inline real
        calculator::MaxwellData::compute_mu_const( const uint aIndex )
        {
            return mMaterial->constant_property( MaterialProperty::mu );
        }

        inline real
        calculator::MaxwellData::compute_mu_h( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::normH, aIndex ) )
            {
                mNormH = norm( this->compute_h( aIndex ));
                this->set( MaxwellDataValue::normH, aIndex );
            }
            return mMaterial->mu( mNormH );
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_b( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::B, aIndex ) )
            {
                this->set( MaxwellDataValue::B, aIndex );
                return ( this->*mFunB )( aIndex );
            }
            return mB ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_b_bulk( const uint aIndex )
        {
            mB = this->compute_mu( aIndex ) * this->compute_h( aIndex ) ;
            return mB ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_b_ts( const uint aIndex )
        {
            this->compute_h( aIndex );
            real mu = this->compute_mu( aIndex );
            mBt = mu * mHt ;
            mBn = mu * mHn ;
            mB = mBt + mBn ;
            return mB ;
        }

        inline real
        calculator::MaxwellData::compute_rho( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::rho, aIndex ) )
            {
                real tRho = ( this->*mFunRho )( aIndex );

                // numerical guard on the power-law output: one clamped value
                // for all consumers ( K, Joule source, element mean ). While
                // clamped, the J-derivative is zero ( consistent tangent )
                mRhoClamped = ( tRho < gRhoMin ) || ( tRho > gRhoMax ) ;
                mRho = mRhoClamped ? std::clamp( tRho, gRhoMin, gRhoMax ) : tRho ;

                this->set( MaxwellDataValue::rho, aIndex );
            }
            return mRho ;
        }


        inline real
        calculator::MaxwellData::compute_drhodj( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::drhodj, aIndex ) )
            {
                // make sure rho — and with it mRhoClamped — is current here
                this->compute_rho( aIndex );

                // consistent tangent: a clamped rho has zero derivative
                mdRhodJ = mRhoClamped ? 0.0 : ( this->*mFundRhodJ )( aIndex );
                this->set( MaxwellDataValue::drhodj, aIndex );
            }
            return mdRhodJ ;
        }

        inline real
        calculator::MaxwellData::compute_drhodT( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::drhodT, aIndex ) )
            {
                // make sure rho — and with it mRhoClamped — is current here
                this->compute_rho( aIndex );

                // refresh mTClamped for THIS index ( cache-hit if the rho
                // law above already computed it ): a cached rho would
                // otherwise leave the flag at whatever index last touched
                // compute_T — same self-sufficiency contract as
                // compute_dcpdT ( Grok phase-3 hardening, 2026-08-13 )
                this->compute_T( aIndex );

                // consistent tangent: a clamped rho has zero derivative, and
                // a clamped T iterate zeroes ALL dT derivatives — same
                // contract as compute_dcpdT / compute_dlambdadT ( Newton
                // must not push on a flat clamp; 2026-08-13 audit )
                mdRhodT = ( mTClamped || mRhoClamped ) ? 0.0 : ( this->*mFundRhodT )( aIndex );
                this->set( MaxwellDataValue::drhodT, aIndex );
            }
            return mdRhodT ;
        }


        inline real
        calculator::MaxwellData::compute_T( const uint aIndex )
        {
            return ( this->*mFunT )( aIndex );
        }

        inline real
        calculator::MaxwellData::compute_cp( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::cp, aIndex ) )
            {
                mCp = mMaterial->cp( this->compute_T( aIndex ) );
                this->set( MaxwellDataValue::cp, aIndex );
            }
            return mCp ;
        }

        inline real
        calculator::MaxwellData::compute_dcpdT( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::dcpdT, aIndex ) )
            {
                real T = this->compute_T( aIndex );

                // consistent tangent: no T-sensitivity while T is clamped
                mdCpdT = mTClamped ? 0.0 : mMaterial->dcpdT( T );
                this->set( MaxwellDataValue::dcpdT, aIndex );
            }
            return mdCpdT ;
        }

        inline real
        calculator::MaxwellData::compute_lambda( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::lambda, aIndex ) )
            {
                mLambda = ( this->*mFunLambda)( aIndex );
                this->set( MaxwellDataValue::lambda, aIndex );
            }
            return mLambda ;
        }

        inline real
        calculator::MaxwellData::compute_dlambdadT( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::dlambdadT, aIndex ) )
            {
                // the variant computes T internally and updates mTClamped
                mdLambdadT = ( this->*mFundLambdadT)( aIndex );

                // consistent tangent: no T-sensitivity while T is clamped
                if ( mTClamped )
                {
                    mdLambdadT = 0.0 ;
                }
                this->set( MaxwellDataValue::dlambdadT, aIndex );
            }
            return mdLambdadT ;
        }

        inline real
        calculator::MaxwellData::compute_T_fem( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::T, aIndex ) )
            {
                real T = dot( mThermalCalculator->Nvec( aIndex ) , mThermalCalculator->q() ) ;

                // transient nonlinear iterates can swing outside the physical
                // table window; clamp on both edges and remember: while
                // clamped, all dT-derivatives are zero ( consistent tangent ).
                // A CONVERGED solution at a clamp is a modeling error — the
                // controller can query T_clamped() for a diagnostic.
                mTClamped = ( T < gTmin ) || ( T > mTmax ) ;
                mT = mTClamped ? std::clamp( T, gTmin, mTmax ) : T ;

                this->set( MaxwellDataValue::T, aIndex );
            }
            return mT ;
        }

        inline real
        calculator::MaxwellData::compute_T_const( const uint aIndex )
        {
            // same contract as compute_T_fem: a user-set bulk temperature
            // outside the material's table window is pinned, and the flag
            // zeroes the dT-derivatives ( consistent tangent )
            mTClamped = ( gTbulk < gTmin ) || ( gTbulk > mTmax );
            mT = mTClamped ? std::clamp( gTbulk, gTmin, mTmax ) : gTbulk ;
            return mT ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_h_bulk_node( const uint aIndex )
        {
            mH = -1. * mMaxwellCalculator->B( aIndex ) * mMaxwellCalculator->q();
            return mH ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_h_bulk_edge( const uint aIndex )
        {
            mH = mMaxwellCalculator->E( aIndex ) * mMaxwellCalculator->q();
            return mH ;
        }

        inline const Vector< real > &
        calculator::MaxwellData::compute_h_ts_edge( const uint aIndex )
        {
            mHn = compute_hn( mMaxwellCalculator, aIndex );
            this->set( MaxwellDataValue::n, aIndex ); // compute_hn writes n
            mHt = mMaxwellCalculator->E( aIndex ) * mMaxwellCalculator->q();

            mH = mHn+mHt ;
            return mH ;
        }

        inline real
        calculator::MaxwellData::compute_lambda_bulk( const uint aIndex )
        {
            return mMaterial->lambda( this->compute_T( aIndex ) );
        }

        inline real
        calculator::MaxwellData::compute_dlambdadT_bulk( const uint aIndex )
        {
            return mMaterial->dlambdadT( this->compute_T( aIndex ) );
        }

        inline real
        calculator::MaxwellData::compute_lambda_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->lambda( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }

        inline real
        calculator::MaxwellData::compute_dlambdadT_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->dlambdadT( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }



        inline real
        calculator::MaxwellData::compute_rho_bulk( const uint aIndex )
        {
            return mMaterial->rho( this->compute_T( aIndex ) );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_bulk( const uint aIndex )
        {
            return mMaterial->drhodT( this->compute_T( aIndex ) );
        }

        inline real
        calculator::MaxwellData::norm_b( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::normB, aIndex ) )
            {
                mNormB = norm( this->compute_b( aIndex ) );
                this->set( MaxwellDataValue::normB, aIndex );
            }
            return mNormB ;
        }

        inline real
        calculator::MaxwellData::norm_j( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::normJ, aIndex ) )
            {
                mNormJ = norm( this->compute_j( aIndex ) );
                this->set( MaxwellDataValue::normJ, aIndex );
            }
            return mNormJ ;
        }

        inline real
        calculator::MaxwellData::compute_rho_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->rho( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }

        // HTS drho/dT family ( T-leg ): preambles mirror the
        // corresponding compute_drhodb_* wrappers — bn_angle block for the
        // thin-shell variants, beta_dummy for bulk, compute_x for defect

        inline real
        calculator::MaxwellData::compute_drhodT_powerlaw_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_powerlaw_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_powerlaw_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_piecewise_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_piecewise_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_piecewise_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_piecewise_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_powerlaw_bulk( const uint aIndex )
        {
            return mMaterial->drho_powerlaw_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_powerlaw_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_piecewise_bulk( const uint aIndex )
        {
            return mMaterial->drho_piecewise_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_piecewise_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_piecewise_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

//------------------------------------------------------------------------------
// riva-law wrappers ( 2026-08-27 ), mirroring the piecewise family
//------------------------------------------------------------------------------

        inline real
        calculator::MaxwellData::compute_rho_riva_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->rho_riva(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_rho_riva_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->rho_riva(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_rho_riva_bulk( const uint aIndex )
        {
            return mMaterial->rho_riva(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_rho_riva_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->rho_riva(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_riva_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_riva_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_riva_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_riva_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_riva_bulk( const uint aIndex )
        {
            return mMaterial->drho_riva_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_riva_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_riva_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_riva_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_riva_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_riva_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_riva_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_riva_bulk( const uint aIndex )
        {
            return mMaterial->drho_riva_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_riva_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_riva_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_riva_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_riva_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_riva_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_riva_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_riva_bulk( const uint aIndex )
        {
            return mMaterial->drho_riva_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_riva_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_riva_dT(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodT_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->drhodT( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodb( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::drhodb, aIndex ) )
            {
                // make sure rho — and with it mRhoClamped — is current here
                this->compute_rho( aIndex );

                // consistent tangent: a clamped rho has zero derivative
                mdRhodB = mRhoClamped ? 0.0 : ( this->*mFundRhodB )( aIndex );
                this->set( MaxwellDataValue::drhodb, aIndex );
            }
            return mdRhodB ;
        }

        inline real
        calculator::MaxwellData::compute_drhodbeta( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::drhodbeta, aIndex ) )
            {
                // make sure rho — and with it mRhoClamped — is current here
                this->compute_rho( aIndex );

                // consistent tangent: a clamped rho has zero derivative
                mdRhodBeta = mRhoClamped ? 0.0 : ( this->*mFundRhodBeta )( aIndex );
                this->set( MaxwellDataValue::drhodbeta, aIndex );
            }
            return mdRhodBeta ;
        }

        inline real
        calculator::MaxwellData::compute_drhodb_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->drhodB( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodbeta_metal( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );
                const Vector< real > & j = this->compute_j( aIndex );

                // the bj_angle function updates NormB and normJ
                mBeta = mMaxwellCalculator->bj_angle( b, j, mNormB, mNormJ );

                this->set( MaxwellDataValue::normB, aIndex );
                this->set( MaxwellDataValue::normJ, aIndex );
                this->set( MaxwellDataValue::beta, aIndex );
            }

            return mMaterial->drhodbeta( this->compute_T( aIndex ), this->norm_b( aIndex ), mBeta );
        }

//------------------------------------------------------------------------------
// HTS drho/d|B| wrappers: each mirrors the preamble of its
// compute_rho_* counterpart exactly, so the cached beta / normB state is
// updated identically whichever of the pair runs first in an iterate.
//------------------------------------------------------------------------------

        inline real
        calculator::MaxwellData::compute_drhodb_powerlaw_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_powerlaw_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_powerlaw_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_piecewise_ts( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_piecewise_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_piecewise_ts_defect( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_piecewise_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_powerlaw_bulk( const uint aIndex )
        {
            return mMaterial->drho_powerlaw_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_powerlaw_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_piecewise_bulk( const uint aIndex )
        {
            return mMaterial->drho_piecewise_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodb_piecewise_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_piecewise_dB(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_rho_powerlaw_bulk( const uint aIndex )
        {
            return mMaterial->rho_powerlaw(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_rho_piecewise_bulk( const uint aIndex )
        {
            return mMaterial->rho_piecewise(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_powerlaw_bulk( const uint aIndex )
        {
            return mMaterial->drho_powerlaw_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_piecewise_bulk( const uint aIndex )
        {
            return mMaterial->drho_piecewise_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy() );
        }

        inline real
        calculator::MaxwellData::compute_rho_powerlaw_ts( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->rho_powerlaw(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_powerlaw_ts( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_powerlaw_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_rho_piecewise_ts( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->rho_piecewise(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_piecewise_ts( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            return mMaterial->drho_piecewise_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta );
        }

        inline real
        calculator::MaxwellData::compute_rho_powerlaw_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->rho_powerlaw(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_powerlaw_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_rho_powerlaw_ts_defect( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->rho_powerlaw(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_powerlaw_ts_defect( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_powerlaw_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_rho_piecewise_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );
            return mMaterial->rho_piecewise(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_drhodj_piecewise_bulk_defect( const uint aIndex )
        {
            this->compute_x( aIndex );
            return mMaterial->drho_piecewise_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                this->beta_dummy(),
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_rho_piecewise_ts_defect( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->rho_piecewise(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime);
        }


        inline real
        calculator::MaxwellData::compute_drhodj_piecewise_ts_defect( const uint aIndex )
        {

            if ( ! this->is_current( MaxwellDataValue::beta, aIndex ) )
            {
                const Vector< real > & b = this->compute_b( aIndex );

                BELFEM_ASSERT( this->is_current( MaxwellDataValue::n, aIndex ), "normal vector has not been updated" );

                mBeta = mMaxwellCalculator->bn_angle( b, mN, mNormB );

                this->set( MaxwellDataValue::beta, aIndex );
                this->set( MaxwellDataValue::normB, aIndex );
            }

            this->compute_x( aIndex );

            return mMaterial->drho_piecewise_dJ(
                this->norm_j( aIndex ),
                this->compute_T( aIndex ),
                this->norm_b( aIndex ),
                mBeta,
                mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::compute_dmudh( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::dmudh, aIndex ) )
            {
                this->set( MaxwellDataValue::dmudh, aIndex );
                return ( this->*mFundMudH ) ( aIndex );
            }
            return mdMudH ;
        }


        inline real
        calculator::MaxwellData::compute_dmu_zero( const uint aIndex )
        {
            mdMudH = 0.0 ;
            return mdMudH ;
        }

        inline real
        calculator::MaxwellData::compute_dmu_material( const uint aIndex )
        {
            if ( ! this->is_current( MaxwellDataValue::normH, aIndex ) )
            {
                mNormH = norm( this->compute_h( aIndex ) );
                this->set( MaxwellDataValue::normH, aIndex );
            }
            mMaterial->dmudH( mNormH, mMu, mdMudH );
            this->set( MaxwellDataValue::mu, aIndex );
            return mdMudH ;
        }

        inline real
        calculator::MaxwellData::compute_volumetric_heatload( const uint aIndex )
        {
            return ( this->*mFunHeat )( aIndex );
        }

        inline real
        calculator::MaxwellData::compute_heatload_user( const uint aIndex )
        {
            this->compute_x( aIndex );

            return mMaterial->volumetric_heatload( mX, mY, mZ, mTime );
        }

        inline real
        calculator::MaxwellData::return_zero( const uint aIndex )
        {
            return 0.0 ;
        }

        inline real
        calculator::MaxwellData::beta_dummy() const
        {
            // bulk HTS has no meaningful tape normal, so the material call
            // takes a dummy beta that must still be the reset() value.
            // beta conventions: ( b, n ) bn_angle for HTS thin shells,
            // ( b, j ) bj_angle for metals — never both in one material
            // ( enforced in the constructor )
            BELFEM_ASSERT( mBeta == 0.5 * constant::pi,
                "bulk-HTS dummy beta was polluted ( expected pi/2 )" );
            return mBeta ;
        }

        inline real
        calculator::MaxwellData::density() const
        {
            BELFEM_ASSERT( ! std::isnan( mDensity ),
                "material %s provides neither density nor ref_density",
                mMaterial->label().c_str() );
            return mDensity ;
        }

        inline bool
        calculator::MaxwellData::T_clamped() const
        {
            return mTClamped ;
        }

        inline bool
        calculator::MaxwellData::rho_clamped() const
        {
            return mRhoClamped ;
        }

        inline Calculator *
        calculator::MaxwellData::maxwell()
        {
            return mMaxwellCalculator ;
        }

        inline Calculator *
        calculator::MaxwellData::thermal()
        {
            return mThermalCalculator ;
        }

        inline calculator::MaxwellData *
        Calculator::maxwell()
        {
            BELFEM_ASSERT( mMaxwellData != nullptr, "MaxwellData has not been initialized" );
            return mMaxwellData ;
        }


    }
}

#endif //BELFEM_CL_FEM_CALCULATOR_HPP
