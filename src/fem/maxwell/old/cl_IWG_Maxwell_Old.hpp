//
// Created by christian on 10/4/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_OLD_HPP
#define BELFEM_CL_IWG_MAXWELL_OLD_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "en_SolverEnums.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_FEM_DomainGroup.hpp"
#include "cl_FEM_BoundaryCondition.hpp"
#include "en_FEM_MagfieldBcType.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * an abstract baseclass for
         * solving the Maxwell equations
         */
        class IWG_Maxwell_Old : public IWG_Timestep
        {

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            const uint mNumberOfDimensions ;

            // flag telling if this IWG uses edge dofs
            const bool mUseEdges ;

            Map< id_t, DomainType > mBlockTypes ;
            Map< id_t, DomainType > mSideSetTypes ;

            // dofs per superconducting block
            Cell< string > mScDofs ;

            // dofs per coil
            Cell< string > mCoilDofs ;

            // dofs per iron block
            Cell< string > mFerroDofs ;

            // dofs per air
            Cell< string > mAirDofs ;

            // dofs along Cut-Interface
            Cell< string > mCutDofs ;

            // dofs along Air-Solid-Interface
            Cell< string > mInterfaceDofsScAir ;
            Cell< string > mInterfaceDofsScFm ;
            Cell< string > mInterfaceDofsFmAir ;
            Cell< string > mInterfaceDofsCut ;


            // special dofs along domain boundary



            Cell< string > mInterfaceDofs ;
            Cell< string > mFarfieldDofs ;
            Cell< string > mSymmetryDofs ;
            Cell< string > mBoundaryDofs ;

            Cell< string > mLambdaDofs ;

            //! surface increment, computed by normal function
            //! or geometry Jacobian function
            real mDomainIncrement = BELFEM_QUIET_NAN ;

            //! list with blocks that contain currents, needed for postprocessing
            Vector< id_t > mCurrentBlocks ;

            // help values
            Cell< string > mAfields = { "ax", "ay", "az" };
            Cell< string > mAfields0 = { "a0x", "a0y", "a0z" };
            Cell< string > mBfields ;

            //! map telling which subtype a domain boundary is
            Map< id_t, MagfieldBcType > mMagfieldTypeMap ;

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            //! work matrix for facet normal
            Vector< real > mWorkN ;

            //! work vector for facet direction
            Vector< real > mWorkS ;
            Vector< real > mWorkT ;

            //! work matrix for facet direction
            Matrix< real > mWorkR;

            // internal work vectors, updated by nabla function
            Vector< real > mNablaXi ;
            Vector< real > mNablaEta ;
            Vector< real > mNablaZeta ;
            Vector< real > mNablaTau ;

            void
            ( IWG_Maxwell_Old::*mFunNabla )(
                    Element * aElement,
                    const uint aIndex  );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunE )(
                    Element * aElement,
                    const uint aIndex  );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunE2D )(
                    Element * aElement,
                    const real aXi, const real aEta );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunE3D )(
                    Element * aElement,
                    const real aXi, const real aEta, const real aZeta );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunCurlH )(
                    Element * aElement,
                    const uint aIndex  );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunCurlA )(
                    Element * aElement,
                    const uint aIndex  );

            const Matrix< real > &
            ( IWG_Maxwell_Old::*mFunGradPhi )(
                    Element * aElement,
                    const uint aIndex  );

            const Vector< real > &
            ( IWG_Maxwell_Old::*mFunNormal )(
                    Element * aElement,
                    const uint aIndex  );

            void
            ( IWG_Maxwell_Old::*mFunScMatrices )(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            ( IWG_Maxwell_Old::*mFunProjectIntpoints )(
                    Element * aElement,
                          Matrix< real > & aIntpointsOnMaster,
                          Matrix< real > & aIntPointsOnSlave );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_Old(
                    const uint aNumberOfDimensions,
                    const IwgType aType,
                    const IwgMode aMode,
                    const SymmetryMode aSymmetryMode,
                    const SideSetDofLinkMode SideSetDofLinkMode,
                    const bool         aUseEdges = false ) ;

//------------------------------------------------------------------------------

            virtual ~IWG_Maxwell_Old() = default ;

//------------------------------------------------------------------------------

            /**
             * save data to HDF5 field
             */
             void
             save( const string & aPath );

//------------------------------------------------------------------------------

            /**
             * load data from HDF5 field
             */
             int
             load( const string & aPath, Mesh * aMesh );

//------------------------------------------------------------------------------

            void
            initialize() ;

//------------------------------------------------------------------------------

            /**
             * makes sure that the mesh is sane
             */
            int
            check_mesh( Mesh * aMesh, const proc_t aMasterRank=0 );

//------------------------------------------------------------------------------

            void
            link_to_group( Group * aGroup );

//------------------------------------------------------------------------------

            virtual void
            link_jacobian_function( Group * aGroup ) ;

//------------------------------------------------------------------------------

            void
            set_blocks(
                    const Vector< id_t >     & aBlockIDs,
                    const Cell< DomainType  > & aBlockTypes );

//------------------------------------------------------------------------------

            void
            set_sidesets(
                    const Vector< id_t >      & aSideSetIDs,
                    const Cell< DomainType  > & aSideSetTypes );

//------------------------------------------------------------------------------

            /**
             * for debugging
             */
            void
            print_dofs( Element * aElement );

//------------------------------------------------------------------------------

            /**
             * this function computes the matrices M and K for the
             * superconducting part, assuming that we have edge functions
             */
            void
            compute_sc_matrices( Element * aElement,
                         Matrix< real > & aM,
                         Matrix< real > & aK );

//------------------------------------------------------------------------------

            /**
             * tell if this IWG has edge based dofs
             */
            bool
            has_edge_dofs() const ;

//------------------------------------------------------------------------------

            void
            set_current_blocks( const Vector< id_t > & aCurrentBlocks );

//------------------------------------------------------------------------------

            const Vector< id_t > &
            current_blocks() const ;

//------------------------------------------------------------------------------

            virtual void
            compute_element_current_2d( Element * aElement, real & aJz );
//------------------------------------------------------------------------------

            virtual void
            compute_element_current_3d( Element * aElement, real & aJx, real & aJy, real & aJz );


//------------------------------------------------------------------------------

            void
            interface_ha_2d( Element * aElement,
                             Matrix< real > & aM, Matrix< real > & aK );

            void
            interface_phia_2d( Element * aElement,
                            Matrix< real > & aM, Matrix< real > & aK );

            void
            interface_hphi_2d( Element * aElement, Matrix< real > & aK, Vector< real > & aRHS );


            void
            set_nan_values();

//------------------------------------------------------------------------------

            void
            add_boundary_condition( BoundaryCondition * aBoundaryCondition );

//------------------------------------------------------------------------------

            void
            update_nabla( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            inline real
            domain_increment() const ;

//------------------------------------------------------------------------------

            /** computes the help matrix for
             *  j = curl h
             */
            const Matrix< real > &
            curl_h( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------

            // edge shape function based on integration point index
            const Matrix< real > &
            E( Element * aElement, const uint aIndex );

            // edge shape function based on node coordinates in 2D
            const Matrix< real > &
            E( Element * aElement, const real aXi, const real aEta );

            // edge shape function based on node coordinates in 3D
            const Matrix< real > &
            E( Element * aElement, const real aXi, const real aEta, const real aZeta );


//------------------------------------------------------------------------------

            // computes the help matrix for
            // b = curl a
            const Matrix< real > &
            curl_a( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            // computes the help matrix for
            // h = -grad phi
            const Matrix< real > &
            grad_phi( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            // computes the normal
            const Vector< real > &
            normal( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * project integration point coordinates from a facet to an element,
             * eg. line to triangle or triangle to tetrahedron
             */
            void
            project_intpoints(
                    Element * aElement,
                          Matrix< real > & aIntpointsOnMaster,
                          Matrix< real > & aIntPointsOnSlave );

//------------------------------------------------------------------------------

            void
            compute_Nxn_2d(
                    const Matrix< real > & aXi,
                    const uint             aIntpoint,
                    const uint             aOffset,
                    const Vector< real > & aNormal,
                          Matrix< real > & aN,
                          Matrix< real > & aNxn );
//------------------------------------------------------------------------------

            void
            compute_B_2d(
                    const uint       aOffset,
                    Matrix< real > & aB );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_2d( Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_3d( Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            projection_jacobian_b2d( Element * aElement, Matrix< real> & aJacobian );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_dof_fields();

//------------------------------------------------------------------------------

             void
             collect_doftypes_per_block_and_sideset();

//------------------------------------------------------------------------------

            void
            nabla_tri3( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            void
            nabla_tet4( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            E_tri3( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            E_tet4( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            E_2d_tri3( Element * aElement, const real aXi, const real aEta );

//------------------------------------------------------------------------------

            const Matrix< real > &
            E_3d_tet4( Element * aElement, const real aXi, const real aEta, const real aZeta  );

//------------------------------------------------------------------------------

            const Matrix< real > &
            curl_h_tri3( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            curl_h_tet4( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            curl_a_tri3( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            curl_a_tet4( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            grad_phi_tri3( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Matrix< real > &
            grad_phi_tet4( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_linear( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            void
            project_intpoints_tri(
                    Element * aElement,
                          Matrix< real > & aIntpointsOnMaster,
                          Matrix< real > & aIntPointsOnSlave );

//------------------------------------------------------------------------------

            void
            project_intpoints_tet(
                    Element * aElement,
                          Matrix< real > & aIntpointsOnMaster,
                          Matrix< real > & aIntPointsOnSlave );

//------------------------------------------------------------------------------

            void
            sc_matrices_const(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            sc_matrices_powerlaw_ej(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            sc_matrices_powerlaw_ejt(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            sc_matrices_powerlaw_ejb(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            sc_matrices_powerlaw_ejbt(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------


        inline bool
        IWG_Maxwell_Old::has_edge_dofs() const
        {
            return mUseEdges ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell_Old::domain_increment() const
        {
            return mDomainIncrement ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_Old::update_nabla( Element * aElement, const uint aIndex )
        {
            ( this->*mFunNabla )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E( Element * aElement, const uint aIndex )
        {
            return ( this->*mFunE )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E( Element * aElement, const real aXi, const real aEta )
        {
            return ( this->*mFunE2D )( aElement, aXi, aEta );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E( Element * aElement, const real aXi, const real aEta, const real aZeta )
        {
            return ( this->*mFunE3D )( aElement, aXi, aEta, aZeta );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::curl_h( Element * aElement, const uint aIndex )
        {
            return ( this->*mFunCurlH )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::curl_a( Element * aElement, const uint aIndex )
        {
            return ( this->*mFunCurlA )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::grad_phi( Element * aElement, const uint aIndex )
        {
            return ( this->*mFunGradPhi )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell_Old::normal( Element * aElement, const uint aIndex )
        {
            return ( this->*mFunNormal )( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_Old::compute_sc_matrices( Element * aElement,
                                  Matrix< real > & aM,
                                  Matrix< real > & aK )
        {
            ( this->*mFunScMatrices )( aElement, aM, aK );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_Old::project_intpoints(
                Element * aElement,
                Matrix< real > & aIntpointsOnMaster,
                Matrix< real > & aIntPointsOnSlave )
        {
            ( this->*mFunProjectIntpoints )(
                    aElement,
                    aIntpointsOnMaster,
                    aIntPointsOnSlave ) ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E_tri3( Element * aElement, const uint aIndex )
        {
            const Matrix< real > & tPoints = mGroup->integration_points();

            return this->E_2d_tri3( aElement,
                                 tPoints( 0, aIndex ),
                                 tPoints( 1, aIndex ) );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E_2d_tri3( Element * aElement, const real aXi, const real aEta )
        {
            real tZeta = 1. - aXi - aEta ;

            // get the work matrix
            Matrix< real > & aE = mGroup->work_E();

            // populate data
            if( aElement->edge_direction( 0 ) )
            {
                aE( 0, 0 ) =  aXi * mNablaEta( 0 ) - aEta * mNablaXi( 0 );
                aE( 1, 0 ) =  aXi * mNablaEta( 1 ) - aEta * mNablaXi( 1 );
            }
            else
            {
                aE( 0, 0 ) =  aEta * mNablaXi( 0 ) - aXi * mNablaEta( 0 );
                aE( 1, 0 ) =  aEta * mNablaXi( 1 ) - aXi * mNablaEta( 1 ) ;
            }

            if( aElement->edge_direction( 1 ) )
            {
                aE( 0, 1 ) = aEta * mNablaZeta( 0 ) - tZeta * mNablaEta( 0 );
                aE( 1, 1 ) = aEta * mNablaZeta( 1 ) - tZeta * mNablaEta( 1 );
            }
            else
            {
                aE( 0, 1 ) = tZeta * mNablaEta( 0 ) - aEta * mNablaZeta( 0 ) ;
                aE( 1, 1 ) = tZeta * mNablaEta( 1 ) - aEta * mNablaZeta( 1 ) ;
            }

            if( aElement->edge_direction( 2 ) )
            {
                aE( 0, 2 ) = tZeta * mNablaXi( 0 ) - aXi * mNablaZeta( 0 ) ;
                aE( 1, 2 ) = tZeta * mNablaXi( 1 ) - aXi * mNablaZeta( 1 ) ;
            }
            else
            {
                aE( 0, 2 ) = aXi * mNablaZeta( 0 ) - tZeta * mNablaXi( 0 ) ;
                aE( 1, 2 ) = aXi * mNablaZeta( 1 ) - tZeta * mNablaXi( 1 ) ;
            }
            return aE ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E_tet4( Element * aElement, const uint aIndex )
        {
            const Matrix< real > & tPoints = mGroup->integration_points();

            return this->E_3d_tet4( aElement,
                                 tPoints( 0, aIndex ),
                                 tPoints( 1, aIndex ),
                                 tPoints( 2, aIndex ) );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        IWG_Maxwell_Old::E_3d_tet4( Element * aElement, const real aXi, const real aEta, const real aZeta )
        {
            real tTau  = 1. - aXi - aEta - aZeta ;

            // get work matrix
            Matrix< real > & aE = mGroup->work_E();

            // populate data
            if( aElement->edge_direction( 0 ) )
            {
                aE( 0, 0 ) = aXi * mNablaEta( 0 ) - aEta * mNablaXi( 0 ) ;
                aE( 1, 0 ) = aXi * mNablaEta( 1 ) - aEta * mNablaXi( 1 ) ;
                aE( 2, 0 ) = aXi * mNablaEta( 2 ) - aEta * mNablaXi( 2 ) ;
            }
            else
            {
                aE( 0, 0 ) = aEta * mNablaXi( 0 ) - aXi * mNablaEta( 0 ) ;
                aE( 1, 0 ) = aEta * mNablaXi( 1 ) - aXi * mNablaEta( 1 ) ;
                aE( 2, 0 ) = aEta * mNablaXi( 2 ) - aXi * mNablaEta( 2 ) ;
            }
            if( aElement->edge_direction( 1 ) )
            {
                aE( 0, 1 ) = aEta * mNablaZeta( 0 ) - aZeta * mNablaEta( 0 ) ;
                aE( 1, 1 ) = aEta * mNablaZeta( 1 ) - aZeta * mNablaEta( 1 ) ;
                aE( 2, 1 ) = aEta * mNablaZeta( 2 ) - aZeta * mNablaEta( 2 ) ;
            }
            else
            {
                aE( 0, 1 ) = aZeta * mNablaEta( 0 ) - aEta * mNablaZeta( 0 )  ;
                aE( 1, 1 ) = aZeta * mNablaEta( 1 ) - aEta * mNablaZeta( 1 )  ;
                aE( 2, 1 ) = aZeta * mNablaEta( 2 ) - aEta * mNablaZeta( 2 )  ;
            }
            if( aElement->edge_direction( 2 ) )
            {
                aE( 0, 2 ) = aZeta * mNablaXi( 0 ) - aXi * mNablaZeta( 0 ) ;
                aE( 1, 2 ) = aZeta * mNablaXi( 1 ) - aXi * mNablaZeta( 1 ) ;
                aE( 2, 2 ) = aZeta * mNablaXi( 2 ) - aXi * mNablaZeta( 2 ) ;
            }
            else
            {
                aE( 0, 2 ) = aXi * mNablaZeta( 0 ) - aZeta * mNablaXi( 0 ) ;
                aE( 1, 2 ) = aXi * mNablaZeta( 1 ) - aZeta * mNablaXi( 1 ) ;
                aE( 2, 2 ) = aXi * mNablaZeta( 2 ) - aZeta * mNablaXi( 2 ) ;
            }
            if( aElement->edge_direction( 3 ) )
            {
                aE( 0, 3 ) = aXi * mNablaTau( 0 ) - tTau * mNablaXi( 0 ) ;
                aE( 1, 3 ) = aXi * mNablaTau( 1 ) - tTau * mNablaXi( 1 ) ;
                aE( 2, 3 ) = aXi * mNablaTau( 2 ) - tTau * mNablaXi( 2 ) ;
            }
            else
            {
                aE( 0, 3 ) = tTau * mNablaXi( 0 ) - aXi * mNablaTau( 0 ) ;
                aE( 1, 3 ) = tTau * mNablaXi( 1 ) - aXi * mNablaTau( 1 ) ;
                aE( 2, 3 ) = tTau * mNablaXi( 2 ) - aXi * mNablaTau( 2 ) ;
            }
            if( aElement->edge_direction( 4 ) )
            {
                aE( 0, 4 ) = aEta * mNablaTau( 0 ) - tTau * mNablaEta( 0 ) ;
                aE( 1, 4 ) = aEta * mNablaTau( 1 ) - tTau * mNablaEta( 1 ) ;
                aE( 2, 4 ) = aEta * mNablaTau( 2 ) - tTau * mNablaEta( 2 ) ;
            }
            else
            {
                aE( 0, 4 ) = tTau * mNablaEta( 0 ) - aEta * mNablaTau( 0 ) ;
                aE( 1, 4 ) = tTau * mNablaEta( 1 ) - aEta * mNablaTau( 1 ) ;
                aE( 2, 4 ) = tTau * mNablaEta( 2 ) - aEta * mNablaTau( 2 ) ;
            }
            if( aElement->edge_direction( 5 ) )
            {
                aE( 0, 5 ) = aZeta * mNablaTau( 0 ) - tTau * mNablaZeta( 0 ) ;
                aE( 1, 5 ) = aZeta * mNablaTau( 1 ) - tTau * mNablaZeta( 1 ) ;
                aE( 2, 5 ) = aZeta * mNablaTau( 2 ) - tTau * mNablaZeta( 2 ) ;
            }
            else
            {
                aE( 0, 5 ) = tTau * mNablaZeta( 0 ) - aZeta * mNablaTau( 0 ) ;
                aE( 1, 5 ) = tTau * mNablaZeta( 1 ) - aZeta * mNablaTau( 1 ) ;
                aE( 2, 5 ) = tTau * mNablaZeta( 2 ) - aZeta * mNablaTau( 2 ) ;
            }

            return aE ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_Old::set_current_blocks( const Vector< id_t > & aCurrentBlocks )
        {
            mCurrentBlocks = aCurrentBlocks ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        IWG_Maxwell_Old::current_blocks() const
        {
            return mCurrentBlocks ;
        }

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_IWG_MAXWELL_OLD_HPP
