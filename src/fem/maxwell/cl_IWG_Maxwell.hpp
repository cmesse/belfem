//
// Created by christian on 12/14/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPP
#define BELFEM_CL_IWG_MAXWELL_HPP
// #define BELFEM_FERROAIR_ENRICHED --> moved to IWG

//#define BELFEM_FERRO_LINEAR
#define BELFEM_FERRO_HPHIA // todo: remove this

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_trans.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
#include "en_SolverEnums.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_FEM_DomainGroup.hpp"
#include "cl_FEM_BoundaryCondition.hpp"
#include "en_FEM_MagfieldBcType.hpp"
#include "cl_EF_EdgeFunction.hpp"
#include "cl_Maxwell_FieldList.hpp"

#include "meshtools.hpp"
#include "en_FEM_NedelecField.hpp"

namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        /**
         * an abstract baseclass for
         * solving the Maxwell equations
         */
        class IWG_Maxwell : public IWG_Timestep
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! element type of block elements
            const ElementType mElementType ;

            //! contains the number of spatial dimensions
            const uint mNumberOfDimensions;

            //! flag telling if this IWG uses edge dofs
            const bool mUseEdges;

            //! list with dofs per entity
            maxwell::FieldList mFields ;

            //! the edge shape function
            EdgeFunction * mEdgeFunction = nullptr ;


            //! edge function for thin shell
            EdgeFunction * mEdgeFunctionTS = nullptr ;


            //! map telling which subtype a domain boundary is
            Map< id_t, MagfieldBcType > mMagfieldTypeMap ;

            //! map telling the gost indices of thin shell elements
            Map< luint, mesh::Element * > mGhostElementMap ;

            Vector< real > mQ0cut = { 0, 0, 0 };

            // this vector is needed for the error estimation in the thin shell mode
            Vector< uint > mNumShellsPerNode ;

            Vector< real > mWorkCurrent = { 0, 0, 0 };
            Vector< real > mWorkCurrentK = { 0, 0, 0 };

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // link to function that collects element data
            void
            ( IWG_Maxwell::*mCollectNedelecData )
            (             Element   * aElement,
                    const NedelecField   aNedelecField,
                     Vector< real > & aNedelecData ) ;

            // link to function that computes superconducting matrices
            void
            ( IWG_Maxwell::*mComputeScMatrices )
            (       Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            // link to function that cuputes superconducting matrices
            void
            ( IWG_Maxwell::*mComputeFerroMatrix )
            (       Element        * aElement,
                    Matrix< real > & aK );

            // integration data for slave element on interface
            // this function is different for triangles and tets
            const IntegrationData *
            ( IWG_Maxwell::*mGetInterfaceData )( Element * aElement );

            // these cells must be consistent with NedelecField. Do not change!
            Cell < string > mNedelecFields = { "edge_h", "edge_h0", "shell_edge_h", "shell_edge_h0" };
            Cell < string > mFaceFields = { "face_h", "face_h0", "face_edge_h", "face_edge_h0" };


            // current evaluation
            real
            ( IWG_Maxwell::*mComputeCurrent )(
                    Element * aElement,
                    const uint aIndex
                    );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell(
                    const ElementType aElementType,
                    const IwgType aType,
                    const IwgMode aMode,
                    const SymmetryMode aSymmetryMode,
                    const SideSetDofLinkMode SideSetDofLinkMode,
                    const bool         aUseEdges = false ) ;

//------------------------------------------------------------------------------

            virtual ~IWG_Maxwell();

//------------------------------------------------------------------------------

            void
            initialize() ;

//------------------------------------------------------------------------------

            /**
             * links the IWG to a group.
             * Must be called before computing matrices.
             * @param aGroup
             */
            void
            link_to_group( Group * aGroup );

            void
            allocate_work_matrices( Group * aGroup );

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            // todo: remove this!
            virtual real
            compute_element_current( Element * aElement );

//------------------------------------------------------------------------------

            bool
            has_edge_dofs() const ;

//------------------------------------------------------------------------------

            // copy dofs into fields from last timestep
            void
            shift_fields();

//------------------------------------------------------------------------------

            /**
             * load data from HDF5 field
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

            /**
             * makes sure that the block and sideset types in the maps
             * are also written into the groups
             */
             void
             update_group_types();

//------------------------------------------------------------------------------

            /**
             * write NaN where node fields are not defined
             */
             void
             set_nan_values();

//------------------------------------------------------------------------------

            // called by factory to create the index map for the thin shells
            void
            create_ghost_map(  Mesh * aMesh,
                              const Vector< id_t >  & aThinShellSideSets,
                              const Vector< proc_t > & aCommTable ) ;

//------------------------------------------------------------------------------

            /**
              * called by factory to specialize Magfield BC
              */
             void
             set_magfield_bc_type( const id_t aSideSetID, const MagfieldBcType aType );

//------------------------------------------------------------------------------

            /**
              * called by factory to support error estimator
              */
            Vector< uint > &
            num_shells_per_node() ;

//------------------------------------------------------------------------------

            void
            compute_thin_shell_error_2d();

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * The function is selected based on the chosen group.
             * Can only be done by child class
             * @param aGroup
             */
            virtual void
            link_jacobian_function( Group * aGroup ) ;

//------------------------------------------------------------------------------

            /**
             * collect the edge data
             */
             void
             collect_nedelec_data( Element        * aElement,
                                   const NedelecField  aNedelecField,
                                   Vector< real > & aNedelecData );

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs_superconductor(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            /**
             * computes the curl of a-field
             * @param aIndex    index of integration point
             * @param aA        matrix populated with [ ax, ay, az ]
             * @return
             */
            inline const Vector< real > &
            curl_a( const uint aIndex, const Matrix< real > & aA );

//------------------------------------------------------------------------------

            /**
             * returns the vector [ ax0, ay0, az0, ax1, ay1, az1, ... ]
             */
            inline const Vector< real > &
            last_a( Element        * aElement );

//------------------------------------------------------------------------------

            const IntegrationData *
            interface_data( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_straight_2d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_curved_2d( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_straight_3d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            normal_curved_3d( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_hphi_2d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_hphi_3d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_ha_2d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_ha_3d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_aphi_2d( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_aphi_3d( Element * aElement );

//------------------------------------------------------------------------------

            void
            compute_interface_ha_tri3(
                    Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_interface_aa_tri3(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_interface_aa_tri6(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_interface_ha_tri6(
                    Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_interface_ha_tri6_tri3(
                    Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_interface_ha_tet4(
                    Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_interface_ha_tet10(
                    Element        * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_air_phi_linear(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_air_phi_higher_order(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs_cut( Element        * aElement,
                                          Matrix< real > & aJacobian,
                                          Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmmetry_a_tri3( Element        * aElement,
                                          Matrix< real > & aJacobian,
                                          Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmmetry_a_tri3( Element        * aElement,
                                                       Matrix< real > & aJacobian,
                                                       Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmmetry_a_tri6( Element        * aElement,
                                                       Matrix< real > & aJacobian,
                                                       Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmmetry_a_tri6( Element        * aElement,
                                                           Matrix< real > & aJacobian,
                                                           Vector< real > & aRHS );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            link_collect_nedelec_data_function( Group * aGroup ) ;

//------------------------------------------------------------------------------

            void
            collect_nedelec_data_linear(
                    Element   * aElement,
                    const NedelecField   aNedelecField,
                    Vector< real > & aNedelecData );

//------------------------------------------------------------------------------

            void
            collect_nedelec_data_qaudratic(
                    Element   * aElement,
                    const NedelecField   aNedelecField,
                    Vector< real > & aNedelecData );

//------------------------------------------------------------------------------

            void
            link_sc_matrix_function( Group * aGroup ) ;

//------------------------------------------------------------------------------

            void
            link_ferro_matrix_function( Group * aGroup ) ;

//------------------------------------------------------------------------------

            void
            sc_matrices_const_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );


//------------------------------------------------------------------------------

            void
            sc_matrices_j_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jt_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jb_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jbt_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_t_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_b_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_bt_first_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            write_element_data_to_mesh(
                    Element * aElement,
                    const real aJ,
                    const real aJJcrit,
                    const real aRho,
                    const real aEJ );

//------------------------------------------------------------------------------

            void
            sc_matrices_const_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_j_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jt_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jb_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_jbt_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_t_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_b_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

            void
            sc_matrices_bt_higher_order(
                    Element * aElement,
                    Matrix< real > & aM,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_tri3(
                    Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_tri6(
                    Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_tet4(
                    Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            ferro_stiffness_tet10(
                    Element * aElement, Matrix< real > & aK );

//------------------------------------------------------------------------------

            const IntegrationData *
            integration_data_tri( Element * aElement );

//------------------------------------------------------------------------------

            const IntegrationData *
            integration_data_tet( Element * aElement );

//------------------------------------------------------------------------------

            real
            compute_current( Element * aElement, const uint aIndex=0 );

//------------------------------------------------------------------------------

            real
            compute_current_2d_first_order( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            real
            compute_current_3d_first_order( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            real
            compute_current_2d_higher_order( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------

            real
            compute_current_3d_higher_order( Element * aElement, const uint aIndex );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline bool
        IWG_Maxwell::has_edge_dofs() const
        {
            return mUseEdges ;
        }

//-----------------------------------------------------------------------------

        inline void
        IWG_Maxwell::collect_nedelec_data(
        Element   * aElement,
        const NedelecField   aNedelecField,
            Vector< real > & aNedelecData )
        {
            (this->*mCollectNedelecData )
                ( aElement, aNedelecField, aNedelecData );
        }


//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell::collect_nedelec_data_linear(
                           Element   * aElement,
                        const NedelecField   aNedelecField,
                        Vector< real > & aNedelecData )
        {
            this->collect_edge_data( aElement,
                                     mNedelecFields( static_cast< uint > ( aNedelecField ) ),
                                     aNedelecData );
       }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell::collect_nedelec_data_qaudratic(
                Element   * aElement,
                const NedelecField   aNedelecField,
                Vector< real > & aNedelecData )
        {

            this->collect_edge_data( aElement,
                                     mNedelecFields( static_cast< uint > ( aNedelecField ) ),
                                     mFaceFields( static_cast< uint > ( aNedelecField ) ),
                                     aNedelecData );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::curl_a( const uint aIndex, const Matrix< real > & aA )
        {
            const Matrix< real > & tB = mEdgeFunction->B( aIndex );

            Vector< real > & aB = mGroup->work_chi() ;

            // bx = dAz/dy - dAy/dz
            aB( 0 ) =     dot( tB.row( 1  ), aA.col( 2 ) )
                             -  dot( tB.row( 2  ), aA.col( 1 ) );

            // by = dAx/dz - dAz/dx
            aB( 1 ) =      dot( tB.row( 2  ), aA.col( 0 ) )
                              -  dot( tB.row( 0  ), aA.col( 2 ) );

            // bz = dAy/dx - dAx/dy
            aB( 2 ) =      dot( tB.row( 0  ), aA.col( 1 ) )
                              -  dot( tB.row( 1  ), aA.col( 0 ) );

            return aB ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::last_a( Element        * aElement )
        {
            // get data from last timestep
            Matrix< real > & tA = mGroup->work_Phi() ;
            this->collect_node_data( aElement, mFields.FerroLast, tA );

            // assemble coordinate vector
            Vector< real > & aQ = mGroup->work_phi() ;

            uint tCount = 0 ;
            for( uint k=0; k<mNumberOfNodesPerElement; ++k )
            {
                for( uint i=0; i<3; ++i )
                {
                    aQ( tCount++ ) = tA( k, i );
                }
            }

            return aQ ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // compute the mass and stiffness matrices
            ( this->*mComputeScMatrices )(
                    aElement,
                    aJacobian,           // work matrix for M
                    mGroup->work_K()     // work matrix for K
                    );

            // get data from last timestep
            this->collect_nedelec_data( aElement,
                                        NedelecField::H0,
                                        mGroup->work_nedelec() );

            // compute the right hand side
            aRHS = aJacobian * mGroup->work_nedelec() ;

            // finalize the Jacobian
            aJacobian += mDeltaTime * mGroup->work_K() ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // compute the stiffness
            ( this->*mComputeFerroMatrix )( aElement, aJacobian );

            // scale stiffness matrix
            aJacobian *= mDeltaTime ;

            // compute the right hand side, is zero since M=0 and theta=0
            aRHS.fill( 0.0 );
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        IWG_Maxwell::interface_data( Element * aElement )
        {
            // link edge function with element
            mEdgeFunction->link( aElement->master() );

            // shape functions on edge element
            mEdgeFunction->precompute( mGroup->master_integration(
                    aElement->facet()->master_index() )->points() );

            // also populate node coordinates
            aElement->master()->get_node_coors( mGroup->work_Xm() );
            aElement->slave()->get_node_coors( mGroup->work_Xs() );
            aElement->get_node_coors( mGroup->work_X() );

            // return the integration function for the nodes
            return ( this->*mGetInterfaceData )( aElement );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_straight_2d( Element * aElement )
        {
            return this->normal_tri3( aElement );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_curved_2d( Element * aElement, const uint aIndex )
        {
            return this->normal_tri6( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_straight_3d( Element * aElement )
        {
            BELFEM_ASSERT( mesh::geometry_type( aElement->element()->type() ) == GeometryType::TRI,
                          "Element must be a triangle" );

            // get the normal vector
            Vector< real > & aN = mNormal3D ;


            aN( 0 ) =
              (aElement->element()->node( 0 )->y() - aElement->element()->node( 1 )->y() )
            * (aElement->element()->node( 0 )->z() - aElement->element()->node( 2 )->z() )
            - (aElement->element()->node( 0 )->y() - aElement->element()->node( 2 )->y() )
            * (aElement->element()->node( 0 )->z() - aElement->element()->node( 1 )->z() );


            aN( 1 ) =
               (aElement->element()->node( 0 )->x() - aElement->element()->node( 2 )->x() )
             * (aElement->element()->node( 0 )->z() - aElement->element()->node( 1 )->z() )
             - (aElement->element()->node( 0 )->x() - aElement->element()->node( 1 )->x() )
             * (aElement->element()->node( 0 )->z() - aElement->element()->node( 2 )->z() );


            aN( 2 ) = (aElement->element()->node( 0 )->x() - aElement->element()->node( 1 )->x() )
                    * (aElement->element()->node( 0 )->y() - aElement->element()->node( 2 )->y() )
                    - (aElement->element()->node( 0 )->x() - aElement->element()->node( 2 )->x() )
                    * (aElement->element()->node( 0 )->y() - aElement->element()->node( 1 )->y() ) ;

            mGroup->work_det_J() = norm( aN );
            aN /= mGroup->work_det_J() ;

            return aN ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_curved_3d( Element * aElement, const uint aIndex )
        {

            BELFEM_ERROR( false, "rewrite this function");

            BELFEM_ASSERT( mesh::geometry_type( aElement->element()->type() ) == GeometryType::TRI,
                          "Element must be a triangle" );

            Vector< real > & aN = mNormal3D ;

            Matrix< real > & tJ = mGroup->work_L() ;
            tJ = mGroup->dNdXi( aIndex ) * mGroup->work_X() ;

			// normal
            aN( 0 ) =   tJ( 0, 1 ) * tJ( 1, 2 )
                      - tJ( 0, 2 ) * tJ( 1, 1 );

            aN( 1 ) =   tJ( 0, 2 ) * tJ( 1, 0 )
                      - tJ( 0, 0 ) * tJ( 1, 2 );

            aN( 2 ) =   tJ( 0, 0 ) * tJ( 1, 1 )
                      - tJ( 0, 1 ) * tJ( 1, 0 );

            mGroup->work_det_J() = norm( aN );
            aN /= mGroup->work_det_J() ;
            return  aN ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        IWG_Maxwell::integration_data_tri( Element * aElement )
        {
            return mGroup->slave_integration( aElement->facet()->slave_index() );
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        IWG_Maxwell::integration_data_tet( Element * aElement )
        {
            return mGroup->slave_integration(
                    aElement->facet()->slave_index() * 3
                + aElement->facet()->orientation_on_slave()  );
        }

//------------------------------------------------------------------------------

        inline Vector< uint > &
        IWG_Maxwell::num_shells_per_node()
        {
            return mNumShellsPerNode ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell::compute_current( Element * aElement, const uint aIndex )
        {
            return ( this->*mComputeCurrent ) ( aElement, aIndex );
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell::compute_current_2d_first_order( Element * aElement, const uint aIndex  )
        {
            real aJz = dot( mEdgeFunction->C( aIndex ).row( 0 ), mGroup->work_nedelec() );
            mMesh->field_data( "elementJz")( aElement->element()->index() )  = aJz ;

            return aJz ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell::compute_current_3d_first_order( Element * aElement, const uint aIndex  )
        {
            mWorkCurrent = mEdgeFunction->C( aIndex ).row( 0 ) * mGroup->work_nedelec();

            mMesh->field_data( "elementJx" )( aElement->element()->index()) = mWorkCurrent( 0 );
            mMesh->field_data( "elementJy" )( aElement->element()->index()) = mWorkCurrent( 1 );
            mMesh->field_data( "elementJz" )( aElement->element()->index()) = mWorkCurrent( 2 );

            return norm( mWorkCurrent );
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell::compute_current_2d_higher_order( Element * aElement, const uint aIndex  )
        {
            mWorkCurrentK( 2 ) = dot( mEdgeFunction->C( aIndex ).row( 0 ), mGroup->work_nedelec() );

            if( aIndex == 0 )
            {
                mWorkCurrent( 2 ) = mGroup->integration_weights()( aIndex ) * mWorkCurrentK( 2 ) ;
            }
            else if( aIndex == mNumberOfIntegrationPoints )
            {
                mWorkCurrent( 2 ) += mGroup->integration_weights()( aIndex ) * mWorkCurrentK( 2 ) ;
                mWorkCurrent( 2 ) /= mEdgeFunction->sum_w() ;

                mMesh->field_data( "elementJz")( aElement->element()->index() )  = mWorkCurrent( 2 ) ;
            }
            else
            {
                mWorkCurrent( 2 ) += mGroup->integration_weights()( aIndex ) * mWorkCurrentK( 2 ) ;
            }

            return mWorkCurrent( 2 );
        }

//------------------------------------------------------------------------------

        inline real
        IWG_Maxwell::compute_current_3d_higher_order( Element * aElement, const uint aIndex  )
        {
            mWorkCurrentK = mEdgeFunction->C( aIndex ).row( 0 ) * mGroup->work_nedelec() * mGroup->integration_weights()( aIndex );

            if( aIndex == 0 )
            {
                mWorkCurrent = mGroup->integration_weights()( aIndex ) * mWorkCurrentK ;
            }
            else if( aIndex == mNumberOfIntegrationPoints )
            {
                mWorkCurrent +=  mGroup->integration_weights()( aIndex ) * mWorkCurrentK ;
                mWorkCurrent /= mEdgeFunction->sum_w() ;

                mMesh->field_data( "elementJx")( aElement->element()->index() )  = mWorkCurrent( 0 ) ;
                mMesh->field_data( "elementJy")( aElement->element()->index() )  = mWorkCurrent( 1 ) ;
                mMesh->field_data( "elementJz")( aElement->element()->index() )  = mWorkCurrent( 2 ) ;
            }
            else
            {
                mWorkCurrent +=  mGroup->integration_weights()( aIndex ) * mWorkCurrentK ;
            }

            return norm( mWorkCurrentK );
        }

//------------------------------------------------------------------------------
    }
}
//------------------------------------------------------------------------------
#endif //BELFEM_CL_IWG_MAXWELL_HPP
