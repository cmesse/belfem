//
// Created by christian on 12/14/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPP
#define BELFEM_CL_IWG_MAXWELL_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

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

namespace belfem
{
    namespace fem
    {
        enum class NedelecField
        {
            H         = 0 ,
            H0        = 1 ,
            TSH       = 2 ,
            TSH0      = 3 ,
            UNDEFINED = 4
        };

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

            //! flag telling if this IWG has sheet dofs
            //! determined during initialization
            bool mHaveShellDofs = false ;

            //! the edge shape function
            EdgeFunction * mEdgeFunction = nullptr ;

            //! map telling which subtype a domain boundary is
            Map< id_t, MagfieldBcType > mMagfieldTypeMap ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            // link to function that collects element data
            void
            ( IWG_Maxwell::*mCollectNedelecData )
            (             Element   * aElement,
                    const NedelecField   aNedelecField,
                     Vector< real > & aNedelecData ) ;

            // link to function that cuputes superconducting matrices
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

            Vector< real > mNormal2D = { 0., 0., };
            Vector< real > mNormal3D = { 0., 0., 0. };

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
            compute_interface_ha_tri6(
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
        private:
//------------------------------------------------------------------------------

            void
            link_collect_nedelec_data_function( Group * aGroup ) ;

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
            aElement->slave()->get_node_coors( mGroup->node_coords() );

            // these data are needed for the normal
            if( aElement->element()->is_curved() )
            {
                aElement->get_node_coors( mGroup->work_X() );
            }

            // return the integration function for the nodes
            return ( this->*mGetInterfaceData )( aElement );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_straight_2d( Element * aElement )
        {

            BELFEM_ASSERT( mesh::geometry_type( aElement->element()->type() ) == GeometryType::LINE,
                          "Element must be a line" );

            // get the normal vector
            Vector< real > & aN = mNormal2D ;
            aN( 0 ) = aElement->element()->node( 1 )->y()  - aElement->element()->node( 0 )->y() ;
            aN( 1 ) = aElement->element()->node( 0 )->x()  - aElement->element()->node( 1 )->x() ;

            mGroup->work_det_J() = norm( aN );
            aN /= mGroup->work_det_J() ;

            // scale determinant
            mGroup->work_det_J() *= 0.5 ;

            return aN ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell::normal_curved_2d( Element * aElement, const uint aIndex )
        {
            Vector< real > & aN = mNormal2D ;

            aN.fill( 0.0 );

            // the direction vector
            aN( 0 ) = dot( mGroup->dNdXi( aIndex ).row( 0 ),
                                  mGroup->work_X().col( 1 ) );

            aN( 1 ) = -dot( mGroup->dNdXi( aIndex ).row( 0 ),
                            mGroup->work_X().col( 0 ) );

            // scale normal and integration scale
            mGroup->work_det_J() = norm( aN );

            aN /= mGroup->work_det_J() ;

            return aN ;
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
    }
}
//------------------------------------------------------------------------------
#endif //BELFEM_CL_IWG_MAXWELL_HPP
