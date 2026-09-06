/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_FEM_POSTPROCESSOR_HPP
#define BELFEM_CL_FEM_POSTPROCESSOR_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
        class Postprocessor
        {

        protected:
            const proc_t mCommRank ;
            const proc_t mCommSize ;
            const uint   mNumDimensions ;
            const uint   mOrder ;
        private:

            // block IDs used on this post processor
            Cell< id_t > mBlockIDs ;

            Cell< index_t > mMyNodeIndices ;
            Cell< index_t > mMyOwnedNodeIndices ;
            Cell< index_t > mMyElementIndices ;
            Cell< index_t > mMyEdgeIndices ;
            Cell< index_t > mMyFaceIndices ;

            Cell< Vector< index_t > > mAllNodeIndices ;
            Cell< Vector< index_t > > mAllOwnedNodeIndices ;

            Cell< Vector< index_t > > mAllEdgeIndices ;
            Cell< Vector< index_t > > mAllFaceIndices ;

            // vector with polynomial coefficients
            // we use a matrix so that we can use the trans() function
            Matrix< real > mPoly ;

            DynamicBitset * mNodeBitset = nullptr ;
            DynamicBitset * mEdgeBitset = nullptr ;
            DynamicBitset * mFaceBitset = nullptr ;
            DynamicBitset * mElementBitset = nullptr ;

            void
            ( Postprocessor::*mFunComputePoly )( const Vector< real > & aX );

            void
            ( Postprocessor::*mFunPoly2D )( const real x, const real y );

            void
            ( Postprocessor::*mFunPoly3D )( const real x, const real y, const real z );

            bool mIsInitialized = false ;

            bool mHaveEdges = false ;
            bool mHaveFaces = false ;

            ElementType mElementType = ElementType::EMPTY ;
            id_t mLastElementID = gNoID ;
            id_t mLastBlockID   = gNoID ;
            Matrix< real > mVandermonde ;
            Map< mesh::Node * , Matrix< real > * > mNodeMatrices ;

            //Cell< Matrix< real > > mNodeCoefficients ;

            uint mNumCoefficients = 0 ;

            Vector< real > mNull ;

        protected:

            DomainType mDomainType = DomainType::Default ;


            Cell< string > mTargetFields ;
            Cell< string > mSourceFields ;
            uint mNumSourceFields = 0 ;
            uint mNumTargetFields = 0 ;

            Kernel     * mKernel   = nullptr;
            Mesh       * mMesh     = nullptr;
            DofManager * mField    = nullptr;
            IWG        * mEquation = nullptr;
            Material   * mMaterial = nullptr;
            Element    * mElement    = nullptr ;
            Block      * mBlock      = nullptr ;
            Calculator * mCalculator = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Postprocessor( Kernel * aKernel , DofManager * aField = nullptr );

            virtual ~Postprocessor() ;

            void
            set_field( DofManager * aField );

            void
            set_source_fields( const Cell< string > & aFields );

            void
            set_target_fields( const Cell< string > & aFields );

            const Cell< string > &
            source_fields() const ;

            const Cell< string > &
            target_fields() const ;

            void
            set_block_ids( const Vector< id_t > & aBlockIDs );

            void
            set_block_ids( const Cell< id_t > & aBlockIDs );

            virtual void
            initialize();

            void
            synch_source_fields();

            void
            synch_target_fields( Matrix< real > & aData );

            virtual void
            run();

            DofManager *
            projector() ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            const Matrix< real > &
            compute_poly( const Vector< real > & aX );

            void
            recover_fields();

            virtual const Vector< real > &
            compute( const uint aK );

            virtual void
            update_element_dofs();

            const Cell< id_t > &
            block_ids() const ;

            // needed for maxwell postproc
            // to create element fields
            const Cell< index_t > &
            my_element_indices() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_target_fields();

            void
            check_if_we_have_edges_and_faces();

            void
            select_elements_and_owned_nodes();

            void
            select_all_relevant_nodes();

            void
            synch_node_indices( const bool aOwnedOnly );

            void
            select_edges();

            void
            select_faces();

            void
            select_polynomial();

            // funciton needed by constructor
            uint
            get_interpolation_order( const Mesh * aMesh ) const;

            void
            compute_node_matrices();

            void
            compute_element_coeffs(
               mesh::Node * aNode,
               Matrix< real > & aNodeCoords,
               Map< mesh::Element*, Matrix< real > * > & aElementData );


            void
            compute_element_coeffs(
                mesh::Element * aElement,
                Matrix< real > & aNodeCoords,
                Map< mesh::Element*, Matrix< real > * > & aElementData );

            void
            check_element_coeffs_done( mesh::Element * aElement, Map< mesh::Element*, Matrix< real > * > & aElementData );

            void
            synch_source_field( const string & aField );

//------------------------------------------------------------------------------
// EVALUATION ROUTINES
//------------------------------------------------------------------------------

            void
            compute_poly_2d( const Vector< real > & aX );

            void
            compute_poly_3d( const Vector< real > & aX );

            void
            poly1_2d( const real x, const real y );

            void
            poly2_2d( const real x, const real y );

            void
            poly3_2d( const real x, const real y );

            void
            poly4_2d( const real x, const real y );

            void
            poly1_3d( const real x, const real y, const real z  );

            void
            poly2_3d( const real x, const real y, const real z  );

            void
            poly3_3d( const real x, const real y, const real z  );

            void
            poly4_3d( const real x, const real y, const real z  );
        };

        inline
        DofManager *
        Postprocessor::projector()
        {
            return mField;
        }

        inline void
        Postprocessor::set_source_fields( const Cell< string > & aFields )
        {
            mSourceFields = aFields ;
            mNumSourceFields = aFields.size() ;
        }

        inline void
        Postprocessor::set_target_fields( const Cell< string > & aFields )
        {
            mTargetFields = aFields ;
            mNumTargetFields = aFields.size() ;
        }

        inline const Cell< string > &
        Postprocessor::source_fields() const
        {
            return mSourceFields;
        }

        inline const Cell< string > &
        Postprocessor::target_fields() const
        {
            return mTargetFields;
        }

        inline const Matrix< real > &
        Postprocessor::compute_poly( const Vector< real > & aX )
        {
            (this->*mFunComputePoly)( aX );
            return mPoly;
        }

        inline void
        Postprocessor::compute_poly_2d( const Vector< real > & aX )
        {
            (this->*mFunPoly2D )( aX( 0 ), aX( 1 ) );
        }

        inline void
        Postprocessor::compute_poly_3d( const Vector< real > & aX )
        {
            (this->*mFunPoly3D )( aX( 0 ), aX( 1 ), aX( 2 ) );
        }

        inline void
        Postprocessor::poly1_2d( const real x, const real y )
        {
            mPoly( 0, 0 ) = 1.0 ;
            mPoly( 1, 0 ) = x ;
            mPoly( 2, 0 ) = y ;
        }

        inline void
        Postprocessor::poly2_2d( const real x, const real y )
        {
            mPoly( 0, 0 ) = 1.0 ;
            mPoly( 1, 0 ) = x ;
            mPoly( 2, 0 ) = y ;
            mPoly( 3, 0 ) = x*x ;
            mPoly( 4, 0 ) = x*y ;
            mPoly( 5, 0 ) = y*y ;
        }

        inline void
        Postprocessor::poly3_2d( const real x, const real y )
        {
            mPoly( 0, 0 ) = 1.0 ;
            mPoly( 1, 0 ) = x ;
            mPoly( 2, 0 ) = y ;
            mPoly( 3, 0 ) = x*x ;
            mPoly( 4, 0 ) = x*y ;
            mPoly( 5, 0 ) = y*y ;
            mPoly( 6, 0 ) = mPoly( 3, 0 )*x ;
            mPoly( 7, 0 ) = mPoly( 3, 0 )*y ;
            mPoly( 8, 0 ) = x * mPoly( 5, 0 );
            mPoly( 9, 0 ) = y * mPoly( 5, 0 );
        }

        inline void
        Postprocessor::poly4_2d( const real x, const real y )
        {
            mPoly(  0, 0 ) = 1.0 ;
            mPoly(  1, 0 ) = x ;
            mPoly(  2, 0 ) = y ;
            mPoly(  3, 0 ) = x*x ;
            mPoly(  4, 0 ) = x*y ;
            mPoly(  5, 0 ) = y*y ;
            mPoly(  6, 0 ) = mPoly( 3, 0 )*x ;  // x^3
            mPoly(  7, 0 ) = mPoly( 3, 0 )*y ;  // x^2 * y
            mPoly(  8, 0 ) = x * mPoly( 5, 0 ); // x * y^2
            mPoly(  9, 0 ) = y * mPoly( 5, 0 ); // y^3
            mPoly( 10, 0 ) = mPoly( 3, 0 ) * mPoly( 3, 0 ); // x^4
            mPoly( 11, 0 ) = mPoly( 6, 0 ) * y ; // x^3 * y
            mPoly( 12, 0 ) = mPoly( 3, 0 ) * mPoly( 5, 0 ) ; // x^2*y^2
            mPoly( 13, 0 ) = x * mPoly( 9, 0 ) ; // x * y^3
            mPoly( 14, 0 ) =  mPoly( 5, 0 ) * mPoly( 5, 0 ) ; // y^4
        }

        inline void
        Postprocessor::poly1_3d( const real x, const real y, const real z )
        {
            mPoly( 0, 0 ) = 1.0 ;
            mPoly( 1, 0 ) = x ;
            mPoly( 2, 0 ) = y ;
            mPoly( 3, 0 ) = z ;
        }

        inline void
        Postprocessor::poly2_3d( const real x, const real y, const real z )
        {
            mPoly(  0, 0 ) = 1.0 ;
            mPoly(  1, 0 ) = x ;
            mPoly(  2, 0 ) = y ;
            mPoly(  3, 0 ) = z ;
            mPoly(  4, 0 ) = x*x ;
            mPoly(  5, 0 ) = x*y ;
            mPoly(  6, 0 ) = y*y ;
            mPoly(  7, 0 ) = y*z ;
            mPoly(  8, 0 ) = z*z ;
            mPoly(  9, 0 ) = z*x ;
        }

        inline void
        Postprocessor::poly3_3d( const real x, const real y, const real z )
        {
            mPoly(  0, 0 ) = 1.0 ;
            mPoly(  1, 0 ) = x ;
            mPoly(  2, 0 ) = y ;
            mPoly(  3, 0 ) = z ;
            mPoly(  4, 0 ) = x*x ;
            mPoly(  5, 0 ) = x*y ;
            mPoly(  6, 0 ) = y*y ;
            mPoly(  7, 0 ) = y*z ;
            mPoly(  8, 0 ) = z*z ;
            mPoly(  9, 0 ) = z*x ;
            mPoly( 10, 0 ) = mPoly( 4, 0 ) * x ; // x^3
            mPoly( 11, 0 ) = mPoly( 4, 0 ) * y ; // x^2 * y
            mPoly( 12, 0 ) = x * mPoly( 6, 0 ) ; // x * y^2 ;
            mPoly( 13, 0 ) = y * mPoly( 6, 0 ) ; // y^3
            mPoly( 14, 0 ) = z * mPoly( 6, 0 ) ; // y^2 * z
            mPoly( 15, 0 ) = y * mPoly( 8, 0 ) ; // y * z^2
            mPoly( 16, 0 ) = z * mPoly( 8, 0 ) ; // z^3
            mPoly( 17, 0 ) = x * mPoly( 8, 0 ) ; // z^2 * x
            mPoly( 18, 0 ) = mPoly( 4, 0 ) * z ; // x^2 * z
            mPoly( 19, 0 ) = x * y * z ;
        }

        inline void
        Postprocessor::poly4_3d( const real x, const real y, const real z )
        {
            mPoly(  0, 0 ) = 1.0 ;
            mPoly(  1, 0 ) = x ;
            mPoly(  2, 0 ) = y ;
            mPoly(  3, 0 ) = z ;
            mPoly(  4, 0 ) = x*x ;
            mPoly(  5, 0 ) = x*y ;
            mPoly(  6, 0 ) = y*y ;
            mPoly(  7, 0 ) = y*z ;
            mPoly(  8, 0 ) = z*z ;
            mPoly(  9, 0 ) = z*x ;
            mPoly( 10, 0 ) = mPoly( 4, 0 ) * x ; // x^3
            mPoly( 11, 0 ) = mPoly( 4, 0 ) * y ; // x^2 * y
            mPoly( 12, 0 ) = x * mPoly( 6, 0 ) ; // x * y^2 ;
            mPoly( 13, 0 ) = y * mPoly( 6, 0 ) ; // y^3
            mPoly( 14, 0 ) = z * mPoly( 6, 0 ) ; // y^2 * z
            mPoly( 15, 0 ) = y * mPoly( 8, 0 ) ; // y * z^2
            mPoly( 16, 0 ) = z * mPoly( 8, 0 ) ; // z^3
            mPoly( 17, 0 ) = x * mPoly( 8, 0 ) ; // z^2 * x
            mPoly( 18, 0 ) = mPoly( 4, 0 ) * z ; // x^2 * z
            mPoly( 19, 0 ) = x * y * z ;
            mPoly( 20, 0 ) = mPoly( 4, 0 ) * mPoly( 4, 0 ) ; // x^4
            mPoly( 21, 0 ) = mPoly( 4, 0 ) * mPoly( 5, 0 ) ; // x^3 * y
            mPoly( 22, 0 ) = mPoly( 4, 0 ) * mPoly( 6, 0 ) ; // x^2 * y^2
            mPoly( 23, 0 ) = mPoly( 5, 0 ) * mPoly( 6, 0 ) ; // x * y^3
            mPoly( 24, 0 ) = mPoly( 6, 0 ) * mPoly( 6, 0 ) ; // y^4
            mPoly( 25, 0 ) = mPoly( 6, 0 ) * mPoly( 7, 0 ) ; // y^3 * z
            mPoly( 26, 0 ) = mPoly( 6, 0 ) * mPoly( 8, 0 ) ; // y^2 * z^2
            mPoly( 27, 0 ) = mPoly( 7, 0 ) * mPoly( 8, 0 ) ; // y * z^3
            mPoly( 28, 0 ) = mPoly( 8, 0 ) * mPoly( 8, 0 ) ; // z^4
            mPoly( 29, 0 ) = mPoly( 9, 0 ) * mPoly( 8, 0 ) ; // z^3 * x
            mPoly( 30, 0 ) = mPoly( 4, 0 ) * mPoly( 8, 0 ) ; // x^2 * z^2
            mPoly( 31, 0 ) = mPoly( 9, 0 ) * mPoly( 4, 0 ) ; // z * x^3
            mPoly( 32, 0 ) = mPoly(  4, 0 ) * mPoly( 7, 0  ) ; // x^2 * y * z
            mPoly( 33, 0 ) = mPoly( 6, 0 ) * mPoly( 9, 0 ) ; // x * y^2 * z
            mPoly( 34, 0 ) = mPoly( 5, 0 ) * mPoly( 8, 0 ) ; // x * y * z^2
        }

        inline const Cell< id_t > &
        Postprocessor::block_ids() const
        {
            return mBlockIDs ;
        }

        inline const Cell< index_t > &
        Postprocessor::my_element_indices() const
        {
            return mMyElementIndices ;
        }


    }
}
#endif //BELFEM_CL_FEM_POSTPROCESSOR_HPP