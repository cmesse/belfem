//
// Created by christian on 6/9/23.
//

#ifndef BELFEM_CL_FEM_CALCULATOR_HPP
#define BELFEM_CL_FEM_CALCULATOR_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"

#include "cl_Mesh.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "en_IWGs.hpp"

namespace belfem
{
    namespace fem
    {
        class Group ;
        class Element ;

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

        } // end namespace calculator

//------------------------------------------------------------------------------

        class Calculator
        {
            // link to group
            Group * mGroup ;

            // link to mesh
            Mesh  * mMesh ;

            // link to current element
            Element * mElement = nullptr ;

            const ModelDimensionality mDimensionality ;

            // switch telling if we are allocated
            bool mIsAllocated = false ;

            uint mNumberOfNodes = BELFEM_UINT_MAX ;
            uint mNumberOfCornerNodes = BELFEM_UINT_MAX ;
            uint mNumberOfIntegrationPoints = 0 ;

            uint mMasterIndex = BELFEM_UINT_MAX ;

                  IntegrationData * mDomainIntegration    = nullptr ;
                  IntegrationData * mLinearIntegration    = nullptr ;
                  IntegrationData * mThinShellIntegration = nullptr ;
                  IntegrationData * mMasterIntegration    = nullptr ;
                  IntegrationData * mSlaveIntegration     = nullptr ;

                  // flag telling if element is straight or curved
                  bool mIsCurved = false ;

            //! stiffness matrix
            Matrix< real > mK ;

            //! mass matrix
            Matrix< real > mM ;

            //! load vector
            Vector< real > mf ;

            //! dof vector for last timestep
            Vector< real > mq0 ;

            //! current dof vector for last timestep
            Vector< real > mq ;


            //! normal vector
            Vector< real > mNormal ;
            uint mNormalIndex      = BELFEM_UINT_MAX ;
            real mSurfaceIncrement = BELFEM_QUIET_NAN ;

            //! dof labels for next timestep
            Cell< string >     mDofLabels ;

            //! dof labels for old timestep
            Cell< string >     mOldDofLabels ;

            //! custom matrices
            Cell< calculator::VectorData * > mVectors ;
            Cell< calculator::MatrixData * > mMatrices ;

            //! maps
            Map< string , calculator::VectorData * > mVectorMap ;
            Map< string , calculator::MatrixData * > mMatrixMap ;


            // pointers for faster access
            Matrix< real > mX  ;  // node coordinates
            Matrix< real > mXc ;  // node coordinates at corners
            calculator::MatrixData * mJ    = nullptr ;  // jacobian matrix
            calculator::MatrixData * mInvJ = nullptr ;  // jacobian matrix
            calculator::MatrixData * mN    = nullptr ;  // node interpolatoin operator
            calculator::MatrixData * mdN   = nullptr ;  // derivative function
            calculator::MatrixData * mB    = nullptr ;  // gradient operator
            calculator::MatrixData * mC    = nullptr ;  // curl operator
            calculator::MatrixData * mE    = nullptr ;  // edge/face operator

            // for faces
            Matrix< real > mXm ;
            calculator::MatrixData * mJm = nullptr ;
            calculator::MatrixData * mBm = nullptr ;
            calculator::MatrixData * mCm = nullptr ;
            calculator::MatrixData * mEm = nullptr ;

            Matrix< real > mXs ;
            calculator::MatrixData * mJs = nullptr ;
            calculator::MatrixData * mBs = nullptr ;
            calculator::MatrixData * mCs = nullptr ;
            calculator::MatrixData * mEs = nullptr ;

            real mDetJ      = BELFEM_QUIET_NAN ;
            uint mDetJIndex = BELFEM_UINT_MAX ;
            real mRadius    = BELFEM_QUIET_NAN ; // only needed if axisymmetric

            // function to link the slave integration data
            IntegrationData * ( Calculator::*mFunSlaveIntegration )( const Element * aElement ) ;

            real ( * mFunInvertJ )( const Matrix< real > & aA, Matrix< real > & aB );

            const Vector< real > & ( Calculator::*mFunNormal )( const uint aIndex );

            // function for node interpolator
            const Matrix< real > & (  Calculator::*mFunN )( const uint aIndex );

            // function for gradient operator
            const Matrix< real > & (  Calculator::*mFunB )( const uint aIndex );

            // volume increment
            real ( Calculator::*mFundV )( const uint aIndex );

            // surface increment
            real ( Calculator::*mFundS )( const uint aIndex );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Calculator( Group * aGroup, const ModelDimensionality aDimensionality );

//------------------------------------------------------------------------------

            ~Calculator() ;

//------------------------------------------------------------------------------

            /**
             * called by dof manager
             */
            void
            allocate();

//------------------------------------------------------------------------------

            void
            link( Element * aElement );

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
             * return the load vector
             */
            Vector< real > &
            f() ;

//------------------------------------------------------------------------------

            /**
             * return the dof vector at current timestep
             */
            Vector< real > &
            q() ;

//------------------------------------------------------------------------------

            /**
             * return the dof vector at the last timestep
             */
            Vector< real > &
            q0() ;

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
             const Vector< real > &
             node_data( const string & aNodeField );

//------------------------------------------------------------------------------

            /**
             * return a vector object
             */
            const Vector< real > &
            nedelec_data_linear( const string & aEdgeField );

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
            Vector< real > &
            vector( const string & aLabel );

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
             * node interpolator, but as vector
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
            initialize_integration( const ElementType aElementType,
                                         const InterpolationType aInterpolationType );

//------------------------------------------------------------------------------

            void
            set_integration_order( const uint aOrder );

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

            uint
            num_intpoints() const ;

//------------------------------------------------------------------------------
        private:
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

            IntegrationData *
            slave_integration_2d( const Element * aElement ) ;

            IntegrationData *
            slave_integration_tet( const Element * aElement ) ;

            IntegrationData *
            slave_integration_hex( const Element * aElement ) ;

//------------------------------------------------------------------------------

            // node interpolator for scalar fields
            const Matrix< real > &
            Nscalar( const uint aIndex ) ;

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
            normal_hex( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_2D_3D( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_axsymmx( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dV_axsymmy( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_line( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_axsymmx( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_axsymmy( const uint aIndex );

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

        inline Vector< real > &
        Calculator::vector( const string & aLabel )
        {
            return mVectorMap( aLabel )->vector() ;
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
            if ( aIndex != mInvJ->index() )
            {
                // remember index
                mInvJ->set_index( aIndex );
                mDetJIndex = aIndex ;

                // compute inverse and remember determinant
                mDetJ = ( * mFunInvertJ )( this->J( aIndex ), mInvJ->matrix() );
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

        inline real
        Calculator::node_interp( const uint aIndex, const Vector< real > & aNodeValues ) const
        {
            return dot( mDomainIntegration->phi( aIndex ),  aNodeValues );
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
            return ( this->*mFundS )( aIndex );
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV( const uint aIndex )
        {
            return ( this->*mFundV)( aIndex );
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dV_2D_3D( const uint aIndex )
        {
            if( mDetJIndex != aIndex )
            {
                mDetJIndex = aIndex ;
                mDetJ = det( this->J( aIndex ) );
            }
            return mDetJ ;
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

                mDetJ *= det( this->J( aIndex ) );
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

                mDetJ *= det( this->J( aIndex ) );
            }
            return mDetJ ;
        }

//------------------------------------------------------------------------------


        inline real
        Calculator::dS_line( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
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
                        mXm.col( 0 ) )
                    * mSurfaceIncrement ;

        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dS_axsymmy( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            return dot( mMasterIntegration->phi( aIndex ).vector_data(),
                        mXm.col( 1 ) )
                   * mSurfaceIncrement ;

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

        inline Vector< real > &
        Calculator::f()
        {
            return mf ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Calculator::q()
        {
            return mq ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Calculator::q0()
        {
            return mq0 ;
        }

//------------------------------------------------------------------------------

        inline uint
        Calculator::num_intpoints() const
        {
            return mNumberOfIntegrationPoints ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_CALCULATOR_HPP