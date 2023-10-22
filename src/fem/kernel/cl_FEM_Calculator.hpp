//
// Created by christian on 6/9/23.
//

#ifndef BELFEM_CL_FEM_CALCULATOR_HPP
#define BELFEM_CL_FEM_CALCULATOR_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "fn_dot.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"

#include "cl_Mesh.hpp"
#include "cl_IF_IntegrationData.hpp"


namespace belfem
{
    //
    enum class Nfunction
    {
        Scalar = 0,
        Vec2d  = 1,
        Vec3d  = 2,
        UNDEFINED = 3
    };

    enum class Bfunction
    {
        Gradient     = 0,
        Planestress  = 1,
        Voigt        = 2,
        UNDEFINED = 3
    };

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

            Vector< real > mIntegrationWeights ;

            // switch telling if we are allocated
            bool mIsAllocated = false ;

            uint mNumberOfNodes = BELFEM_UINT_MAX ;

            uint mMasterIndex = BELFEM_UINT_MAX ;

                  IntegrationData *   mDomainIntegration    = nullptr ;
                  IntegrationData *   mThinShellIntegration = nullptr ;
            Cell< IntegrationData * > mMasterIntegration ;
            Cell< IntegrationData * > mSlaveIntegration ;

            //! stiffness matrix
            Matrix< real > mK ;

            //! load vector
            Vector< real > mf ;

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
            Matrix< real > mX ;  // node coordinates
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

            // function to link the slave integration data
            const IntegrationData * ( Calculator::*mFunSlaveIntegration )( const Element * aElement ) ;

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

            Calculator( Group * aGroup );

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
             * return the integation weights
             */
            const Vector< real > &
            weights() const ;

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
            initialize_master_integration( const ElementType aElementType,
                                           const InterpolationType aInterpolationType );

//------------------------------------------------------------------------------

            void
            initialize_slave_integration( const ElementType aElementType,
                                           const InterpolationType aInterpolationType );

//------------------------------------------------------------------------------

            void
            allocate_memory();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------



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

            const IntegrationData *
            slave_integration_2d( const Element * aElement ) ;

            const IntegrationData *
            slave_integration_tet( const Element * aElement ) ;

            const IntegrationData *
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
            dS_line( const uint aIndex );

//------------------------------------------------------------------------------

            real
            dS_triangle( const uint aIndex );

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

        inline const Vector< real > &
        Calculator::weights() const
        {
            return mIntegrationWeights ;
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
                mJ->matrix().matrix_data() =
                          mIntegration->dNdXi( aIndex ) * mX ;
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
            return mIntegration->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Calculator::N2D( const uint aIndex )
        {
            if( mN->index() != aIndex )
            {
                // precomputed data
                const Vector< real > & tPhi = mIntegration->phi( aIndex );

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
                const Vector< real > & tPhi = mIntegration->phi( aIndex );

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
                        this->invJ( aIndex ) * mIntegration->dNdXi( aIndex );
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
                tdN = this->invJ( aIndex ) * mIntegration->dNdXi( aIndex );

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
                tdN = this->invJ( aIndex ) * mIntegration->dNdXi( aIndex );

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
            return dot( mIntegration->phi( aIndex ),  aNodeValues );
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
        Calculator::dS_line( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            return mSurfaceIncrement ;
        }

//------------------------------------------------------------------------------

        inline real
        Calculator::dS_triangle( const uint aIndex )
        {
            // compute the normal if it hasn't been computed so far
            ( this->*mFunNormal ) ( aIndex );

            // return the surface increment
            return mSurfaceIncrement ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_CALCULATOR_HPP
