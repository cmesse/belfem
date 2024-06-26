//
// Created by christian on 12/22/21.
//

#ifndef BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP
#define BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP

#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_HPhi_Tri6 : public IWG_Maxwell
        {
            // work vector for thin shell stiffnesses
            real mWork[ 16 ];

            // tangential components of H per layer
            Vector< real > mHt ;

            // curl operator for thin shell
            Matrix< real > mC ;
            Vector< real > mdCdx ;
            Vector< real > mdCdy ;

            // link to function
            void
            ( IWG_Maxwell_HPhi_Tri6::*mFunJacobian )
            (       Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );


            // matrix for element values that need to be integrated
            // over thin shells
            Matrix< real > mLayerData ;

            // work vector for normal component of H
            Vector< real > mH ;

            // edge interpolation vector
            Vector< real > mE ;

            // edge index vector
            Vector< uint > mU ;

            Vector< real > mHn ;
            Vector< real > mW ;

            Matrix< real > mAHn;
            Vector< real > mBHn;

            //Matrix< real > mArho ;
            //Vector< real > mBrho ;
            Matrix< real > mNxi ;
            Vector< real > mCrho ;
            Matrix< real > mRho ;
            Matrix< real > mdRhodx ;
            Matrix< real > mdRhody ;
            Matrix< real > mJz ;
            Matrix< real > mJc ;
            Vector< int >  mPivot ;
            Matrix< real > mL ; // obsolete
            Matrix< real > mLm ;
            Matrix< real > mLs ;

            Vector< real > mN ;
            Vector< real > mdNdeta ;
            Vector< real > md2Ndeta2 ;
            Vector< real > mdEdxi ;
            Vector< real > mdEdy ;
            Vector< real > md2Edy2 ;
            Vector< real > md2Edxdy ;


            Vector< real >  mXi ;

            Vector< real >  mXi0 ;
            Vector< real >  mXi1 ;

            Vector< real >  mEta ;
            Vector< real >  mPhi ;
            Vector< real >  mPsi ;
            Vector< real >  mSigma ;
            Vector< real >  mX ;
            Vector< real >  mY ;

            // transformation matrix
            Matrix< real > mTm ;
            Matrix< real > mTs ;

            Matrix< real > mK ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_HPhi_Tri6();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_HPhi_Tri6();

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            real
            compute_element_current( Element * aElement ) ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            link_jacobian_function( Group * aGroup );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_superconductor(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_scair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_scfm(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_fmfm(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_fmair(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_wave(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmetry_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmetry_air(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_symmetry_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_antisymmetry_ferro(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_cut(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs_thinshell(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS ) ;

//------------------------------------------------------------------------------

            void
            compute_layer_mass(
                    const uint aLayer,
                    const real aE0,
                    const real aE1,
                    Matrix< real > & aM );

//------------------------------------------------------------------------------

            void
            compute_layer_stiffness(
                    const uint aLayer,
                    const uint aIntPoint,
                    const real aE0,
                    const real aE1,
                    const Vector< real > & aHt,
                    const real aHn,
                    const real adHndx,
                    const real aXLength,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            void
            compute_thinshell_interface(
                    Element * aElement,
                    const bool aMaster,
                    const uint aLayer,
                    const uint aIntPoint,
                    const real aE0,
                    const real aE1,
                    Matrix< real > & aL );

//------------------------------------------------------------------------------

            void
            compute_layer_stabilizer(
                    Element * aElement,
                    const Vector< real > & aHt,
                    const real aXLength,
                    const real aYLength,
                    const real aNminus1,
                    Matrix< real > & aM ,
                    Matrix< real > & aK );

//------------------------------------------------------------------------------

            const Vector< real > &
            collect_q0_thinshell( Element * aElement );

//------------------------------------------------------------------------------

            const Vector< real > &
            eval_quad9( const real aX, const real aY );

//------------------------------------------------------------------------------

            const Vector< real > &
            eval_quad9_dx( const real aX, const real aY );

//------------------------------------------------------------------------------

            const Vector< real > &
            eval_quad9_dy( const real aX, const real aY );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            ( this->*mFunJacobian ) ( aElement, aJacobian, aRHS ) ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_superconductor(
                aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_air_phi_higher_order(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_ferro(
            Element        * aElement,
            Matrix< real > & aJacobian,
            Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_ferro(
                    aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_scfm(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // make sure that the facet is correctly oriented
            BELFEM_ASSERT( static_cast< DomainType >( aElement->master()->element()->physical_tag() ) == DomainType::Conductor,
                           "Master element of facet %lu must be of DomainType::Conductor", ( long unsigned int ) aElement->id() );

            BELFEM_ASSERT( static_cast< DomainType >( aElement->slave()->element()->physical_tag() ) == DomainType::Ferro,
                           "Master element of facet %lu must be of DomainType::Ferro", ( long unsigned int ) aElement->id() );


#ifdef BELFEM_FERRO_LINEAR
            this->compute_interface_ha_tri6_tri3( aElement, aJacobian, mGroup->work_K() );
#else
            this->compute_interface_ha_tri6( aElement, aJacobian, mGroup->work_K() );
#endif

            // compute RHS
            aRHS = aJacobian * this->collect_q0_ha_2d( aElement );

            // compone matrices
            aJacobian += mGroup->work_K() * this->timestep() ;
        }
//------------------------------------------------------------------------------


        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmfm(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

#ifdef BELFEM_FERRO_LINEAR
            this->compute_interface_aa_tri3( aElement, aJacobian, aRHS );
#else
            this->compute_interface_aa_tri6( aElement, aJacobian, aRHS );
#endif

        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            IWG_Maxwell::compute_jacobian_and_rhs_cut( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell_HPhi_Tri6::eval_quad9( const real aX, const real aY )
        {
            mCrho( 0 ) = 1.0 ;
            mCrho( 1 ) = aX ;
            mCrho( 2 ) = aY ;
            mCrho( 3 ) = aX*aX ;
            mCrho( 4 ) = aX*aY ;
            mCrho( 5 ) = aY*aY ;
            mCrho( 6 ) = mCrho( 3 )*aY ;
            mCrho( 7 ) = aX*mCrho( 5 ) ;
            mCrho( 8 ) = mCrho( 3 )*mCrho( 5 );
            return mCrho ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell_HPhi_Tri6::eval_quad9_dx( const real aX, const real aY )
        {
            mCrho( 0 ) = 0.0 ;
            mCrho( 1 ) = 1.0 ;
            mCrho( 2 ) = 0.0 ;
            mCrho( 3 ) = aX + aX ;
            mCrho( 4 ) = aY ;
            mCrho( 5 ) = 0.0 ;
            mCrho( 6 ) = mCrho( 3 )*aY ;
            mCrho( 7 ) = aY*aY ;
            mCrho( 8 ) = mCrho( 3 )*mCrho( 7 );
            return mCrho ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        IWG_Maxwell_HPhi_Tri6::eval_quad9_dy( const real aX, const real aY )
        {
            mCrho( 0 ) = 0.0 ;
            mCrho( 1 ) = 0.0 ;
            mCrho( 2 ) = 1.0 ;
            mCrho( 3 ) = 0.0 ;
            mCrho( 4 ) = aX ;
            mCrho( 5 ) = aY + aY ;
            mCrho( 6 ) = aX * aX ;
            mCrho( 7 ) = mCrho( 5 )*aY ;
            mCrho( 8 ) = mCrho( 5 )*mCrho( 6 );
            return mCrho ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWG_MAXWELL_HPHI_TRI6_HPP
