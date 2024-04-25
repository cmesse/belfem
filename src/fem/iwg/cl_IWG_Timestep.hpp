//
// Created by christian on 7/18/23.
//

#ifndef BELFEM_CL_IWG_TIMESTEP_HPP
#define BELFEM_CL_IWG_TIMESTEP_HPP

#include "cl_IWG.hpp"
#include "en_SolverEnums.hpp"


namespace belfem
{
    namespace fem
    {
        class Element ;

        class IWG_Timestep: public IWG
        {

            EulerMethod mMethod;

            // timestepping pointer
            void
            ( IWG_Timestep:: * mTimestep )(
                    Matrix< real > & aM,
                    Matrix< real > & aK,
                    Vector< real > & aF,
                    Vector< real > & aQ0  );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            real mDeltaTime = 1.0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Timestep( const IwgType aType,
                          const ModelDimensionality aDimensionality,
                          const IwgMode aMode=IwgMode::Iterative,
                          const SymmetryMode aSymmetryMode=SymmetryMode::PositiveDefiniteSymmetric,
                          const DofMode      aDofMode=DofMode::AllBlocksEqual,
                          const SideSetDofLinkMode aSideSetDofLinkMode=SideSetDofLinkMode::FacetOnly );


            virtual ~IWG_Timestep() = default ;

            // set the timestepping method
            void
            set_timestepping_method( const EulerMethod aMethod, const bool aHaveStiffness = true );

            // return the timestepping method
            EulerMethod
            method() const;

//------------------------------------------------------------------------------

            virtual void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
             * a default interface to compute the matrices of a timest epping
             * scheme of shape of
             * \f$ M \, \dot q + K ,\, q = f \f$
             *
             * @param aElement
             * @param aM    : mass matrix
             * @param aK    : stiffness matrix
             * @param aF    : load vector
             */
            virtual void
            compute_mkf( Element * aElement,
                         Matrix< real > & aM,
                         Matrix< real > & aK,
                         Vector< real > & aF,
                         const Vector< real > & aQ );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            static_subfield( Matrix< real > & aM, Matrix< real > & aK, Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------
            void
            explicit_euler( Matrix< real > & aM, Matrix< real > & aK, Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            crank_nicolson( Matrix< real > & aM, Matrix< real > & aK, Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            galerkin( Matrix< real > & aM, Matrix< real > & aK, Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            backwards_euler( Matrix< real > & aM,
                             Matrix< real > & aK,
                             Vector< real > & aF,
                             Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            derivative( Matrix< real > & aM, Matrix< real > & aK,
                        Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            euler_no_stiffness( Matrix< real > & aM, Matrix< real > & aK,
                                Vector< real > & aF, Vector< real > & aQ0  );

//------------------------------------------------------------------------------

            void
            derivative_no_stiffness( Matrix< real > & aM, Matrix< real > & aK,
                                     Vector< real > & aF, Vector< real > & aQ0  );


//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

                                                inline void
        IWG_Timestep::explicit_euler( Matrix< real > & aM,
                                      Matrix< real > & aK,
                                      Vector< real > & aF,
                                      Vector< real > & aQ0  )
        {
            aF.vector_data() *= mDeltaTime ;
            aK.matrix_data() *= mDeltaTime ;
            aF.vector_data() += ( aM.matrix_data() - aK.matrix_data() ) * aQ0.vector_data() ;
        }

        inline void
        IWG_Timestep::crank_nicolson( Matrix< real > & aM,
                                      Matrix< real > & aK,
                                      Vector< real > & aF,
                                      Vector< real > & aQ0  )
        {
            aF.vector_data() *=       mDeltaTime ;
            aK.matrix_data() *= 0.5 * mDeltaTime ;
            aF.vector_data() += ( aM.matrix_data() - aK.matrix_data() ) * aQ0.vector_data() ;
            aM.matrix_data() += aK.matrix_data() ;
        }

        inline void
        IWG_Timestep::galerkin( Matrix< real > & aM,
                                Matrix< real > & aK,
                                Vector< real > & aF,
                                Vector< real > & aQ0 )
        {
            aF.vector_data() *= mDeltaTime ;
            aK.matrix_data() *= mDeltaTime / 3. ;
            aF.vector_data() += ( aM.matrix_data() - aK.matrix_data() ) * aQ0.vector_data() ;
            aM.matrix_data() += aK.matrix_data() ;
            aM.matrix_data() += aK.matrix_data() ;
        }

        inline void
        IWG_Timestep::backwards_euler( Matrix< real > & aM,
                                       Matrix< real > & aK,
                                       Vector< real > & aF,
                                       Vector< real > & aQ0  )
        {
            aF.vector_data() *= mDeltaTime ;
            aK.matrix_data() *= mDeltaTime ;
            aF.vector_data() += aM.matrix_data() * aQ0.vector_data() ;
            aM.matrix_data() += aK.matrix_data() ;
        }


        inline void
        IWG_Timestep::derivative( Matrix< real > & aM,
                                  Matrix< real > & aK,
                                  Vector< real > & aF,
                                  Vector< real > & aQ0  )
        {
            aF -= aK * aQ0 ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Timestep::euler_no_stiffness( Matrix< real > & aM,
                                      Matrix< real > & aK,
                                      Vector< real > & aF,
                                      Vector< real > & aQ0  )
        {
            aF.vector_data() *= mDeltaTime ;
            aF.vector_data() += aM.matrix_data() * aQ0.vector_data() ;
        }

        inline void
        IWG_Timestep::derivative_no_stiffness( Matrix< real > & aM,
                                               Matrix< real > & aK,
                                               Vector< real > & aF,
                                               Vector< real > & aQ0  )
        {
            // do nothing
        }

        inline void
        IWG_Timestep::static_subfield( Matrix< real > & aM,
                                       Matrix< real > & aK,
                                       Vector< real > & aF,
                                       Vector< real > & aQ0  )
        {
            aM.matrix_data() = aK.matrix_data() * mDeltaTime ;
            aF.vector_data() *= mDeltaTime ;
        }

//------------------------------------------------------------------------------

        // return the timestepping method
        inline EulerMethod
        IWG_Timestep::method() const
        {
            return mMethod ;
        }


//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IWG_TIMESTEP_HPP
