//
// Created by christian on 6/9/23.
//

#ifndef BELFEM_CL_FEM_SWAP_HPP
#define BELFEM_CL_FEM_SWAP_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace fem
    {
        class Group ;
        class Swap
        {
            Group * mGroup ;

            bool mHaveThermal     = false ;
            bool mHaveMechanical  = false ;
            bool mHaveMaxwellPhi = false ;
            bool mHaveMaxwellA   = false ;
            bool mHaveMaxwellH   = false ;
            bool mHaveMaxwellJ   = false ;
            bool mHaveLagrange   = false ;

            // link to mesh
            Mesh  * mMesh ;

            //! Jacobian
            Matrix< real > mJ ;

            //! inverse of Jacobian
            Matrix< real > mInvJ ;

            //! coordinates for the element
            Matrix< real > mX ;

            //! coordinates for the master element
            Matrix< real > mXm ;

            //! coordinates for the slave element
            Matrix< real > mXs ;

            //! help matrix for gradient operator
            Matrix< real > mdN ;

            //! gradient operator
            Matrix< real > mB ;

            //! gradient operator of master element
            Matrix< real > mBm ;

            //! gradient operator of slave element
            Matrix< real > mBs ;

            //! matrix for curl operator
            Matrix< real > mC ;

            //! vector for generalized degrees of freedom
            Vector< real > mq ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // TEMPERATURE FIELDS
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            //! vector for nodal temperatures
            Vector< real > mT ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // DISPLACEMENT FIELDS
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            //! nodal displacement in x-direction
            Vector< real > mux ;

            //! nodal displacement in y-direction
            Vector< real > muy ;

            //! nodal displacement in x-direction
            Vector< real > muz ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Maxwell Fields
            // - - - - - - - - - - - - - - - - - - - - - - - - - - -

            //! scalar potential on nodes
            Vector< real > mphi ;

            //! scalar potential on master nodes
            Vector< real > mphim ;

            //! scalar potential on slave nodes
            Vector< real > mphis ;

            // -   -   -   -   -   -   -   -   -   -   -   -   -   -

            //! vector potential on nodes in x-direction
            Vector< real > max ;

            //! vector potential on nodes in y-direction
            Vector< real > may ;

            //! vector potential on nodes in z-direction
            Vector< real > maz ;

            //! vector potential
            Vector< real > ma ;

            //! vector potential on master
            Vector< real > mam ;

            //! vector potential on slave
            Vector< real > mas ;

            // -   -   -   -   -   -   -   -   -   -   -   -   -   -


            //! magnetic field on edges
            Vector< real > mh ;

            //! magnetic field on master
            Vector< real > mhm ;

            //! magnetic field on slave
            Vector< real > mhs ;

            // -   -   -   -   -   -   -   -   -   -   -   -   -   -

            //! current density in x-direction
            Vector< real > mjx ;

            //! current density in y-direction
            Vector< real > mjy ;

            //! current density in z-direction
            Vector< real > mjz ;

            // -   -   -   -   -   -   -   -   -   -   -   -   -   -

            //! Lagrange multiplier
            Vector< real > mlambda ;

            //! Lagrange multiplier master
            Vector< real > mlambdam ;

            //! Lagrange multiplier slave
            Vector< real > mlambdas ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Swap( Group * aGroup );

//------------------------------------------------------------------------------

            ~Swap() = default ;

//------------------------------------------------------------------------------

            /***
             * the Jacobian Matrix for the Geometry
             * @return
             */
            Matrix< real > &
            J();

//------------------------------------------------------------------------------

            /***
             * the inverted Jacobian for the Geometry
             * @return
             */
            Matrix< real > &
            invJ();

//------------------------------------------------------------------------------

            /**
             * Node coordinates for the elment
             * @return
             */
            Matrix< real > &
            X();

//------------------------------------------------------------------------------

            /**
             * Node coordinates for the master element
             * ( facet only )
             */
            Matrix< real > &
            Xm();

//------------------------------------------------------------------------------

            /**
            * Node coordinates of the slave element
            * ( facet only )
            */
            Matrix< real > &
            Xs();

//------------------------------------------------------------------------------

            /**
             * matrix for gradient operator
             */
             Matrix< real > &
             B();

//------------------------------------------------------------------------------

            /**
             * matrix for gradient operator of master
             * ( facet only )
             */
            Matrix< real > &
            Bm();

//------------------------------------------------------------------------------

            /**
             * matrix for gradient operator of slave
             * ( facet only )
             */
            Matrix< real > &
            Bs();

//------------------------------------------------------------------------------

            /**
             * matrix for curl operator,
             * alternatively used for material property tensor
             */
            Matrix< real > &
            C();

//------------------------------------------------------------------------------
//      DOFS
//------------------------------------------------------------------------------
            /**
             * vector for generalized degree of freedom
             */
             Vector< real> &
             q();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      THERMAL
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            /**
             * vector for nodal temperatures
             */
            Vector< real> &
            T();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      MECHANICAL
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * vector for displacement in x-direction
             */
            Vector< real> &
            ux();

            /**
             * vector for displacement in y-direction
             */
            Vector< real> &
            uy();

            /**
             * vector for displacement in z-direction
             */
            Vector< real> &
            uz();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      MAXWELL PHI
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * vector for magnetic scalar field
             */
            Vector< real> &
            phi();

            /**
             * vector for magnetic scalar field on master
             */
            Vector< real> &
            phim();

            /**
            * vector for magnetic scalar field on slave
            */
            Vector< real> &
            phis();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      MAXWELL A
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * vector for magnetic scalar field
             */
            Vector< real> &
            a();

            /**
             * vector for magnetic scalar field on master
             */
            Vector< real> &
            am();

            /**
            * vector for magnetic scalar field on slave
            */
            Vector< real> &
            as();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      MAXWELL H
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * magnetic field
             */
            Vector< real> &
            h();

            /**
             * magnetic field on master
             */
            Vector< real> &
            hm();

            /**
            * magnetic field on slave
            */
            Vector< real> &
            hs();

// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      MAXWELL multipliers
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * Lagrange multiplier
             */
            Vector< real> &
            lambda();

            /**
             * Lagrange multiplier from master
             */
            Vector< real> &
            lambdam();

            /**
            * Lagrange multiplier from slave
            */
            Vector< real> &
            lambdas();

//------------------------------------------------------------------------------
        /**
         * this subroutine allocates the memory of the matrices
         */
        void
        allocate();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * detects which fields we are using
             */
            void
            set_switches();

//------------------------------------------------------------------------------
        };


        inline Matrix< real > &
        Swap::J()
        {
            return mJ ;
        }

        inline Matrix< real > &
        Swap::invJ()
        {
            return mInvJ ;
        }

        inline Matrix< real > &
        Swap::X()
        {
            return mX ;
        }

        inline Matrix< real > &
        Swap::Xm()
        {
            return mXm ;
        }


        inline Matrix< real > &
        Swap::Xs()
        {
            return mXs ;
        }


        inline Matrix< real > &
        Swap::Swap::B()
        {
            return mB ;
        }


        inline Matrix< real > &
        Swap::Bm()
        {
            return mB ;
        }


        inline Matrix< real > &
        Swap::Bs()
        {
            return mBs ;
        }


        inline Matrix< real > &
        Swap::C()
        {
            return mC ;
        }


        inline Vector< real> &
        Swap::q()
        {
            return mq ;
        }

        inline Vector< real> &
        Swap::T()
        {
            return mT ;
        }

        inline Vector< real> &
        Swap::ux()
        {
            return mux ;
        }

        inline Vector< real> &
        Swap::uy()
        {
            return muy ;
        }

        inline Vector< real> &
        Swap::uz()
        {
            return muz ;
        }

        inline Vector< real> &
        Swap::phi()
        {
            return mphi ;
        }

        inline Vector< real> &
        Swap::phim()
        {
            return mphim ;
        }

        inline Vector< real> &
        Swap::phis()
        {
            return mphis ;
        }

        inline Vector< real> &
        Swap::a()
        {
            return ma ;
        }

        inline Vector< real> &
        Swap::am()
        {
            return mam ;
        }

        inline Vector< real> &
        Swap::as()
        {
            return mas ;
        }

        inline Vector< real> &
        Swap::h()
        {
            return mh ;
        }

        inline Vector< real> &
        Swap::hm()
        {
            return mhm ;
        }

        inline Vector< real> &
        Swap::hs()
        {
            return mhs ;
        }

        inline Vector< real> &
        Swap::lambda()
        {
            return mlambda ;
        }

        inline Vector< real> &
        Swap::lambdam()
        {
            return mlambdam ;
        }

        inline Vector< real> &
        Swap::lambdas()
        {
            return mlambdas ;
        }
    }
}

#endif //BELFEM_CL_FEM_SWAP_HPP
