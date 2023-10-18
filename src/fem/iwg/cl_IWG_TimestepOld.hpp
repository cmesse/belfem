//
// Created by christian on 7/28/21.
//

#ifndef BELFEM_CL_IWG_TimestepOld_HPP
#define BELFEM_CL_IWG_TimestepOld_HPP

#include "cl_IWG.hpp"
#include "en_SolverEnums.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------
        class Element ;
        class Group ;
        /**
         * an abstract baseclass for timestepping problems
         */
        class IWG_TimestepOld: public IWG
        {
//------------------------------------------------------------------------------
            protected:
//------------------------------------------------------------------------------
            /**
             * this is the timestepping parameter.
             * Explicit:           theta = 0.0 ( not recommended )
             * Crank-Nicolson :    theta = 1/2 ( based on finite difference approach )
             * Galerkin:           theta = 2/3 ( based on finite element over time )
             * Backwards Implicit: theta = 1.0 ( recommended setting )
             */
            real mTheta     = 1.0;

            /**
             * this function computes the temperature for the current position
             * we use a function pointer to circumvent an unneccessary
             * if statement
             */
            real
            ( IWG_TimestepOld:: * mComputeT ) ( const uint aK );

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

            IWG_TimestepOld( const IwgType aType,
                          const IwgMode aMode=IwgMode::Iterative,
                          const SymmetryMode aSymmetryMode=SymmetryMode::PositiveDefiniteSymmetric,
                          const DofMode      aDofMode=DofMode::AllBlocksEqual,
                          const SideSetDofLinkMode aSideSetDofLinkMode=SideSetDofLinkMode::FacetOnly );

//------------------------------------------------------------------------------

            virtual ~IWG_TimestepOld() = default ;

//------------------------------------------------------------------------------


            /**
             * chooses the timestepping method. By default, we use
             * Euler Implicit, since it is very stable
             *
             * @param aEulerMethod
             */
            void
            set_euler_method( const EulerMethod aEulerMethod );

            void
            set_euler_method( const real aTheta );

//------------------------------------------------------------------------------

            /**
             * sets the timestep of the timestepping methd
             *
             * @param aDeltaTime timestep in seconds
             */
            void
            set_timestep( const real aDeltaTime );

//------------------------------------------------------------------------------

            /**
             * returns the theta-parameter for the timestepping
             */
            const real &
            theta () const;

//------------------------------------------------------------------------------

            /**
            * * these functions compute the temperature.
            * we use function pointers to avoid an if-statement
            *
            * @param aK     index of integration point
            * @return
            */

            real
            compute_T( const uint aK );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            real
            compute_Ttheta( const uint aK );

//------------------------------------------------------------------------------

            real
            compute_T0( const uint aK );

//------------------------------------------------------------------------------

            real
            compute_T1( const uint aK );

//------------------------------------------------------------------------------
        };
//---------------------------------------------------------------------------

        /**
         * returns the theta-parameter for the timestepping
         */
        inline const real &
        IWG_TimestepOld::theta () const
        {
            return mTheta ;
        }

//------------------------------------------------------------------------------

        inline real
        IWG_TimestepOld::compute_T( const uint aK )
        {
            return ( this->*mComputeT )( aK );
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_IWG_TimestepOld_HPP
