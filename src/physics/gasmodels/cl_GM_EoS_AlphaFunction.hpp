//
// Created by Christian Messe on 14.09.19.
//

#ifndef BELFEM_CL_GM_EOS_ALPHAFUNCTION_HPP
#define BELFEM_CL_GM_EOS_ALPHAFUNCTION_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        /**
         * The Alpha function is part of the Cubic gas model
         */
        class AlphaFunction
        {
//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            //! critical temperature
            const real             mTcrit;

            //! inverse of critical temperature
            const real             mInvTcrit;

            //! last used T
            real                   mT[ 3 ] = { 0.0 };

            // constants to be used by the function
            const real             mC1;
            const real             mC2;
            const real             mC3;

            //! 0 : alpha
            //! 1: dalphadT
            //! 2 : d2alphadT
            real                  mAlpha[ 3 ];

            // help values
            real                  mWork[ 4 ];

            // temperature with minimal alpha
            real                  mTmin;
            real                  mAlphaMin;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            //! empty constructor
            AlphaFunction();

//----------------------------------------------------------------------------

            /**
             * a special constructor that initializes coefficients
             */
            AlphaFunction( const real & aTcrit, const Vector< real > & aCoeffs );

//----------------------------------------------------------------------------

            /**
             * a special constructor that initializes coefficients
             */
            AlphaFunction(
                    const real & aTcrit,
                    const real & aC1,
                    const real & aC2,
                    const real & aC3 );

//----------------------------------------------------------------------------

            virtual ~AlphaFunction() = default;

//----------------------------------------------------------------------------


            virtual real
            alpha( const real & aT );

//----------------------------------------------------------------------------

            virtual real
            dalphadT( const real & aT );

//----------------------------------------------------------------------------

            virtual real
            d2alphadT2( const real & aT );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            virtual void
            eval( const real & aT, const int aDeriv );


//----------------------------------------------------------------------------

            void
            find_minimum();

//----------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------

        /**
         * an empty function that does nothing
         */
        class AlphaFunction_Empty : public AlphaFunction
        {
//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_Empty();

//----------------------------------------------------------------------------
            real
            alpha( const real & aT );

//----------------------------------------------------------------------------

            real
            dalphadT( const real & aT );

//----------------------------------------------------------------------------

            real
            d2alphadT2( const real & aT );
        };

//----------------------------------------------------------------------------

        /*
         * The classic alpha function
         *
         */
        class AlphaFunction_Classic : public AlphaFunction
        {

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_Classic(
                    const real & aTcrit,
                    const real & aC1,
                    const real & aC2 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real & aT, const int aDeriv );

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
        /*
        * Coquelet, Chapoy, Richon for SRK
        * based on Mathias and Coepman, but cut off where dalpha/dT = 0
        *
        * 10.1023/B:IJOT.0000022331.46865.2f
        * 10.1016/0378-3812(83)80084-3
        */
        class AlphaFunction_MC : public AlphaFunction
        {
            real mX[4];
            real mF[3];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_MC(
                    const real & aTcrit,
                    const real & aC1,
                    const real & aC2,
                    const real & aC3 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real & aT, const int aDeriv );

        };

//----------------------------------------------------------------------------

        /*
         * Coquelet, Chapoy, Richon for Peng Robinson
         *
         * 10.1023/B:IJOT.0000022331.46865.2f
         *
         */
        class AlphaFunction_CCR : public AlphaFunction
        {
            // c1*( 1-T/Tcrit ) and derivatives
            real mF[2];

            // 1 + c2*(1-sqrt(T/Tcrit))^2 + c3*(1-sqrt(T/Tcrit))^3
            real mG[3];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_CCR(
                    const real & aTcrit,
                    const real & aC1,
                    const real & aC2,
                    const real & aC3 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real & aT, const int aDeriv );

        };
//----------------------------------------------------------------------------

        /**
         * Mahmoodi and Sedigh
         *
         * 10.1016/j.fluid.2016.12.015
         *
         */
        class AlphaFunction_PM : public AlphaFunction
        {
            // 0 : f     ( help function )
            // 1 : dfdT
            // 2 : d2fdT
            real                  mF[ 3 ];

            // sqrt( T / T_crit )
            real                  mX;

            // 0 : Y     ( help function )
            // 1 : dYdT
            // 2 : d2YdTh
            real                  mY[ 3 ];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_PM(
                    const real           & aTcrit,
                    const Vector< real > & aCoeffs );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real & aT, const int aDeriv );

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */

#endif //BELFEM_CL_GM_EOS_ALPHAFUNCTION_HPP
