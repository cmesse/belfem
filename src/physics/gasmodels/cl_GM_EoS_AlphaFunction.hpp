/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

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
            mutable real           mT[ 3 ] = { 0.0 };

            // constants to be used by the function
            const real             mC1;
            const real             mC2;
            const real             mC3;

            //! 0 : alpha
            //! 1: dalphadT
            //! 2 : d2alphadT
            mutable real          mAlpha[ 3 ];

            // help values
            mutable real          mWork[ 4 ];

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
            AlphaFunction( const real aTcrit, const Vector< real > & aCoeffs );

//----------------------------------------------------------------------------

            /**
             * a special constructor that initializes coefficients
             */
            AlphaFunction(
                    const real aTcrit,
                    const real & c1,
                    const real & c2,
                    const real & c3 );

//----------------------------------------------------------------------------

            virtual ~AlphaFunction() = default;

//----------------------------------------------------------------------------


            virtual real
            alpha( const real T ) const;

//----------------------------------------------------------------------------

            virtual real
            dalphadT( const real T ) const;

//----------------------------------------------------------------------------

            virtual real
            d2alphadT2( const real T ) const;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            virtual void
            eval( const real T, const int aDeriv ) const;


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
            alpha( const real T ) const;

//----------------------------------------------------------------------------

            real
            dalphadT( const real T ) const;

//----------------------------------------------------------------------------

            real
            d2alphadT2( const real T ) const;
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
                    const real aTcrit,
                    const real & c1,
                    const real & c2 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real T, const int aDeriv ) const;

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
            mutable real mX[4];
            mutable real mF[3];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_MC(
                    const real aTcrit,
                    const real & c1,
                    const real & c2,
                    const real & c3 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real T, const int aDeriv ) const;

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
            mutable real mF[2];

            // 1 + c2*(1-sqrt(T/Tcrit))^2 + c3*(1-sqrt(T/Tcrit))^3
            mutable real mG[3];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_CCR(
                    const real aTcrit,
                    const real & c1,
                    const real & c2,
                    const real & c3 );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real T, const int aDeriv ) const;

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
            mutable real          mF[ 3 ];

            // sqrt( T / T_crit )
            mutable real          mX;

            // 0 : Y     ( help function )
            // 1 : dYdT
            // 2 : d2YdTh
            mutable real          mY[ 3 ];

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            AlphaFunction_PM(
                    const real             aTcrit,
                    const Vector< real > & aCoeffs );

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            void
            eval( const real T, const int aDeriv ) const;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */

#endif //BELFEM_CL_GM_EOS_ALPHAFUNCTION_HPP
