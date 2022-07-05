//
// Created by Christian Messe on 05.09.19.
//

#include "GT_globals.hpp"

#ifndef BELFEM_CL_GM_STATAVALS_HPP
#define BELFEM_CL_GM_STATAVALS_HPP

#define BELFEM_STATEVAL_T       0
#define BELFEM_STATEVAL_P       1
#define BELFEM_STATEVAL_V       2

#define BELFEM_STATEVAL_M       3
#define BELFEM_STATEVAL_R       4

#define BELFEM_STATEVAL_ALPHA   5
#define BELFEM_STATEVAL_BETA    6
#define BELFEM_STATEVAL_KAPPA   7

#define BELFEM_STATEVAL_U       8
#define BELFEM_STATEVAL_H       9
#define BELFEM_STATEVAL_S      10
#define BELFEM_STATEVAL_DSDT   11
#define BELFEM_STATEVAL_DSDP   12
#define BELFEM_STATEVAL_CP     13
#define BELFEM_STATEVAL_DHDP   14  // special for tablegas
#define BELFEM_STATEVAL_DCPDT  15
#define BELFEM_STATEVAL_CV     16
#define BELFEM_STATEVAL_GAMMA  17
#define BELFEM_STATEVAL_C      18
#define BELFEM_STATEVAL_HVAP   19
#define BELFEM_STATEVAL_MU     20
#define BELFEM_STATEVAL_LAMBDA 21
#define BELFEM_STATEVAL_PR     22

#define BELFEM_STATEVAL_PI     23  // special parameter for tablegas: pi = 1000 * log10( p )
#define BELFEM_STATEVAL_DMDT   24  // special parameter for tablegas
#define BELFEM_STATEVAL_DMDP   25  // special parameter for tablegas

#define BELFEM_NSTATEVALS      26

#define BELFEM_T_REF 298.15
#define BELFEM_P_REF 100000.0
#define BELFEM_EPSILON_T 1e-7
#define BELFEM_EPSILON_P 1e-4
#define BELFEM_EPSILON_S 1e-4
#define BELFEM_EPSILON_H 1e-4   // for total state
#define BELFEM_EPSILON_X 1e-6   // for molar fraction
#define BELFEM_MAXIT 1000

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Bitset.hpp"

namespace belfem
{
    namespace gasmodels
    {
        class Statevals
        {
            real mStatevals[BELFEM_NSTATEVALS] = { 0.0 };
            Bitset<BELFEM_NSTATEVALS> mStatebits;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Statevals() = default;

//----------------------------------------------------------------------------

            ~Statevals() = default;

//----------------------------------------------------------------------------

            /**
             * update temperature and pressure
             */
            inline void
            update_Tp( const real & aT, const real & aP );

//----------------------------------------------------------------------------

            /**
             * update temperature and inverse density
             */
            inline void
            update_Tv( const real & aT, const real & aV );

//----------------------------------------------------------------------------

            /**
             * update pressure and inverse density
             */
            inline void
            update_pv( const real & aP, const real & aV );

//----------------------------------------------------------------------------

            /**
             * reset all flags of this state
             */
            inline void
            reset();

//----------------------------------------------------------------------------

            /**
             * test if a value is up do date
             */
            inline bool
            test( const index_t & aIndex ) const;

//----------------------------------------------------------------------------

            /**
             * get a value of the state
             */
            inline const real &
            get( const index_t & aIndex ) const;

//----------------------------------------------------------------------------

            /**
             * set a value of the state
             */
            inline void set( const index_t aIndex, const real & aValue );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        bool
        Statevals::test( const index_t & aIndex ) const
        {
            return mStatebits.test( aIndex );
        }

//------------------------------------------------------------------------------
        void
        Statevals::update_Tp( const real & aT, const real & aP )
        {
            BELFEM_ASSERT( aT > 0, "Invalid Temprerature" );
            BELFEM_ASSERT( aP > 0, "Invalid Pressure" );
            if(       aT != mStatevals[ BELFEM_STATEVAL_T ]
                   || aP != mStatevals[ BELFEM_STATEVAL_P ] )
            {
                mStatebits.reset();

                mStatevals[ BELFEM_STATEVAL_T ] = aT;
                mStatevals[ BELFEM_STATEVAL_P ] = aP;

                mStatebits.set( BELFEM_STATEVAL_T );
                mStatebits.set( BELFEM_STATEVAL_P );
            }
        }

//----------------------------------------------------------------------------

        void
        Statevals::update_Tv( const real & aT, const real & aV )
        {
            if(       aT != mStatevals[ BELFEM_STATEVAL_T ]
                   || aV != mStatevals[ BELFEM_STATEVAL_V ] )
            {
                mStatebits.reset();

                mStatevals[ BELFEM_STATEVAL_T ] = aT;
                mStatevals[ BELFEM_STATEVAL_V ] = aV;

                mStatebits.set( BELFEM_STATEVAL_T );
                mStatebits.set( BELFEM_STATEVAL_V );
            }
        }

//----------------------------------------------------------------------------

        void
        Statevals::update_pv( const real & aP, const real & aV )
        {
            if(       aP != mStatevals[ BELFEM_STATEVAL_P ]
                   || aV != mStatevals[ BELFEM_STATEVAL_V ] )
            {
                mStatebits.reset();

                mStatevals[ BELFEM_STATEVAL_P ] = aP;
                mStatevals[ BELFEM_STATEVAL_V ] = aV;

                mStatebits.set( BELFEM_STATEVAL_P );
                mStatebits.set( BELFEM_STATEVAL_V );
            }
        }

//----------------------------------------------------------------------------

        void
        Statevals::reset()
        {
            mStatebits.reset();
        }

//----------------------------------------------------------------------------

        const real &
        Statevals::get( const index_t & aIndex ) const
        {
            BELFEM_ASSERT( aIndex < BELFEM_NSTATEVALS ,
                          "Invalid State index: %u", ( unsigned int ) aIndex );

            return mStatevals[ aIndex ];
        }

//----------------------------------------------------------------------------

        void
        Statevals::set( const index_t aIndex, const real & aValue )
        {
            BELFEM_ASSERT( aIndex < BELFEM_NSTATEVALS ,
                    "Invalid State index: %u", ( unsigned int ) aIndex );

            // set value
            mStatevals[ aIndex ] = aValue;

            // update flag
            mStatebits.set( aIndex );
        }

//----------------------------------------------------------------------------
    } /* namespace gasmodels */
} /* namespace belfem */
#endif //BELFEM_CL_GM_STATAVALS_HPP
