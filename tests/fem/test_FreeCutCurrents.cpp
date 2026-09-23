/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * IWG_Maxwell::set_currents fixes abstract node dof k to aI( k ) and
 * nothing else. An abstract dof beyond the length of aI stays FREE: on a
 * periodic conductor the surplus generator is a legitimate free axial-MMF
 * dof, and fixing it to zero is a modeling error, not a tidy-up. The loop
 * is also tolerant of a list SHORTER than aI, which is the worker-rank
 * shape ( abstract nodes are collected on the master only ).
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Node.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "en_Maxwell_Formulations.hpp"

using namespace belfem;

namespace
{
    // IWG_Maxwell with its abstract-node dof list seeded directly; the list
    // is non-owning, as in the IWG
    class AbstractDofIwg : public fem::IWG_Maxwell
    {
    public:

        AbstractDofIwg() :
            fem::IWG_Maxwell( fem::maxwell::Formulation::HPhi,
                              ModelDimensionality::ThreeD,
                              false, false )
        {
        }

        void
        seed( Cell< fem::Dof * > & aDofs )
        {
            mAbstractNodeDofs = aDofs ;
        }
    };

    // owns the nodes and dofs the list points to; every dof starts free
    class AbstractDofs
    {
        Cell< mesh::Node * > mNodes ;
        Cell< fem::Dof * >   mDofs ;

    public:

        explicit AbstractDofs( const uint aCount )
        {
            mNodes.set_size( aCount, nullptr );
            mDofs.set_size( aCount, nullptr );
            for ( uint k = 0; k < aCount; ++k )
            {
                mNodes( k ) = new mesh::Node( k + 1 );
                mDofs( k ) = new fem::Dof( k + 1, 0, mNodes( k ) );
            }
        }

        AbstractDofs( const AbstractDofs & ) = delete ;
        AbstractDofs( AbstractDofs && ) = delete ;
        AbstractDofs & operator=( const AbstractDofs & ) = delete ;
        AbstractDofs & operator=( AbstractDofs && ) = delete ;

        ~AbstractDofs()
        {
            for ( fem::Dof * tDof : mDofs )
            {
                delete tDof ;
            }
            for ( mesh::Node * tNode : mNodes )
            {
                delete tNode ;
            }
        }

        Cell< fem::Dof * > &
        dofs()
        {
            return mDofs ;
        }

        fem::Dof *
        dof( const uint aIndex )
        {
            return mDofs( aIndex );
        }
    };

    void
    three_currents( Vector< real > & aI )
    {
        aI.set_size( 3 );
        aI( 0 ) = 10.0 ;
        aI( 1 ) = -20.0 ;
        aI( 2 ) = 30.0 ;
    }
}

//------------------------------------------------------------------------------

TEST( FreeCutCurrents, SurplusAbstractDofStaysFree )
{
    AbstractDofs tDofs( 4 );
    AbstractDofIwg tIWG ;
    tIWG.seed( tDofs.dofs() );

    Vector< real > tI ;
    three_currents( tI );
    tIWG.set_currents( tI );

    for ( uint k = 0; k < 3; ++k )
    {
        EXPECT_TRUE( tDofs.dof( k )->is_fixed() ) << "dof " << k ;
        EXPECT_EQ( tDofs.dof( k )->value(), tI( k ) ) << "dof " << k ;
    }

    // the fourth generator: free, and its value untouched
    EXPECT_FALSE( tDofs.dof( 3 )->is_fixed() );
    EXPECT_EQ( tDofs.dof( 3 )->value(), 0.0 );
}

//------------------------------------------------------------------------------

TEST( FreeCutCurrents, ShortDofListFixesNothingExtra )
{
    AbstractDofs tDofs( 1 );
    AbstractDofIwg tIWG ;
    tIWG.seed( tDofs.dofs() );

    Vector< real > tI ;
    three_currents( tI );
    EXPECT_NO_THROW( tIWG.set_currents( tI ) );

    EXPECT_TRUE( tDofs.dof( 0 )->is_fixed() );
    EXPECT_EQ( tDofs.dof( 0 )->value(), tI( 0 ) );
}

//------------------------------------------------------------------------------

TEST( FreeCutCurrents, EmptyCurrentVectorIsANoOp )
{
    AbstractDofs tDofs( 4 );
    AbstractDofIwg tIWG ;
    tIWG.seed( tDofs.dofs() );

    Vector< real > tI ;
    EXPECT_NO_THROW( tIWG.set_currents( tI ) );

    for ( uint k = 0; k < 4; ++k )
    {
        EXPECT_FALSE( tDofs.dof( k )->is_fixed() ) << "dof " << k ;
    }
}
