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

#ifndef CL_CURVE_HPP
#define CL_CURVE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Segment.hpp"
#include "cl_SideSet.hpp"


namespace belfem
{
    namespace mesh
    {
        class Curve
        {
            const id_t mID ;

            string mLabel = "curve";

            SideSet * mSideSetA ;
            SideSet * mSideSetB ;

            const ElementType mElementType ;

            Cell< Segment * > mSegments ;
            Cell< Node * >    mNodes;
            Cell< Edge * >    mEdges ;

            // flag telling if this loop is closed
            bool mIsClosed ;

            // length coordinate associated with each node
            Vector< real > mS ;

            friend class BfmFile ;

        public :

            Curve( const id_t aID, SideSet * aSideSetA, SideSet * aSideSetB );

            Curve( const id_t aID, const ElementType aElementType );

            ~Curve();

            id_t
            id() const ;

            ElementType
            element_type() const ;

            Cell< Node * > &
            nodes() ;

            const Cell< Node * > &
            nodes() const ;

            Node *
            first() ;

            const Node *
            first() const ;

            Node *
            last() ;

            const Node *
            last() const ;

            Cell< Segment * > &
            segments() ;

            const Cell< Segment * > &
            segments() const ;

            SideSet *
            sideset_a() ;

            SideSet *
            sideset_b() ;

            void
            sideset_a( SideSet * aSideset) ;

            void
            sideset_b( SideSet * aSideset) ;

            bool
            is_closed() const ;

            bool
            is_terminal() const ;

            void
            set_closed_flag( const bool aSwitch );

            void
            save( const std::string & aFileName );

            Vector< real > &
            arclength();

            const Vector< real > &
            arclength() const ;

            real
            arclength( const index_t aIndex ) const ;

            real
            length() const ;

            void
            reverse();

            string &
            label() ;

            const string &
            label() const ;

            Cell< Edge * > &
            edges() ;

            void
            assign_edges();

            size_t
            memory() const;
        };

        inline id_t
        Curve::id() const
        {
            return mID ;
        }

        inline ElementType
        Curve::element_type() const
        {
            return mElementType ;
        }

        inline Cell< Node * > &
        Curve::nodes()
        {
            return mNodes ;
        }

        inline const Cell< Node * > &
        Curve::nodes() const
        {
            return mNodes ;
        }

        inline Node *
        Curve::first()
        {
            return mNodes( 0 );
        }

        inline const Node *
        Curve::first() const
        {
            return mNodes( 0 );
        }

        inline Node *
        Curve::last()
        {
            return mNodes( mNodes.size() -1 );
        }

        inline const Node *
        Curve::last() const
        {
            return mNodes( mNodes.size() -1 );
        }

        inline Cell< Segment * > &
        Curve::segments()
        {
            return mSegments ;
        }

        inline
        const Cell< Segment * > &
        Curve::segments() const
        {
            return mSegments ;
        }

        inline SideSet *
        Curve::sideset_a()
        {
            return mSideSetA ;
        }

        inline SideSet *
        Curve::sideset_b()
        {
            return mSideSetB ;
        }

        inline bool
        Curve::is_closed() const
        {
            return mIsClosed ;
        }

        inline bool
        Curve::is_terminal() const
        {
            return mSideSetA != nullptr && mSideSetB != nullptr ;
        }

        inline void
        Curve::set_closed_flag( const bool aSwitch )
        {
            mIsClosed = aSwitch ;
        }

        inline Vector< real > &
        Curve::arclength()
        {
            return mS ;
        }

        inline const Vector< real > &
        Curve::arclength() const
        {
            return mS ;
        }

        inline real
        Curve::arclength( const index_t aIndex ) const
        {
            return mS( aIndex );
        }

        inline real
        Curve::length() const
        {
            return mS( mNodes.size() -1 ) - mS( 0 );
        }

        inline string &
        Curve::label()
        {
            return mLabel ;
        }

        inline const string &
        Curve::label() const
        {
            return mLabel ;
        }

        inline Cell< Edge * > &
        Curve::edges()
        {
            return mEdges ;
        }

        inline size_t Curve::memory() const
        {
            return sizeof( Curve )
                + mNodes.size() * sizeof( Node * )
                + mEdges.size() * sizeof( Edge * )
                + mSegments.size() * sizeof( Segment * )
                + mS.length() * sizeof( real )
                + mLabel.capacity() * sizeof( char ) ;
        }

    }
}
#endif //CL_CURVE_HPP
