//
// Created by Christian Messe on 10.07.20.
//
#include "commtools.hpp"
#include "cl_SolverWrapper.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        Wrapper::Wrapper( const string & aLabel ) :
            mMyRank(   gComm.rank() ),
            mCommSize( gComm.size() ),
            mLabel( aLabel )
        {

        }

//------------------------------------------------------------------------------
        void
        Wrapper::initialize()
        {
            // make sure that this wrapper has not been initialized yet
            BELFEM_ERROR( ! mIsInitialized,
                    "Wrapper for %s has already been initialized",
                         mLabel.c_str() );

            // set the initialized flag
            mIsInitialized = true ;
        }

//------------------------------------------------------------------------------

        void
        Wrapper::initialize( SpMatrix & aMatrix,
                    const SymmetryMode aSymmetryMode,
                    const int aNumRhsColumns )
        {
            BELFEM_ERROR( false, "initialize() is not implemented for %s",
                         mLabel.c_str() );
        }

//------------------------------------------------------------------------------

        void
        Wrapper::free()
        {
            // set the initialized flag
            mIsInitialized = false ;
        }

//------------------------------------------------------------------------------

        const string &
        Wrapper::label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        void
        Wrapper::solve( SpMatrix & aMatrix,
                           Vector <real> & aLHS,
                           Vector <real> & aRHS )
        {
            BELFEM_ERROR( false, "solve() for vectors is not implemented for %s",
                         mLabel.c_str() );
        }

//------------------------------------------------------------------------------

        void
        Wrapper::solve( SpMatrix & aMatrix,
                           Matrix <real> & aLHS,
                           Matrix <real> & aRHS )
        {
            BELFEM_ERROR( false, "solve() for matrices is not implemented for %s",
                         mLabel.c_str() );
        }

//------------------------------------------------------------------------------

        void
        Wrapper::mat2vec( const Matrix< real > & aM,
                 Vector< real > & aV )
        {
            uint n = aM.n_rows();
            uint m = aM.n_cols();
            uint c = 0;

            aV.set_size( n*m );

            for( uint j=0; j<m; ++j )
            {
                for( uint i=0; i<n; ++i )
                {
                    aV( c++ ) = aM( i, j );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Wrapper::vec2mat( const Vector< real > & aV,
                 Matrix< real > & aM )
        {
            uint n = aM.n_rows();
            uint m = aM.n_cols();
            uint c = 0;

            for( uint j=0; j<m; ++j )
            {
                for( uint i=0; i<n; ++i )
                {
                    aM( i, j ) = aV( c++ );
                }
            }
        }

//------------------------------------------------------------------------------

        real
        Wrapper::get_determinant() const
        {
            BELFEM_ERROR( false,
                    "get_determinant() is not implemented for solver %s",
                mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------
    }
}
