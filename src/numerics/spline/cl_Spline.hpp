//
// Created by Christian Messe on 2019-01-26.
//

#ifndef BELFEM_SPLINE_HPP
#define BELFEM_SPLINE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "filetools.hpp"

#include "HDF5_Tools.hpp"


namespace belfem
{
    namespace spline
    {
        void
        create_helpmatrix(
                const real & aSize,
                const real & aDeltaX,
                  SpMatrix & aA );
    }

//------------------------------------------------------------------------------

    class Spline
    {
        index_t mNumberOfPoints;
        index_t mNumberOfIntervals;

        real    mXmin;
        real    mXmax;
        real    mDeltaX;
        real    mInvDeltaX;

        Matrix< real >  mData;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        // constructor with helpmatrix
        Spline( const Vector< real > & aX,
                const Vector< real > & aY,
                      SpMatrix       & aA,
                const         real     aXref = 0,
                const         real     aSref = 0,
                const         proc_t   aMasterProc = gNoOwner );

//------------------------------------------------------------------------------

        // parallel constructir
        Spline( const proc_t aMasterProc );

//------------------------------------------------------------------------------

        // constructor that creates empty container
        Spline( const index_t & aN, const real aXmin, const real aXmax );

//------------------------------------------------------------------------------

        inline index_t
        n() const;

//------------------------------------------------------------------------------

        inline real
        x_min() const;

//------------------------------------------------------------------------------

        inline real
        x_max() const;

//------------------------------------------------------------------------------

        inline real
        delta_x() const;

//------------------------------------------------------------------------------

        inline Matrix< real > &
        coefficients();

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        coefficients() const;

//------------------------------------------------------------------------------

        /**
         * interpolate the function
         */
        inline real
        eval( const real aX ) const;

//------------------------------------------------------------------------------
        /**
         * interpolate first derivative
         */
        inline real
        deval( const real aX ) const;

//------------------------------------------------------------------------------
        /**
         * interpolate second derivative
         */
        inline real
        ddeval( const real aX ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate entropy
         */
        inline real
        entropy( const real aX ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate entropy derivative
         */
        inline real
        dentropy( const real aX ) const;

//------------------------------------------------------------------------------

        /**
         * expose matrix
         */
        inline Matrix< real > &
        matrix_data()
        {
            return mData;
        }

//------------------------------------------------------------------------------

        /**
         * expose matrix
         */
        inline const Matrix< real > &
        matrix_data() const
        {
            return mData;
        }

//------------------------------------------------------------------------------

        void
        save(   const string            & aLabel,
                 const string            & aPath,
                 const enum FileMode   aMode=FileMode::NEW );

//------------------------------------------------------------------------------

        void
        save(   hid_t        & aGroup,
                herr_t       & aStatus );

//------------------------------------------------------------------------------

        void
        save_to_database( const string & aDatabase, const string & aLabel );

//------------------------------------------------------------------------------

        void
        update_data(
                SpMatrix             & aHelpMatrix,
                const Vector< real > & aValues,
                const         real     aXref=0.0,
                const         real     aSref=0.0  );

//------------------------------------------------------------------------------

        void
        initialize( const Vector< real > & aX,
                    const Vector< real > & aY,
                          SpMatrix       & aA,
                    const         real     aXref = 0,
                    const         real     aSref = 0 );


//------------------------------------------------------------------------------

        inline index_t
        find_col( const real aX ) const;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        /**
         * check input and return stepsize
         */
        real
        check_input( const Vector< real > & aX,
                     const Vector< real > & aY,
                           SpMatrix       & aA);
//------------------------------------------------------------------------------

        /**
         * right hand side for equation system
         */
        void
        create_rhs(
                const Vector< real > & aY,
                      Vector< real > & aB );

//------------------------------------------------------------------------------

        /**
         * create main coefficients from solution
         */
         void
         create_coeffs(
                const Vector< real > & aX,
                const Vector< real > & aY,
                const Vector< real > & aDYDX );

//------------------------------------------------------------------------------

        /**
         * create main coefficients from solution
         */
        void
        create_coeffs(
                const Vector< real > & aY,
                const Vector< real > & aDYDX );


//------------------------------------------------------------------------------

        /**
         * special case for heat polynomials
         */
         void
         create_entropy( const real aXref, const real & aSref );

//------------------------------------------------------------------------------

        /**
         * synchronize data between procs
         */
        void
        synchronize( const proc_t aMasterProc );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------

    inline index_t
    Spline::n() const
    {
        return mNumberOfPoints;
    }

//------------------------------------------------------------------------------

    inline real
    Spline::x_min() const
    {
        return mXmin;
    }

//------------------------------------------------------------------------------

    inline real
    Spline::x_max() const
    {
        return mXmax;
    }

//------------------------------------------------------------------------------

    inline real
    Spline::delta_x() const
    {
        return mDeltaX;
    }

//------------------------------------------------------------------------------

    Matrix< real > &
    Spline::coefficients()
    {
        return mData;
    }

//------------------------------------------------------------------------------

    const Matrix< real > &
    Spline::coefficients() const
    {
        return mData;
    }

//------------------------------------------------------------------------------

    index_t
    Spline::find_col( const belfem::real aX ) const
    {
        // make sure that x is correct
        BELFEM_ASSERT(
                mXmin <= aX && aX <= mXmax,
                "Value X=%f out of bounds [ %f, %f ]",
                aX, mXmin, mXmax );

        // return column
        return std::floor( ( aX - mXmin ) * mInvDeltaX );
    }

//------------------------------------------------------------------------------

    real
    Spline::eval( const real aX ) const
    {
	    index_t tCol = find_col( aX );

        return (   ( mData( 0, tCol )   * aX
                   + mData( 1, tCol ) ) * aX
                   + mData( 2, tCol ) ) * aX
                   + mData( 3, tCol );
    }

//------------------------------------------------------------------------------

    real
    Spline::deval( const real aX ) const
    {
        auto tCol = find_col( aX );

        return ( ( 3.0 * mData( 0, tCol )   * aX
                 + 2.0 * mData( 1, tCol ) ) * aX
                 +       mData( 2, tCol ) );
    }

//------------------------------------------------------------------------------

    real
    Spline::ddeval( const real aX ) const
    {
        auto tCol = find_col( aX );

        return 6.0 * mData( 0, tCol ) * aX + 2.0 * mData( 1, tCol );
    }

//------------------------------------------------------------------------------

    real
    Spline::entropy( const real aX ) const
    {
        auto tCol = find_col( aX );
        return      ( 1.5 * mData( 0, tCol )   * aX
                   +  2.0 * mData( 1, tCol ) ) * aX
                   +        mData( 2, tCol )   * std::log( aX )
                   +        mData( 4, tCol );
    }

//------------------------------------------------------------------------------

    real
    Spline::dentropy( const real aX ) const
    {
        auto tCol = find_col( aX );
        return   3.0 * mData( 0, tCol ) * aX
               + 2.0 * mData( 1, tCol )
               +       mData( 2, tCol )/aX;
    }

//------------------------------------------------------------------------------
}

#endif // BELFEM_SPLINE_HPP
