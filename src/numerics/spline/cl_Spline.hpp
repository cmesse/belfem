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

#ifndef BELFEM_SPLINE_HPP
#define BELFEM_SPLINE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "filetools.hpp"

#include "hdf5_tools.hpp"
#include "Spline_Enums.hpp"

namespace belfem
{
    namespace spline
    {
        void
       create_helpmatrix(
              const real & aSize,
              const real & aDeltaX,
              SpMatrix & aA,
              const SplineBC aStartBC = SplineBC::NoCurvature,
              const SplineBC aEndBC = SplineBC::NoCurvature );
    }

//------------------------------------------------------------------------------

    /**
     * @brief Cubic spline on a uniform grid, C2, with natural, parabolic or clamped boundary conditions.
     *
     * @ingroup grp_numerics_spline
     * @see @ref numerics_spline_spline_usage_guide
     */
    class Spline
    {
        const proc_t  mCommRank ;

        index_t mNumberOfPoints;
        index_t mNumberOfIntervals;

        real    mXmin;
        real    mXmax;
        real    mDeltaX;
        real    mInvDeltaX;

        spline::ExtraMode mExtraMode = spline::ExtraMode::None;

        Matrix< real >  mData;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        // constructor with helpmatrix, without boundary conditions
        //
        // aXref carries no default on purpose. With one, a three argument call
        // matches this and the overload below equally well and does not compile.
        // The two agree there anyway: this one forwards NoCurvature at both ends
        // and zero slopes, which is exactly what the other one defaults to.
        Spline( const Vector< real > & aX,
                const Vector< real > & aY,
                      SpMatrix       & aA,
                const         real     aXref,
                const         real     aSref = 0.,
                const         proc_t   aMasterProc = gNoOwner );

        // constructor with helpmatrix
        Spline( const Vector< real > & aX,
                const Vector< real > & aY,
                      SpMatrix       & aA,
                const spline::SplineBC aStartBC = spline::SplineBC::NoCurvature,
                const spline::SplineBC aEndBC = spline::SplineBC::NoCurvature,
                const         real     adYdX0 = 0.,
                const         real     adYdX1 = 0.,
                const         real     aXref = 0.,
                const         real     aSref = 0.,
                const         proc_t   aMasterProc = gNoOwner );

//------------------------------------------------------------------------------

        Spline(
            const string & aFile,
            const string & aLabel,
            const proc_t aMaster=0 );

//------------------------------------------------------------------------------

        Spline( const hid_t aGroup,  const proc_t aMaster );

//------------------------------------------------------------------------------

        // parallel constructor
        Spline( const proc_t aMasterProc );

//------------------------------------------------------------------------------

        // constructor that creates empty container
        Spline( const index_t & aN, const real aXmin, const real aXmax );

//------------------------------------------------------------------------------

        index_t
        n() const;

//------------------------------------------------------------------------------

        real
        x_min() const;

//------------------------------------------------------------------------------

        real
        x_max() const;

//------------------------------------------------------------------------------

        real
        delta_x() const;

//------------------------------------------------------------------------------

        Matrix< real > &
        coefficients();

//------------------------------------------------------------------------------

        const Matrix< real > &
        coefficients() const;

//------------------------------------------------------------------------------

        /**
         * interpolate the function
         */
        real
        eval( const real aX ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate first derivative
         */
        real
        deval( const real aX ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate the function
         */
        real
        eval( const real aX, const index_t aCol ) const;

//------------------------------------------------------------------------------
        /**
         * interpolate first derivative
         */
        real
        deval( const real aX, const index_t aCol ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate second derivative
         */
        real
        ddeval( const real aX ) const;

        real
        ddeval( const real aX, const index_t aCol ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate entropy
         */
        real
        entropy( const real aX ) const;

        real
        entropy( const real aX, const index_t aCol ) const;

//------------------------------------------------------------------------------

        /**
         * interpolate entropy derivative
         */
        real
        dentropy( const real aX ) const;

        real
        dentropy( const real aX, const index_t aCol ) const;

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
         * Declare what the extra row of the coefficient table holds.
         *
         * initialize() and update_data() set this themselves. It is needed by a
         * caller that writes coefficients straight through matrix_data(), which
         * bypasses them: mixing the tables of several splines mixes the extra
         * row along with the rest, so the result carries the same quantity the
         * parts did, but nothing on that path has said so.
         */
        inline void
        set_extra_mode( const spline::ExtraMode aMode )
        {
            mExtraMode = aMode ;
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

        herr_t
        save(   hid_t  & aGroup );

//------------------------------------------------------------------------------

        herr_t
        load(   hid_t & aGroup );

//------------------------------------------------------------------------------

        void
        save_to_database( const string & aDatabase, const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * Recompute spline coefficients for new y-values on the same grid.
         *
         * @note Rank-0 only: this method executes only on rank 0
         *       (mCommRank == 0). It does not broadcast the new coefficients;
         *       other ranks keep the old table. synchronize() is private, so
         *       if all ranks need the update, re-create the spline through a
         *       constructor that takes aMasterProc. This is by design — not
         *       every spline is constructed on all ranks, so automatic
         *       synchronization would be incorrect in general.
         */
        void
        update_data(
                SpMatrix             & aHelpMatrix,
                const Vector< real > & aValues,
                const spline::SplineBC aStartBC=spline::SplineBC::NoCurvature,
                const spline::SplineBC aEndBC=spline::SplineBC::NoCurvature,
                const         real     adYdX0=0.0,
                const         real     adYdX1=0.0,
                const         real     aXref=0.0,
                const         real     aSref=0.0  );

//------------------------------------------------------------------------------

        void
        initialize( const Vector< real > & aX,
                    const Vector< real > & aY,
                          SpMatrix       & aA,
                    const spline::SplineBC aStartBC=spline::SplineBC::NoCurvature,
                    const spline::SplineBC aEndBC=spline::SplineBC::NoCurvature,
                    const         real     adYdX0=0.0,
                    const         real     adYdX1=0.0,
                    const         real     aXref=0.0,
                    const         real     aSref=0.0  );


//------------------------------------------------------------------------------

        inline index_t
        find_col( const real aX ) const;

//------------------------------------------------------------------------------

        void
        create_integral( const real aXref=BELFEM_QUIET_NAN, const real aYref=0.0 );

//------------------------------------------------------------------------------

        real
        integrate( const real aX ) const;

        real
        integrate( const real aX0, const real aX1 ) const ;

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
                      Vector< real > & aB,
                      const spline::SplineBC aStartBC,
                      const spline::SplineBC aEndBC,
                      const         real     adYdX0,
                      const         real     adYdX1);

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

        void
        add_row_to_data();

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

    inline Matrix< real > &
    Spline::coefficients()
    {
        return mData;
    }

//------------------------------------------------------------------------------

    inline const Matrix< real > &
    Spline::coefficients() const
    {
        return mData;
    }

//------------------------------------------------------------------------------

    inline index_t
    Spline::find_col( const real aX ) const
    {
        // out-of-range x: clamp the interval index only; eval() then extrapolates
        // the edge cubic (see spline_usage_guide.md, "Extrapolation behaviour")
        //BELFEM_ASSERT(mXmin <= aX && aX <= mXmax,
        //              "X=%f out of bounds [%f, %f]", aX, mXmin, mXmax);

        const index_t col = static_cast<index_t>((std::min(std::max(aX,mXmin),mXmax) - mXmin) * mInvDeltaX);
        return std::min(col, mNumberOfIntervals - 1);
    }


//------------------------------------------------------------------------------

    inline real
    Spline::eval( const real aX ) const
    {
	    index_t tCol = find_col( aX );

        return (   ( mData( 0, tCol )   * aX
                   + mData( 1, tCol ) ) * aX
                   + mData( 2, tCol ) ) * aX
                   + mData( 3, tCol );
    }

//------------------------------------------------------------------------------

    inline real
    Spline::deval( const real aX ) const
    {
        auto tCol = find_col( aX );

        return ( ( 3.0 * mData( 0, tCol )   * aX
                 + 2.0 * mData( 1, tCol ) ) * aX
                 +       mData( 2, tCol ) );
    }

//------------------------------------------------------------------------------

    inline real
    Spline::eval( const real aX, const index_t aCol ) const
    {
        BELFEM_ASSERT( aCol < mNumberOfIntervals,
            "Spline column %lu out of range, must be < %lu",
            ( long unsigned int ) aCol, ( long unsigned int ) mNumberOfIntervals );

        return (   ( mData( 0, aCol )   * aX
                   + mData( 1, aCol ) ) * aX
                   + mData( 2, aCol ) ) * aX
                   + mData( 3, aCol );
    }

//------------------------------------------------------------------------------

    inline real
    Spline::deval( const real aX, const index_t aCol ) const
    {
        BELFEM_ASSERT( aCol < mNumberOfIntervals,
            "Spline column %lu out of range, must be < %lu",
            ( long unsigned int ) aCol, ( long unsigned int ) mNumberOfIntervals );

        return ( ( 3.0 * mData( 0, aCol )   * aX
                 + 2.0 * mData( 1, aCol ) ) * aX
                 +       mData( 2, aCol ) );
    }

//------------------------------------------------------------------------------

    inline real
    Spline::ddeval( const real aX ) const
    {
        auto tCol = find_col( aX );

        return 6.0 * mData( 0, tCol ) * aX + 2.0 * mData( 1, tCol );
    }

    inline real
    Spline::ddeval( const real aX, const index_t aCol ) const
    {
        BELFEM_ASSERT( aCol < mNumberOfIntervals,
            "Spline column %lu out of range, must be < %lu",
            ( long unsigned int ) aCol, ( long unsigned int ) mNumberOfIntervals );

        return 6.0 * mData( 0, aCol ) * aX + 2.0 * mData( 1, aCol );
    }

//------------------------------------------------------------------------------

    inline real
    Spline::entropy( const real aX ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Entropy,
            "No entropy tables present" );
        auto tCol = find_col( aX );
        return      ( 1.5 * mData( 0, tCol )   * aX
                   +  2.0 * mData( 1, tCol ) ) * aX
                   +        mData( 2, tCol )   * std::log( aX )
                   +        mData( 4, tCol );
    }

    inline real
    Spline::entropy( const real aX, const index_t aCol ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Entropy,
            "No entropy tables present" );

        BELFEM_ASSERT( aCol < mNumberOfIntervals,
            "Spline column %lu out of range, must be < %lu",
            ( long unsigned int ) aCol, ( long unsigned int ) mNumberOfIntervals );


        return      ( 1.5 * mData( 0, aCol )   * aX
                   +  2.0 * mData( 1, aCol ) ) * aX
                   +        mData( 2, aCol )   * std::log( aX )
                   +        mData( 4, aCol );
    }


//------------------------------------------------------------------------------

    inline real
    Spline::dentropy( const real aX ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Entropy,
            "No entropy tables present" );

        auto tCol = find_col( aX );
        return   3.0 * mData( 0, tCol ) * aX
               + 2.0 * mData( 1, tCol )
               +       mData( 2, tCol )/aX;
    }

    inline real
    Spline::dentropy( const real aX, const index_t aCol ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Entropy,
            "No entropy tables present" );

        BELFEM_ASSERT( aCol < mNumberOfIntervals,
            "Spline column %lu out of range, must be < %lu",
            ( long unsigned int ) aCol, ( long unsigned int ) mNumberOfIntervals );


        return   3.0 * mData( 0, aCol ) * aX
               + 2.0 * mData( 1, aCol )
               +       mData( 2, aCol )/aX;
    }

    inline real
    Spline::integrate( const real aX ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Integral,
            "No integration tables present" );

        index_t tCol = find_col( aX );

        return    ((( 0.25 * mData( 0, tCol ) * aX
       + mData( 1, tCol )/3.0   ) * aX
       + 0.5 * mData( 2, tCol ) ) * aX
       +       mData( 3, tCol ) ) * aX
       +       mData( 4, tCol );
    }

    inline real
    Spline::integrate( const real aX0, const real aX1 ) const
    {
        BELFEM_ASSERT( mExtraMode == spline::ExtraMode::Integral,
           "No integration tables present" );

        return this->integrate( aX1 ) - this->integrate( aX0 );
    }

//------------------------------------------------------------------------------
}

#endif // BELFEM_SPLINE_HPP
