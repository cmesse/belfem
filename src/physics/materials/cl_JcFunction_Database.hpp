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

#ifndef BELFEM_CL_JCFUNCTION_DATABASE_HPP
#define BELFEM_CL_JCFUNCTION_DATABASE_HPP
#include <algorithm>
#include <cmath>
#include "constants.hpp"
#include "commtools.hpp"
#include "cl_JcFunction.hpp"
#include "cl_Database.hpp"
#include "fn_min.hpp"
namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Database-driven Jc function using lookup tables
         *
         * Implements Jc(B, angle, T) or n(B, angle, T) through interpolation
         * of tabulated data stored in HDF5 format.
         *
         * The database uses a structured tensor grid with higher-order finite
         * element interpolation (quad4/9/16 for 2D, hex8/27/64 for 3D).
         *
         * Coordinate system:
         * - Dimension 0: Temperature T [K] (linear scale)
         * - Dimension 1: Magnetic field log10(B) [-] (log scale, B in T)
         * - Dimension 2: Angle θ between n and B [rad]
         * - Values: log10(Jc) or log10(n)
         *
         * The database clamps out-of-range T and high-field inputs to the valid
         * domain and handles periodic angle wrapping (angle ± π); below Bmin a
         * table with a /source group uses the self-field bridge instead of a clamp.
         *
         * Dependencies: normB, angleNxB, T
         *
         * Database file format:
         * - HDF5 file with structured tensor mesh
         * - Dataset name matches material label
         * - See Database class documentation for detailed format
         *
         * Usage:
         * @code
         * JcFunction* jc = factory.create_jc_function(
         *     "path/to/ybco.hdf5", "jc"   // group name; "jc" or "n" for the shipped tables
         * );
         * material->set_jc_function(jc);  // Material takes ownership
         * @endcode
         */
        class JcFunctionDatabase : public JcFunction
        {
            Database * mDatabase ;   //!< Tensor mesh database for interpolation
            const real mTmin ;       //!< Minimum temperature [K]
            const real mTmax ;       //!< Maximum temperature [K]
            const real mBmin ;       //!< Minimum field magnitude [T] (10^log10_min)
            const real mBmax ;       //!< Maximum field magnitude [T] (10^log10_max)
            const real mAngleMin ;   //!< Minimum angle [rad]
            const real mAngleMax ;   //!< Maximum angle [rad]
            const real mInvLn10 = 1.0/std::log( 10. );  //!< 1/ln(10) for log conversion

            //! Measured self-field level, read from the table's embedded
            //! /source rows at B = 0 and angle-averaged. Empty when
            //! the file carries no /source group, in which case the historical
            //! clamp below mBmin is kept and behaviour is unchanged.
            Vector< real > mSelfFieldT ;      //!< temperatures [K], ascending
            Vector< real > mSelfFieldValue ;  //!< angle-free value at B = 0
            bool mHaveSelfField = false ;

        public:

            /**
             * @brief Constructor - loads database from file
             * @param aPath Path to HDF5 database file
             * @param aLabel Group name inside the HDF5 file — "jc" or "n" for the
             *        shipped tables. The self-field bridge keys on this label: only
             *        "jc" reads the icw column and applies the Icw/t_eff calibration;
             *        anything else is treated as an n table.
             *
             * Loads a database-driven Jc or n function from file.
             * The constructor reads the database dimensions and extracts
             * the valid ranges for each coordinate.
             *
             * The database stores:
             * - T in linear scale [K]
             * - log10(B) where B is in [T]
             * - angle in [rad]
             * - log10(Jc) or log10(n) as values
             *
             * Example database creation (Python with h5py):
             * @code
             * # Create 3D grid: T × log10(B) × angle
             * T = np.linspace(4, 90, 50)      # [K]
             * logB = np.linspace(-2, 1, 40)   # log10(B), B in 0.01-10 T
             * angle = np.linspace(0, np.pi, 30)  # [rad]
             * # Store log10(Jc) values on this grid
             * @endcode
             */
            JcFunctionDatabase( const string & aPath, const string & aLabel ) :
                JcFunction(),
                mDatabase( new Database( aPath, aLabel ) ),
                mTmin( mDatabase->min( 0 ) ),
                mTmax( mDatabase->max( 0 ) ),
                mBmin( std::pow( 10., mDatabase->min( 1 ) ) ),
                mBmax( std::pow( 10., mDatabase->max( 1 ) ) ),
                mAngleMin( mDatabase->min( 2 ) ),
                mAngleMax( mDatabase->max( 2 ) )
            {
                // the caller delivers the UNFOLDED angle theta in [ 0, pi ]
                // a folded export covering less than that range
                // would silently wrap or extrapolate — the angle coordinate
                // is not clamped by Database::evaluate. Coverage test, not a
                // width test: [ -pi/2, pi/2 ] has width pi and is still wrong
                BELFEM_ERROR(    mAngleMin <= 1e-6
                              && mAngleMax >= constant::pi - 1e-6,
                    "angular range of database %s spans [ %g, %g ] rad, but the unfolded\n"
                    "bn_angle contract requires coverage of [ 0, pi ]: a folded export\n"
                    "would silently mis-sample theta > pi/2",
                    aLabel.c_str(), ( double ) mAngleMin, ( double ) mAngleMax );

                this->set_dependency( JcParameter::normB );
                this->set_dependency( JcParameter::angleNxB );
                this->set_dependency( JcParameter::T );

                this->load_self_field( aPath, aLabel );
            }

            /**
             * @brief Read the measured B = 0 level from the table's /source
             *
             * The rebuilt tables embed their raw measurement rows as
             * /source/points with /source/columns naming them, and those rows
             * include the self-field scans the log10 B axis structurally
             * cannot hold ( log10( 0 ) is undefined ). Rows at exactly B = 0
             * are binned by temperature and averaged over the stage angles --
             * the average is the whole point, the measured angular scatter
             * there being only a few tenths of a percent.
             *
             * Silently inert on any file without /source: mHaveSelfField stays
             * false and every accessor keeps the historical clamp, so older
             * tables behave exactly as before.
             */
            void
            load_self_field( const string & aPath,
                                                 const string & aLabel )
            {
                // Which raw column carries this function: the jc table is built
                // from the sheet current per width and needs the layer thickness
                // to become a current density; n is dimensionless and is read
                // straight. n is NOT given the jc treatment -- its low-field
                // behaviour is its own.
                const bool tIsJc = aLabel == "jc" ;
                const string tColumn = tIsJc ? "icw" : "n" ;

                Vector< real > tT ;
                Vector< real > tValue ;

    #ifdef BELFEM_HDF5
                if ( comm_rank() == 0 )
                {
                    // raw handle rather than the HDF5 wrapper: /source is a
                    // SIBLING of the group Database opened, and the wrapper
                    // exposes no root handle to probe. Same idiom as
                    // fn_rho_database_is_current.hpp
                    hid_t tFile = H5Fopen( aPath.c_str(),
                                           H5F_ACC_RDONLY, H5P_DEFAULT );

                    if ( tFile >= 0 )
                    {
                        // tri-state, as in rho_database_is_current: only a
                        // definite hit counts, a probe error falls through to
                        // the historical clamp
                        if ( H5Lexists( tFile, "source", H5P_DEFAULT ) > 0 )
                        {
                            hid_t tGroup = H5Gopen2( tFile, "source", H5P_DEFAULT );

                            if ( tGroup >= 0 )
                            {
                                herr_t tStatus = 0 ;

                                if (    hdf5::dataset_exists( tGroup, "points" )
                                     && hdf5::dataset_exists( tGroup, "columns" ) )
                                {
                                    Cell< string > tColumns ;
                                    hdf5::load_strings_from_file(
                                            tGroup, "columns", tColumns, tStatus );

                                    Matrix< real > tPoints ;
                                    hdf5::load_matrix_from_file(
                                            tGroup, "points", tPoints, tStatus );

                                    // columns located by NAME, never by position:
                                    // the raw schema belongs to the data source
                                    index_t tIt = tColumns.size() ;
                                    index_t tIb = tColumns.size() ;
                                    index_t tIv = tColumns.size() ;

                                    for ( index_t k = 0; k < tColumns.size(); ++k )
                                    {
                                        if ( tColumns( k ) == "temperature" ) tIt = k ;
                                        if ( tColumns( k ) == "field" )       tIb = k ;
                                        if ( tColumns( k ) == tColumn )       tIv = k ;
                                    }

                                    // the schema promises columns.size() == points.n_cols,
                                    // but a malformed file must not become out-of-bounds
                                    // access in a release build ( code-audit finding 4 )
                                    if (    tIt < tPoints.n_cols()
                                         && tIb < tPoints.n_cols()
                                         && tIv < tPoints.n_cols()
                                         && tIt < tColumns.size()
                                         && tIb < tColumns.size()
                                         && tIv < tColumns.size() )
                                    {
                                        // reduce UNSCALED first; the calibration
                                        // below needs the curve's own 77.5 K value
                                        this->reduce_self_field_rows(
                                                tPoints, tIt, tIb, tIv,
                                                1.0, tT, tValue );

                                        if ( tIsJc && tT.length() > 0 )
                                        {
                                            // The raw column's UNIT belongs to the data
                                            // vendor, not to us -- sp-ap ships icw in
                                            // A/cm while meta declares A/m, and assuming
                                            // either silently breaks the other
                                            // ( measured: reading A/cm as A/m left J0
                                            // 100x low and the monotonicity guard
                                            // clamped every call ). So the curve is
                                            // CALIBRATED: scaled so its 77.5 K value
                                            // equals the table's own declared level,
                                            // Icw_77p5K_sf / t_eff. No unit assumption
                                            // survives that.
                                            real tScale = 0.0 ;

                                            if ( hdf5::group_exists( tFile, "meta" ) )
                                            {
                                                hid_t tMeta = H5Gopen2( tFile, "meta",
                                                                        H5P_DEFAULT );
                                                if ( tMeta >= 0 )
                                                {
                                                    real tIcw  = 0.0 ;
                                                    real tTeff = 0.0 ;
                                                    if ( hdf5::dataset_exists( tMeta,
                                                            "Icw_77p5K_sf_A_per_m" ) )
                                                    {
                                                        hdf5::load_scalar_from_file(
                                                            tMeta, "Icw_77p5K_sf_A_per_m",
                                                            tIcw, tStatus );
                                                    }
                                                    if ( hdf5::dataset_exists( tMeta,
                                                            "t_eff_m" ) )
                                                    {
                                                        hdf5::load_scalar_from_file(
                                                            tMeta, "t_eff_m",
                                                            tTeff, tStatus );
                                                    }
                                                    H5Gclose( tMeta );

                                                    if ( tIcw > 0.0 && tTeff > 0.0 )
                                                    {
                                                        real tDummy ;
                                                        const real tCurve77p5 =
                                                            interp_curve( tT, tValue,
                                                                          77.5, tDummy );
                                                        if ( tCurve77p5 > 0.0 )
                                                        {
                                                            tScale = ( tIcw / tTeff )
                                                                   / tCurve77p5 ;
                                                        }
                                                    }
                                                }
                                            }

                                            if ( tScale > 0.0 )
                                            {
                                                for ( belfem::index_t k = 0;
                                                      k < tValue.length(); ++k )
                                                {
                                                    tValue( k ) *= tScale ;
                                                }
                                            }
                                            else
                                            {
                                                // no declared level to calibrate against:
                                                // keep the clamp rather than guess a unit
                                                tT.set_size( 0, 0.0 );
                                                tValue.set_size( 0, 0.0 );
                                            }
                                        }
                                    }
                                }
                                H5Gclose( tGroup );
                            }
                        }
                        H5Fclose( tFile );
                    }
                }
    #endif

                // the LENGTH is decided on rank 0 and must reach the others
                // before the payload: a few dozen entries, so broadcast rather
                // than the chunked share/receive pair
                uint tN = ( comm_rank() == 0 ) ? ( uint ) tT.length() : 0 ;
                broadcast( tN );

                if ( tN == 0 ) return ;

                if ( comm_rank() != 0 )
                {
                    tT.set_size( tN, 0.0 );
                    tValue.set_size( tN, 0.0 );
                }
                broadcast( tT );
                broadcast( tValue );

                mSelfFieldT = tT ;
                mSelfFieldValue = tValue ;
                mHaveSelfField = true ;
            }

            /**
             * @brief Bin the B = 0 rows by temperature and average over angle
             *
             * The angular average IS the physics: at zero field there is no
             * field direction, and the measured rows scatter only a few
             * tenths of a percent across the stage angles. A temperature is
             * a distinct bin when it differs by more than the tolerance,
             * which is loose enough to merge one setpoint's jitter
             * ( 77.48-77.51 K ) and tight enough to keep 2.5 K spacing apart.
             */
            void
            reduce_self_field_rows(
                    const Matrix< real > & aPoints,
                    const index_t aColT,
                    const index_t aColB,
                    const index_t aColV,
                    const real aScale,
                    Vector< real > & aT,
                    Vector< real > & aValue ) const
            {
                const index_t tNumRows = aPoints.n_rows() ;

                // one setpoint's temperature jitter is a few hundredths of a
                // kelvin; the setpoint spacing is 2.5 K
                constexpr real tTol = 0.25 ;

                Cell< real > tT ;
                Cell< real > tSum ;
                Cell< uint > tCount ;

                for ( index_t i = 0; i < tNumRows; ++i )
                {
                    const real tB = aPoints( i, aColB ) ;
                    const real tTi = aPoints( i, aColT ) ;
                    const real tVi = aPoints( i, aColV ) ;

                    // exact self-field rows only. The source carries NaNs in
                    // rows where a fit did not converge; those must not enter
                    // an average
                    if ( tB != 0.0 ) continue ;
                    if ( ! std::isfinite( tTi ) ) continue ;
                    if ( ! std::isfinite( tVi ) ) continue ;
                    if ( tVi <= 0.0 ) continue ;

                    index_t tBin = tT.size() ;
                    for ( index_t k = 0; k < tT.size(); ++k )
                    {
                        if ( std::abs( tT( k ) - tTi ) < tTol ) { tBin = k ; break ; }
                    }

                    if ( tBin == tT.size() )
                    {
                        tT.push( tTi );
                        tSum.push( tVi );
                        tCount.push( 1 );
                    }
                    else
                    {
                        tSum( tBin ) += tVi ;
                        ++tCount( tBin ) ;
                    }
                }

                const index_t tNumBins = tT.size() ;
                if ( tNumBins < 2 ) return ;

                // ascending in T, as self_field() assumes
                Vector< index_t > tOrder( tNumBins );
                for ( index_t k = 0; k < tNumBins; ++k ) tOrder( k ) = k ;

                for ( index_t a = 0; a + 1 < tNumBins; ++a )
                {
                    for ( index_t b = a + 1; b < tNumBins; ++b )
                    {
                        if ( tT( tOrder( b ) ) < tT( tOrder( a ) ) )
                        {
                            const index_t tSwap = tOrder( a ) ;
                            tOrder( a ) = tOrder( b ) ;
                            tOrder( b ) = tSwap ;
                        }
                    }
                }

                aT.set_size( tNumBins, 0.0 );
                aValue.set_size( tNumBins, 0.0 );

                for ( index_t k = 0; k < tNumBins; ++k )
                {
                    const index_t j = tOrder( k ) ;
                    aT( k ) = tT( j ) ;
                    aValue( k ) = aScale * tSum( j )
                                / static_cast< real >( tCount( j ) ) ;
                }
            }

            /**
             * @brief Destructor - deletes the database
             */
            ~JcFunctionDatabase() override
            {
                delete mDatabase ;
            }

            //! smallest value the table can produce; values are stored as
            //! log10, so exponentiate the raw minimum
            real
            min_value() const override
            {
                return std::pow( 10., min( mDatabase->values() ) ) ;
            }

            //--------------------------------------------------------------
            // low-field transition
            //--------------------------------------------------------------
            //
            // Below mBmin the table has no data and the historical behaviour
            // was a hard clamp: jc constant in B, d/dB identically zero. Two
            // problems, both measured. The clamp is wrong by up to 23 % near
            // Tc ( the deck's whole self-field range sits under the floor ),
            // and the tangent jumps from 0 to the spline slope the instant
            // |B| crosses mBmin, which is a Newton tangent.
            //
            // The replacement is a quadratic in LINEAR B -- the same move
            // Copper::create_kohler makes for magnetoresistance, where the
            // physical boundary condition is built into the polynomial's form
            // rather than imposed afterwards:
            //
            //     jc(B) = J0 + c B + b B^2
            //
            // anchored at the MEASURED angle-free self-field level J0(T) and
            // matching the table's own value and slope at mBmin, so the join
            // is C1 and nothing at or above mBmin changes.
            //
            // Two properties are not conveniences, they are what the data and
            // the literature require:
            //
            //  * The linear term is KEPT. The Kim family used throughout the
            //    HTS literature -- jc = jc0 / ( 1 + |B_eff|/B0 )^alpha, e.g.
            //    Riva et al. 2023 Eq. 2 and Denis et al. 2026 Eq. 11, with
            //    Messe et al. 2023 §2.6 naming Kim 1962 -- has
            //    d jc / d|B| = -jc0 alpha / B0 at B = 0: FINITE and non-zero.
            //    A form with zero slope at the origin would contradict it.
            //  * J0 is ANGLE-FREE, and by measurement rather than by fiat. In
            //    the Kim family the anisotropy enters only inside |B_eff|, so
            //    it vanishes with the field; the embedded source rows agree,
            //    scattering only 0.25-0.73 % across 55 stage angles at B = 0
            //    while the table spreads ~4 % at mBmin. Carrying the mBmin
            //    angular shape down to the axis was measured to overshoot the
            //    permitted spread by an order of magnitude, which is why J0
            //    cannot be extrapolated from the table and is read instead.
            //
            // Monotone as long as J0 >= jc(mBmin, theta) for every theta,
            // which the measurements satisfy; see the guard in eval().

            //! measured angle-free value at B = 0 and its dT, by linear
            //! interpolation on the stored curve. Clamped at both ends: the
            //! curve spans the measured temperature range, not the table's
            static real
            interp_curve( const Vector< real > & aT,
                          const Vector< real > & aValue,
                          const real T, real & dJ0dT )
            {
                const index_t tN = aT.length() ;

                if ( T <= aT( 0 ) )
                {
                    dJ0dT = 0.0 ;
                    return aValue( 0 ) ;
                }
                if ( T >= aT( tN - 1 ) )
                {
                    dJ0dT = 0.0 ;
                    return aValue( tN - 1 ) ;
                }

                index_t k = 1 ;
                while ( k < tN - 1 && aT( k ) < T ) ++k ;

                const real tDT = aT( k ) - aT( k - 1 ) ;
                dJ0dT = ( aValue( k ) - aValue( k - 1 ) ) / tDT ;

                return aValue( k - 1 )
                     + dJ0dT * ( T - aT( k - 1 ) ) ;
            }

            real
            self_field( const real T, real & dJ0dT ) const
            {
                // every internal caller sits behind mHaveSelfField; a future
                // external caller must not turn an empty curve into an
                // out-of-bounds read ( code-audit, Grok finding 4 )
                BELFEM_ASSERT( mHaveSelfField,
                    "self_field() called on a table without a self-field curve" );
                return interp_curve( mSelfFieldT, mSelfFieldValue, T, dJ0dT ) ;
            }

            //! the table's own value and d/dB at the low-field edge
            void
            edge_value_and_slope( const real angle, const real T,
                                  real & V, real & S ) const
            {
                const real tTc = std::clamp( T, mTmin, mTmax ) ;
                const real tU  = std::log( mBmin ) * mInvLn10 ;
                const real tTh = this->wrap_angle( angle ) ;

                // the VALUE goes through eval() rather than the raw spline so
                // that anything eval() applies to the stored field ( e.g. a
                // lift-factor reference shift ) reaches the bridge too; the
                // log-slope is invariant under such an additive shift, so the
                // derivative may read the spline directly
                V = this->eval( mBmin, angle, T ) ;
                S = V * mDatabase->evaluate_derivy( tTc, tU, tTh ) / mBmin ;
            }

            /**
             * @brief Wrap the caller's angle into the table window
             *
             * theta in [ 0, pi ] passes through untouched; values within
             * the same 1e-6 tolerance the ctor coverage guard uses are
             * snapped onto the boundary — a table whose endpoints sit an
             * ulp inside 0 or pi must not swap poles at exactly theta = 0
             * or pi. Anything further out uses the +-pi periodicity.
             */
            real
            wrap_angle( const real angle ) const
            {
                return angle < mAngleMin - 1e-6 ? angle + constant::pi :
                       angle > mAngleMax + 1e-6 ? angle - constant::pi :
                       std::clamp( angle, mAngleMin, mAngleMax );
            }

            /**
             * @brief Evaluate Jc or n at given field, angle, and temperature
             * @param normB Magnetic field magnitude [T]
             * @param angle Angle between surface normal and field [rad]
             * @param T Temperature [K]
             * @return Jc [A/m²] or n [-] interpolated from database
             *
             * Performs 3D tensor interpolation:
             * 1. Clamps T to [Tmin, Tmax]
             * 2. |B| ≥ Bmin: clamps B to [Bmin, Bmax] and converts to log10(B).
             *    |B| < Bmin: if the file carries a /source group, evaluates the
             *    C1 self-field bridge jc = J0(T) + c·B + b·B² ( see the
             *    low-field transition notes above ); otherwise clamps to Bmin
             * 3. Wraps angle into [angle_min, angle_max] using ±π periodicity;
             *    the caller delivers the unfolded θ ∈ [0, π], which
             *    passes through untouched — θ = π reads the stored 180° node,
             *    honoring the measured asymmetry about 90°
             * 4. Interpolates log10(Jc) or log10(n) using FEM shape functions
             * 5. Returns 10^(interpolated value)
             *
             * Above Bmax and outside [Tmin, Tmax] inputs are clamped to the
             * database bounds; below Bmin only tables without /source clamp.
             */
            real
            eval( const real normB, const real angle, const real T ) const override
            {
                if ( mHaveSelfField && normB < mBmin )
                {
                    real tJ0dT ;
                    const real tJ0 = this->self_field( T, tJ0dT ) ;
                    real tV, tS ;
                    this->edge_value_and_slope( angle, T, tV, tS ) ;

                    // a table whose edge value exceeds the measured
                    // self-field level cannot be bridged monotonically --
                    // that is an inconsistent table, not a physical case
                    if ( tJ0 < tV ) return tV ;

                    const real tD = tV - tJ0 ;
                    real tC = 2. * tD / mBmin - tS ;
                    real tB ;

                    // monotone limiter ( code-audit finding 2 ): with a steep
                    // table slope and a small lift the exact-match quadratic
                    // turns non-monotone ( measured: 60 % of (T,theta) pairs,
                    // worst interior hump 6.1 % near Tc ). Where c would be
                    // positive, fall to the pure parabola: value match kept,
                    // monotone by construction, tangents endpoint-exact; the
                    // cost is a bounded one-sided slope kink at the join, at
                    // angles where the lift is nearly zero anyway
                    if ( tC > 0.0 )
                    {
                        tC = 0.0 ;
                        tB = tD / ( mBmin * mBmin ) ;
                    }
                    else
                    {
                        tB = ( tS - tD / mBmin ) / mBmin ;
                    }

                    return tJ0 + normB * ( tC + normB * tB ) ;
                }

                return std::pow( 10.,
                    mDatabase->evaluate(
                        std::clamp( T, mTmin, mTmax ),
                        std::log( std::clamp( normB, mBmin, mBmax ) ) * mInvLn10,
                        this->wrap_angle( angle ) ) );
            }

            /**
             * @brief d(value)/d|B| from the spline, clamp-consistent
             *
             * The table stores f = log10(value) over ( T, u = log10 B, θ ),
             * so with value = 10^f and du/dB = 1/( B ln10 ):
             *
             *     d(value)/dB = value · ln10 · (∂f/∂u) · du/dB
             *                 = value · (∂f/∂u) / B          — ln10 CANCELS.
             *
             * ( Verified independently by both audit voices 2026-08-13; an
             * earlier plan note carried a spurious ln10 here. ) When |B| is
             * outside the table window, eval() returns the CLAMPED value,
             * which is constant in |B| — the consistent tangent is exactly
             * zero, and the derivative implements that clamp decision
             * itself rather than trusting a caller to.
             */
            real
            deval_dB( const real normB, const real angle, const real T ) const override
            {
                if ( mHaveSelfField && normB < mBmin )
                {
                    real tJ0dT ;
                    const real tJ0 = this->self_field( T, tJ0dT ) ;
                    real tV, tS ;
                    this->edge_value_and_slope( angle, T, tV, tS ) ;
                    if ( tJ0 < tV ) return 0.0 ;

                    const real tD = tV - tJ0 ;
                    real tC = 2. * tD / mBmin - tS ;
                    real tB ;

                    // monotone limiter -- same branch as eval(), so the
                    // tangent differentiates the returned value exactly
                    if ( tC > 0.0 )
                    {
                        tC = 0.0 ;
                        tB = tD / ( mBmin * mBmin ) ;
                    }
                    else
                    {
                        tB = ( tS - tD / mBmin ) / mBmin ;
                    }

                    // -> tS at mBmin ( C1 ) on the exact branch, -> tC
                    // finite at B = 0 ( Kim ) where the limiter is idle
                    return tC + 2. * tB * normB ;
                }

                if ( normB <= mBmin || normB >= mBmax ) return 0.0 ;

                return this->eval( normB, angle, T )
                    * mDatabase->evaluate_derivy(
                          std::clamp( T, mTmin, mTmax ),
                          std::log( normB ) * mInvLn10,
                          this->wrap_angle( angle ) )
                    / normB ;
            }

            /**
             * @brief d(value)/dθ from the spline
             *
             * θ is stored linearly, so the log10 storage contributes the
             * full ln10: d(value)/dθ = value · ln10 · ∂f/∂θ. The ±π wrap is
             * a shift ( dθ_wrapped/dθ = 1 ), applied identically to eval().
             */
            real
            deval_dbeta( const real normB, const real angle, const real T ) const override
            {
                if ( mHaveSelfField && normB < mBmin )
                {
                    // J0 is angle-free, so only c and b carry theta:
                    //   d/dtheta = Vth ( 2 beta - beta^2 )
                    //            + Sth mBmin ( beta^2 - beta ),   beta = B/mBmin
                    // EXACT at both ends -- the second bracket vanishes at
                    // beta = 0 and beta = 1, so this reproduces 0 at zero
                    // field and the table's own dV/dtheta at mBmin.
                    //
                    // DECLARED APPROXIMATION: Sth = d/dtheta ( dV/dB ) needs a
                    // MIXED spline partial that Database does not expose, and
                    // is taken as zero. The dropped term is bounded by
                    // |Sth| mBmin/4 ( the extremum of beta^2 - beta ), i.e. it
                    // is worst mid-interval and exactly zero where the tangent
                    // has to match. Supplying mixed partials would remove it.
                    real tJ0dT ;
                    const real tJ0 = this->self_field( T, tJ0dT ) ;
                    real tV, tS ;
                    this->edge_value_and_slope( angle, T, tV, tS ) ;

                    const real tVth0 = tV * mDatabase->evaluate_derivz(
                            std::clamp( T, mTmin, mTmax ),
                            std::log( mBmin ) * mInvLn10,
                            this->wrap_angle( angle ) ) / mInvLn10 ;

                    // inconsistent table: eval() degrades to the historical
                    // clamp, so the tangent must be the CLAMP's tangent --
                    // the value still varies with theta there, and a zero
                    // here would hand Newton a derivative that does not
                    // differentiate the value ( code-audit finding 1 )
                    if ( tJ0 < tV ) return tVth0 ;

                    const real tBeta = normB / mBmin ;

                    // limiter branch: value = J0 + D(theta) beta^2, so the
                    // theta shape factor is beta^2 instead of beta( 2-beta );
                    // both are endpoint-exact ( 0 at B = 0, tVth0 at mBmin )
                    const real tD = tV - tJ0 ;
                    const bool tLimited = ( 2. * tD / mBmin - tS ) > 0.0 ;

                    return tVth0 * ( tLimited ? tBeta * tBeta
                                              : tBeta * ( 2. - tBeta ) ) ;
                }

                return this->eval( normB, angle, T )
                    * mDatabase->evaluate_derivz(
                          std::clamp( T, mTmin, mTmax ),
                          std::log( std::clamp( normB, mBmin, mBmax ) ) * mInvLn10,
                          this->wrap_angle( angle ) )
                    / mInvLn10 ;
            }

            /**
             * @brief d(value)/dT from the spline, clamp-consistent
             *
             * T is stored linearly: d(value)/dT = value · ln10 · ∂f/∂T.
             * Outside the temperature window the clamped value is constant
             * in T — the consistent tangent is zero.
             */
            real
            deval_dT( const real normB, const real angle, const real T ) const override
            {
                if ( mHaveSelfField && normB < mBmin )
                {
                    //   d/dT = J0' + ( VT - J0' )( 2 beta - beta^2 )
                    //               + ST mBmin ( beta^2 - beta )
                    // exact at beta = 0 ( -> J0', the measured self-field
                    // slope ) and at beta = 1 ( -> VT, the table's own dV/dT ).
                    // ST is dropped on the same bound as Sth in deval_dbeta.
                    real tJ0dT ;
                    const real tJ0 = this->self_field( T, tJ0dT ) ;
                    real tV, tS ;
                    this->edge_value_and_slope( angle, T, tV, tS ) ;

                    const real tVT = ( T <= mTmin || T >= mTmax ) ? 0.0 :
                          tV * mDatabase->evaluate_derivx(
                              T,
                              std::log( mBmin ) * mInvLn10,
                              this->wrap_angle( angle ) ) / mInvLn10 ;

                    // inconsistent table: clamp value, clamp tangent -- see
                    // deval_dbeta. dT is a live thermal Newton tangent
                    if ( tJ0 < tV ) return tVT ;

                    const real tBeta = normB / mBmin ;

                    // shape factor follows eval()'s limiter branch
                    const real tD = tV - tJ0 ;
                    const bool tLimited = ( 2. * tD / mBmin - tS ) > 0.0 ;
                    const real tShape = tLimited ? tBeta * tBeta
                                                 : tBeta * ( 2. - tBeta ) ;

                    return tJ0dT + ( tVT - tJ0dT ) * tShape ;
                }

                if ( T <= mTmin || T >= mTmax ) return 0.0 ;

                return this->eval( normB, angle, T )
                    * mDatabase->evaluate_derivx(
                          T,
                          std::log( std::clamp( normB, mBmin, mBmax ) ) * mInvLn10,
                          this->wrap_angle( angle ) )
                    / mInvLn10 ;
            }
        };
    }
}
#endif //BELFEM_CL_JCFUNCTIONDATABASE_HPP