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
#include <memory>
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
         * domain and handles periodic angle wrapping (angle ± π). Below Bmin,
         * n and a jc table with a measured self-field curve continue the table
         * C1 and non-increasing in |B| ( see the low-field notes below ); a jc
         * table without that curve keeps the clamp.
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
            std::unique_ptr< Database > mDatabase ;   //!< Tensor mesh database for interpolation
            const real mTmin ;       //!< Minimum temperature [K]
            const real mTmax ;       //!< Maximum temperature [K]
            const real mBmin ;       //!< Minimum field magnitude [T] (10^log10_min)
            const real mBmax ;       //!< Maximum field magnitude [T] (10^log10_max)
            const real mAngleMin ;   //!< Minimum angle [rad]
            const real mAngleMax ;   //!< Maximum angle [rad]
            const real mInvLn10 = 1.0/std::log( 10. );  //!< 1/ln(10) for log conversion

            //! Measured self-field level of a jc table, read from the table's
            //! embedded /source rows at B = 0 and angle-averaged, projected to
            //! be non-increasing in T. Empty when the file carries no usable
            //! /source rows; such a jc table keeps the clamp below mBmin.
            Vector< real > mSelfFieldT ;      //!< temperatures [K], ascending
            Vector< real > mSelfFieldValue ;  //!< angle-free value at B = 0
            Vector< real > mSelfFieldSlope ;  //!< its dT at the knots ( Fritsch & Butland 1984, Eq. 5 )
            bool mHaveSelfField = false ;

            //! any label but "jc" is an n table
            bool mIsN = false ;

            //! width of the smooth maximum that forms the jc anchor, relative
            //! to the requirement b
            static constexpr real mSmoothWidth = 1e-3 ;

            //! the table and its derivatives on the face log10 B = log10 Bmin,
            //! at the clamped temperature and the wrapped angle. d/dT terms are
            //! zero outside [ Tmin, Tmax ], where the value is held
            struct Edge
            {
                real V   = 0.0 ;   //!< value at Bmin
                real S   = 0.0 ;   //!< d value / dB at Bmin
                real V_T = 0.0 ;
                real V_A = 0.0 ;   //!< d / d angle
                real S_T = 0.0 ;
                real S_A = 0.0 ;
            };


        public:

            /**
             * @brief Constructor - loads database from file
             * @param aPath Path to HDF5 database file
             * @param aLabel Group name inside the HDF5 file — "jc" or "n" for the
             *        shipped tables. Only "jc" reads the measured self-field curve
             *        ( icw column, Icw/t_eff calibration ); anything else is an n
             *        table, continued below Bmin from the table alone.
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
                mDatabase( open_checked( aPath, aLabel ) ),
                mTmin( mDatabase->min( 0 ) ),
                mTmax( mDatabase->max( 0 ) ),
                mBmin( std::pow( 10., mDatabase->min( 1 ) ) ),
                mBmax( std::pow( 10., mDatabase->max( 1 ) ) ),
                mAngleMin( mDatabase->min( 2 ) ),
                mAngleMax( mDatabase->max( 2 ) )
            {
                this->set_dependency( JcParameter::normB );
                this->set_dependency( JcParameter::angleNxB );
                this->set_dependency( JcParameter::T );

                mIsN = aLabel != "jc" ;

                // n continues from the table alone; jc anchors on the measured
                // self-field curve where the file carries one
                if ( ! mIsN )
                {
                    this->load_self_field( aPath, aLabel );
                    if ( mHaveSelfField )
                    {
                        this->finalize_self_field();
                    }
                }
            }

            /**
             * @brief Open the table and refuse one the reader cannot serve
             *
             * Runs before any member reads the database, and owns it until every
             * check has passed, so a refusal does not leak it. Every rank holds
             * the same values ( Database broadcasts them ), so every rank refuses
             * together.
             */
            static Database *
            open_checked( const string & aPath, const string & aLabel )
            {
                std::unique_ptr< Database > tDatabase( new Database( aPath, aLabel ) );

                BELFEM_ERROR( tDatabase->dimension() == 3,
                    "database %s in %s has %u axes; a jc or n table needs T, log10 B and theta",
                    aLabel.c_str(), aPath.c_str(), ( unsigned int ) tDatabase->dimension() );

                // the caller delivers the UNFOLDED angle theta in [ 0, pi ]; a
                // folded export covering less than that range would silently
                // wrap or extrapolate. Coverage test, not a width test:
                // [ -pi/2, pi/2 ] has width pi and is still wrong
                BELFEM_ERROR(    tDatabase->min( 2 ) <= 1e-6
                              && tDatabase->max( 2 ) >= constant::pi - 1e-6,
                    "angular range of database %s spans [ %g, %g ] rad, but the unfolded\n"
                    "bn_angle contract requires coverage of [ 0, pi ]: a folded export\n"
                    "would silently mis-sample theta > pi/2",
                    aLabel.c_str(), ( double ) tDatabase->min( 2 ), ( double ) tDatabase->max( 2 ) );

                // the continuation below Bmin is monotone only if the table does
                // not rise with |B| on its lowest-field face. Certified over the
                // whole face ( Bernstein bound ), not sampled; the tolerance
                // accepts the rounding residue of a zero slope
                const real tSlopeBound = tDatabase->upper_bound_derivy_on_lower_y_face() ;
                BELFEM_ERROR( tSlopeBound <= 1e-10,
                    "database %s in %s cannot be certified not to rise with |B| at its "
                    "lowest field ( bound on d log10 value / d log10 B: %g ), so it "
                    "cannot be continued monotonically below Bmin",
                    aLabel.c_str(), aPath.c_str(), ( double ) tSlopeBound );

                if ( aLabel != "jc" )
                {
                    // n below Bmin is at least its value at Bmin, which must
                    // keep the power law a power law
                    const real tFloor = tDatabase->lower_bound_on_lower_y_face() ;
                    BELFEM_ERROR( tFloor >= std::log10( 1.01 ),
                        "database %s in %s reaches n = %g at its lowest field; n must be at least 1.01",
                        aLabel.c_str(), aPath.c_str(), ( double ) std::pow( 10., tFloor ) );
                }

                return tDatabase.release() ;
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
             * false and the jc table keeps the clamp below Bmin.
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
                                    // access in a release build
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
            ~JcFunctionDatabase() override = default ;

            //! smallest value the table can produce; values are stored as
            //! log10, so exponentiate the raw minimum
            real
            min_value() const override
            {
                return std::pow( 10., min( mDatabase->values() ) ) ;
            }


            //--------------------------------------------------------------
            // below Bmin
            //--------------------------------------------------------------
            //
            // The table's field axis is log10 B and ends at Bmin; B = 0 is out
            // of its reach. Below Bmin the value is continued in LINEAR B from
            // the table's own value V and slope S = dV/dB at Bmin, so the join
            // is C1 and nothing at or above Bmin changes. The table does not
            // rise with |B| on that face ( certified at load ), so S <= 0.
            //
            // n: the parabola with value V and slope S at Bmin and slope 0 at
            // B = 0, n = V + S ( B^2 - Bmin^2 ) / ( 2 Bmin ). n depends on |B|,
            // so a zero slope at the origin is the smooth choice.
            //
            // jc: a cubic Hermite from an anchor J at B = 0 to ( V, S ). The
            // anchor is the MEASURED angle-free self-field level a( T ) of the
            // table's /source rows. In the Kim family used throughout the HTS
            // literature, jc = jc0 / ( 1 + |B_eff| / B0 )^alpha ( Riva et al.
            // 2023 Eq. 2, Denis et al. 2026 Eq. 11; Messe et al. 2023 §2.6 names
            // Kim 1962 ), the slope at B = 0 is finite and non-zero, and the
            // anisotropy vanishes with the field. A monotone C1 join needs
            // J >= b = V - S Bmin / 3. The measured level sits within a few
            // percent of b over most of the range and below it at some angles,
            // and outside the measured temperatures its linear continuation
            // can fall well short ( at low T and near Tc ), so the anchor is
            // J = smax( a, b ): the measured level where it clears b, lifted
            // smoothly to b elsewhere, and angle-dependent only where lifted.
            // With Delta = ( V - J ) / Bmin < 0 and beta = S / Delta in
            // [ 0, 3 ), the slope at B = 0 is s0 = 2 Delta - S for beta <= 2
            // ( the exact quadratic, a finite Kim slope ) and 0 for beta > 2,
            // which is common wherever the anchor is lifted. Either way
            // ( alpha, beta ) = ( s0 / Delta, beta ) lies in the square
            // 0 <= alpha, beta <= 3, which suffices for a monotone cubic
            // ( Fritsch & Carlson 1980, §4 ), so the join does not rise with
            // |B| at any point. At beta = 2 the value is continuous but its T
            // and angle tangents are not. A jc table without a measured curve
            // keeps the clamp below Bmin.

            //! linear interpolation on the raw curve, clamped at both ends;
            //! used only to calibrate the curve at 77.5 K
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

            //! isotonic projection ( the running maximum from high to low T )
            //! of the measured curve, and its knot slopes ( Fritsch & Butland
            //! 1984, Eq. 5; one-sided secants at the ends ), which keep every segment
            //! monotone ( Fritsch & Carlson 1980, §4 )
            void
            finalize_self_field()
            {
                const index_t tN = mSelfFieldT.length() ;

                for ( index_t k = tN - 1; k > 0; --k )
                {
                    mSelfFieldValue( k - 1 ) = std::max( mSelfFieldValue( k - 1 ), mSelfFieldValue( k ) );
                }

                mSelfFieldSlope.set_size( tN, 0.0 );
                for ( index_t k = 1; k + 1 < tN; ++k )
                {
                    const real h0 = mSelfFieldT( k ) - mSelfFieldT( k - 1 ) ;
                    const real h1 = mSelfFieldT( k + 1 ) - mSelfFieldT( k ) ;
                    const real d0 = ( mSelfFieldValue( k ) - mSelfFieldValue( k - 1 ) ) / h0 ;
                    const real d1 = ( mSelfFieldValue( k + 1 ) - mSelfFieldValue( k ) ) / h1 ;

                    // zero where the secants differ in sign or one vanishes,
                    // else the weighted harmonic mean
                    if ( d0 * d1 > 0.0 )
                    {
                        const real w0 = 2.0 * h1 + h0 ;
                        const real w1 = h1 + 2.0 * h0 ;
                        mSelfFieldSlope( k ) = ( w0 + w1 ) / ( w0 / d0 + w1 / d1 ) ;
                    }
                }
                mSelfFieldSlope( 0 ) = ( mSelfFieldValue( 1 ) - mSelfFieldValue( 0 ) )
                                     / ( mSelfFieldT( 1 ) - mSelfFieldT( 0 ) ) ;
                mSelfFieldSlope( tN - 1 ) = ( mSelfFieldValue( tN - 1 ) - mSelfFieldValue( tN - 2 ) )
                                          / ( mSelfFieldT( tN - 1 ) - mSelfFieldT( tN - 2 ) ) ;
            }

            //! the measured anchor and its dT: the C1 Hermite of the knots,
            //! continued linearly from the end knots outside them
            real
            self_field( const real T, real & aDT ) const
            {
                const index_t tN = mSelfFieldT.length() ;

                if ( T <= mSelfFieldT( 0 ) )
                {
                    aDT = mSelfFieldSlope( 0 ) ;
                    return mSelfFieldValue( 0 ) + aDT * ( T - mSelfFieldT( 0 ) ) ;
                }
                if ( T >= mSelfFieldT( tN - 1 ) )
                {
                    aDT = mSelfFieldSlope( tN - 1 ) ;
                    return mSelfFieldValue( tN - 1 ) + aDT * ( T - mSelfFieldT( tN - 1 ) ) ;
                }

                index_t k = 1 ;
                while ( k < tN - 1 && mSelfFieldT( k ) < T ) ++k ;

                const real h = mSelfFieldT( k ) - mSelfFieldT( k - 1 ) ;
                const real s = ( T - mSelfFieldT( k - 1 ) ) / h ;
                const real y0 = mSelfFieldValue( k - 1 ) ;
                const real y1 = mSelfFieldValue( k ) ;
                const real m0 = mSelfFieldSlope( k - 1 ) * h ;
                const real m1 = mSelfFieldSlope( k ) * h ;

                aDT = ( ( 6. * s * s - 6. * s ) * y0 + ( 3. * s * s - 4. * s + 1. ) * m0
                      + ( -6. * s * s + 6. * s ) * y1 + ( 3. * s * s - 2. * s ) * m1 ) / h ;

                return ( 2. * s * s * s - 3. * s * s + 1. ) * y0
                     + ( s * s * s - 2. * s * s + s ) * m0
                     + ( -2. * s * s * s + 3. * s * s ) * y1
                     + ( s * s * s - s * s ) * m1 ;
            }

            //! the table on the face log10 B = log10 Bmin; d/dT and d/dangle
            //! terms only when asked for
            void
            edge( const real angle, const real T,
                  const bool aWithT, const bool aWithAngle, Edge & aEdge ) const
            {
                const real tT = std::clamp( T, mTmin, mTmax ) ;
                const real tU = mDatabase->min( 1 ) ;
                const real tA = this->wrap_angle( angle ) ;
                const real tLn10 = 1.0 / mInvLn10 ;

                const real f_u = mDatabase->evaluate_derivy( tT, tU, tA ) ;

                aEdge.V = std::pow( 10., mDatabase->evaluate( tT, tU, tA ) ) ;
                aEdge.S = aEdge.V * f_u / mBmin ;

                if ( aWithT && T > mTmin && T < mTmax )
                {
                    aEdge.V_T = aEdge.V * tLn10 * mDatabase->evaluate_derivx( tT, tU, tA ) ;
                    aEdge.S_T = ( aEdge.V_T * f_u
                                + aEdge.V * mDatabase->evaluate_derivxy( tT, tU, tA ) ) / mBmin ;
                }
                if ( aWithAngle )
                {
                    aEdge.V_A = aEdge.V * tLn10 * mDatabase->evaluate_derivz( tT, tU, tA ) ;
                    aEdge.S_A = ( aEdge.V_A * f_u
                                + aEdge.V * mDatabase->evaluate_derivyz( tT, tU, tA ) ) / mBmin ;
                }
            }

            //! true if B lies on the continuation of this label
            bool
            continued( const real normB ) const
            {
                return normB < mBmin && ( mIsN || mHaveSelfField ) ;
            }

            //! the jc anchor J = smax( a, b ) and its derivative along one of
            //! T or angle; aA_q, aV_q, aS_q are that derivative of a, V, S
            real
            anchor( const real a, const real V, const real S,
                    const real aA_q, const real aV_q, const real aS_q,
                    real & aJ_q ) const
            {
                const real b   = V - S * mBmin / 3. ;
                const real b_q = aV_q - aS_q * mBmin / 3. ;
                // the width scales with b, which is at least V > 0: a width
                // built on |a| would kink where the continued curve crosses 0
                const real e   = mSmoothWidth * b ;
                const real e_q = mSmoothWidth * b_q ;
                const real w   = std::sqrt( ( a - b ) * ( a - b ) + e * e ) ;

                aJ_q = 0.5 * ( aA_q + b_q + ( ( a - b ) * ( aA_q - b_q ) + e * e_q ) / w ) ;
                return 0.5 * ( a + b + w ) ;
            }

            //! the jc Hermite on [ 0, Bmin ]: value, and the slope at B = 0
            //! with its derivative along T or angle
            real
            jc_join( const real x, const real J, const real V, const real S,
                     const real aJ_q, const real aV_q, const real aS_q,
                     real & aValue_q, real & aDB ) const
            {
                const real tDelta = ( V - J ) / mBmin ;
                const real tBeta  = S / tDelta ;
                const bool tExact = tBeta <= 2.0 ;

                const real s0   = tExact ? 2. * tDelta - S : 0.0 ;
                const real s0_q = tExact ? 2. * ( aV_q - aJ_q ) / mBmin - aS_q : 0.0 ;

                const real h00 = 2. * x * x * x - 3. * x * x + 1. ;
                const real h10 = x * x * x - 2. * x * x + x ;
                const real h01 = -2. * x * x * x + 3. * x * x ;
                const real h11 = x * x * x - x * x ;

                aValue_q = h00 * aJ_q + h10 * mBmin * s0_q + h01 * aV_q + h11 * mBmin * aS_q ;

                aDB = ( ( 6. * x * x - 6. * x ) * J
                      + ( 3. * x * x - 4. * x + 1. ) * mBmin * s0
                      + ( -6. * x * x + 6. * x ) * V
                      + ( 3. * x * x - 2. * x ) * mBmin * S ) / mBmin ;

                return h00 * J + h10 * mBmin * s0 + h01 * V + h11 * mBmin * S ;
            }

            //! value below Bmin and one derivative of it: aMode 0 d/dB,
            //! 1 d/dT, 2 d/dangle
            real
            continuation( const real normB, const real angle, const real T,
                          const uint aMode, real & aDerivative ) const
            {
                Edge tEdge ;
                this->edge( angle, T, aMode == 1, aMode == 2, tEdge ) ;

                const real tV_q = aMode == 1 ? tEdge.V_T : tEdge.V_A ;
                const real tS_q = aMode == 1 ? tEdge.S_T : tEdge.S_A ;

                if ( mIsN )
                {
                    const real tShape = ( normB * normB - mBmin * mBmin ) / ( 2. * mBmin ) ;
                    aDerivative = aMode == 0 ? tEdge.S * normB / mBmin
                                             : tV_q + tS_q * tShape ;
                    return tEdge.V + tEdge.S * tShape ;
                }

                const bool tInside = T > mTmin && T < mTmax ;
                real a_T ;
                const real a = this->self_field( std::clamp( T, mTmin, mTmax ), a_T ) ;
                const real tA_q = ( aMode == 1 && tInside ) ? a_T : 0.0 ;

                real tJ_q ;
                const real J = this->anchor( a, tEdge.V, tEdge.S, tA_q, tV_q, tS_q, tJ_q ) ;

                real tValue_q, tDB ;
                const real tValue = this->jc_join( normB / mBmin, J, tEdge.V, tEdge.S,
                                                   tJ_q, tV_q, tS_q, tValue_q, tDB ) ;

                aDerivative = aMode == 0 ? tDB : tValue_q ;
                return tValue ;
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
            //! true where wrap_angle holds the angle on a window end, so the
            //! value does not change with it
            bool
            angle_is_snapped( const real angle ) const
            {
                return ( angle >= mAngleMin - 1e-6 && angle < mAngleMin )
                    || ( angle >  mAngleMax && angle <= mAngleMax + 1e-6 ) ;
            }

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
             * @return Jc [A/m²] or n [-]
             *
             * T is clamped to [ Tmin, Tmax ] and B above Bmax to Bmax; the
             * angle is wrapped into the table window. Below Bmin the value is
             * the continuation described above ( a jc table without a measured
             * self-field curve keeps the clamp ).
             */
            real
            eval( const real normB, const real angle, const real T ) const override
            {
                if ( this->continued( normB ) )
                {
                    real tUnused ;
                    return this->continuation( normB, angle, T, 0, tUnused ) ;
                }

                return std::pow( 10.,
                    mDatabase->evaluate(
                        std::clamp( T, mTmin, mTmax ),
                        std::log( std::clamp( normB, mBmin, mBmax ) ) * mInvLn10,
                        this->wrap_angle( angle ) ) );
            }

            /**
             * @brief d(value)/d|B|, the exact derivative of eval
             *
             * The table stores f = log10( value ) over u = log10 B, so
             * d(value)/dB = value · (∂f/∂u) / B — ln10 cancels. Zero where
             * eval clamps B ( above Bmax, and below Bmin for a jc table
             * without a measured self-field curve ); at Bmin itself the
             * continuation and the table share the slope.
             */
            real
            deval_dB( const real normB, const real angle, const real T ) const override
            {
                if ( this->continued( normB ) )
                {
                    real tDB ;
                    this->continuation( normB, angle, T, 0, tDB ) ;
                    return tDB ;
                }

                if ( normB >= mBmax ) return 0.0 ;
                if ( normB <= mBmin && ! ( mIsN || mHaveSelfField ) ) return 0.0 ;

                return this->eval( normB, angle, T )
                    * mDatabase->evaluate_derivy(
                          std::clamp( T, mTmin, mTmax ),
                          std::log( normB ) * mInvLn10,
                          this->wrap_angle( angle ) )
                    / normB ;
            }

            /**
             * @brief d(value)/dθ, the exact derivative of eval
             *
             * θ is stored linearly, so d(value)/dθ = value · ln10 · ∂f/∂θ; the
             * ±π wrap is a shift ( dθ_wrapped/dθ = 1 ), and the snap onto a
             * window end holds the value ( zero ).
             */
            real
            deval_dbeta( const real normB, const real angle, const real T ) const override
            {
                if ( this->angle_is_snapped( angle ) ) return 0.0 ;

                if ( this->continued( normB ) )
                {
                    real tDA ;
                    this->continuation( normB, angle, T, 2, tDA ) ;
                    return tDA ;
                }

                return this->eval( normB, angle, T )
                    * mDatabase->evaluate_derivz(
                          std::clamp( T, mTmin, mTmax ),
                          std::log( std::clamp( normB, mBmin, mBmax ) ) * mInvLn10,
                          this->wrap_angle( angle ) )
                    / mInvLn10 ;
            }

            /**
             * @brief d(value)/dT, the exact derivative of eval
             *
             * T is stored linearly: d(value)/dT = value · ln10 · ∂f/∂T. Outside
             * the temperature window the value is held, so the tangent is zero.
             */
            real
            deval_dT( const real normB, const real angle, const real T ) const override
            {
                if ( T <= mTmin || T >= mTmax ) return 0.0 ;

                if ( this->continued( normB ) )
                {
                    real tDT ;
                    this->continuation( normB, angle, T, 1, tDT ) ;
                    return tDT ;
                }

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