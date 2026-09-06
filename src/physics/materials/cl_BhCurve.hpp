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

#ifndef CL_BHCURVE_HPP
#define CL_BHCURVE_HPP
#include "typedefs.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief B-H curve for ferromagnetic materials
         *
         * Represents the magnetization curve B(H) and its inverse H(B) for
         * ferromagnetic materials using cubic spline interpolation of
         * tabulated data.
         *
         * The curve provides:
         * - ν(B) = H/B: Reluctivity [A/(T·m)]
         * - μ(H) = B/H: Permeability [T·m/A]
         * - dμ/dH: Derivative of permeability for Newton iterations
         * - B(H): Magnetic flux density from field intensity
         * - H(B): Magnetic field intensity from flux density
         *
         * For fields above saturation (B > Bsat or H > Hsat), the material
         * behaves linearly with an offset due to saturation magnetization:
         * - B = μ₀·( H + Msat )  for H > Hsat
         * - H = ν₀·B − Msat      for B > Bsat  ( Msat = ν₀·Bsat − Hsat, in A/m )
         *
         * Data file format (HDF5 database):
         * - aPath points to an HDF5 file that holds one group per material.
         * - aLabel selects the group (e.g. "RoxieIron", "SAE1010").
         * - Each group stores two serialized cubic splines, "bnur" ( ν(B) ) and
         *   "hnur" ( 1/μ(H) ), plus the scalars "bsat" and "hsat".
         *
         * Usage (BhCurve is abstract; construct the concrete BhSplineCurve):
         * @code
         * BhCurve* bh = new BhSplineCurve( "bhdata.hdf5", "RoxieIron" );
         * material->load_bh_curve(bh);  // Material takes ownership
         *
         * // Later, in FEM assembly:
         * real H = material->H(B);  // field intensity at |B| through the curve
         *                           // ( reluctivity itself is bh->nu(B); Material::nu() is Poisson's ratio )
         * @endcode
         *
         * @ingroup grp_physics_materials
         * @see @ref physics_materials_materials_usage_guide
         */
        class BhCurve
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            const proc_t mCommRank ;   //!< MPI rank
            const string mPath ;      //!< path where data file is located
            const string mLabel ;     //!< Material label for error messages



//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * @brief Constructor - records the database path and group label
             * @param aPath  Path to the HDF5 database file (e.g. "bhdata.hdf5")
             * @param aLabel Group name to load (e.g. "RoxieIron")
             *
             * The concrete BhSplineCurve opens the HDF5 file, selects the group
             * named aLabel and reads the pre-computed cubic splines "bnur"
             * ( ν(B) ) and "hnur" ( 1/μ(H) ) together with the saturation
             * scalars "bsat" and "hsat"; Msat is derived from these.
             */
            BhCurve( const string & aPath, const string & aLabel );

            /**
             * @brief Virtual destructor ( the concrete curve owns and deletes its splines )
             */
            virtual ~BhCurve() = default ;

            /**
             * @brief Reluctivity ν(B) = H/B
             * @param B Magnetic flux density magnitude [T]
             * @return Reluctivity ν [A/(T·m)]
             *
             * For B < Bsat: Uses cubic spline interpolation
             * For B ≥ Bsat: Returns ν₀ - Msat/B (linear + offset)
             */
            virtual real
            nu( const real B ) const;

            /**
             * @brief Permeability μ(H) = B/H
             * @param H Magnetic field intensity magnitude [A/m]
             * @return Permeability μ [T·m/A]
             *
             * For H < Hsat: Uses cubic spline interpolation
             * For H ≥ Hsat: Returns μ₀·(1 + Msat/H)
             */
            virtual real
            mu( const real H ) const ;

            /**
             * @brief Permeability and its derivative with respect to H
             * @param H Magnetic field intensity magnitude [A/m]
             * @param[out] mu Permeability μ(H) [T·m/A]
             * @param[out] dmudH Derivative dμ/dH [T·m/A²]
             *
             * Computes both μ and dμ/dH efficiently for Newton-Raphson
             * iterations in nonlinear magnetic FEM solvers.
             *
             * For H < Hsat: Uses spline evaluation and derivative
             * For H ≥ Hsat: Analytical formulas for linear regime
             */
            virtual void
            dmudH( const real H, real & mu, real & dmudH ) const ;

            /**
             * @brief Magnetic flux density from field intensity
             * @param H Magnetic field intensity magnitude [A/m]
             * @return Magnetic flux density B [T]
             *
             * Convenience function: B = μ(H) · H
             */
            real
            B( const real H ) const;

            /**
             * @brief Magnetic field intensity from flux density
             * @param B Magnetic flux density magnitude [T]
             * @return Magnetic field intensity H [A/m]
             *
             * Convenience function: H = ν(B) · B
             */
            real
            H( const real B ) const ;

        };

//------------------------------------------------------------------------------

        inline real
        BhCurve::B( const real H ) const
        {
            return this->mu( H ) * H ;
        }

        inline real
        BhCurve::H( const real B ) const
        {
            return this->nu( B ) * B ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //CL_BHCURVE_HPP
