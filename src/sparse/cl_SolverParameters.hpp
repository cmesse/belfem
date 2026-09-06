/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_SOLVERPARAMETERS_HPP
#define BELFEM_CL_SOLVERPARAMETERS_HPP
#include "en_SolverEnums.hpp"
#include "cl_Input_Section.hpp"
#include "commtools.hpp"

namespace belfem
{
    /**
     * @brief Configuration for a Solver.
     *
     * @ingroup grp_sparse
     * @see @ref sparse_sparse_usage_guide
     */
    class SolverParameters
    {
        const proc_t     mCommRank ;
        const SolverType mSolverType = gDefaultSolver ;

        // CSR is recommended: STRUMPACK's set_distributed_csr_matrix() is
        // its primary and most-tested MPI interface. The MPIAIJ path
        // (set_MPIAIJ_matrix) causes hangs during distributed multifrontal
        // factorization and lacks upstream C++ test coverage.
        DistributedMatrixType mDistributedMatrixType = DistributedMatrixType::CSR ;

        ReorderingMethod  mReorderingMethod   = ReorderingMethod::AUTOMATIC ;

        CompressionMethod mCompressionMethod  = CompressionMethod::AUTOMATIC ;

        // BLR compression tolerance, consumed ONLY when the deck says
        // compression scheme : blr. One number, two strictnesses:
        // STRUMPACK applies it as a RELATIVE rel_tol, MUMPS as the
        // ABSOLUTE CNTL(7) dropping parameter ( MUMPS 5.7.3 §5.19 ).
        // See src/sparse/doc/solver_memory_and_compression.md §4
        real mCompressionCutoff = 1e-8 ;

        // true once the deck ( or a programmatic caller ) stated a
        // cutoff — distinguishes a chosen trade from the inherited
        // default in the headroom check
        bool mHaveCompressionCutoff = false ;

        // MUMPS only: per-process working-memory cap in MB, handed to
        // ICNTL(23) before the first factorization. 0 = not stated; the
        // wrapper then measures the machine itself in initialize() ( see
        // MUMPS::memory_budget_mb ) and applies that number as ICNTL(23)
        // on the first out-of-workspace failure. An explicit value wins
        // over that measurement
        uint mMemoryBudget = 0 ;
        bool mHaveMemoryBudget = false ;

        // MC64 matrix matching, only used by STRUMPACK. Its permutation is
        // load-bearing for the coupled h-phi Jacobian: without it the Newton
        // stage hits exact zero pivots that replace_tiny_pivots cannot absorb
        // (observed 2026-07-05, serial hphirun, ZERO_PIVOT at Newton 17).
        // Expensive in MPI mode ( STRUMPACK gathers the full matrix on rank 0,
        // once per initialize ) -- ( matching : off ) is an opt-out for cases
        // verified to tolerate it.
        bool mUseMatrixMatching = true ;

        // STRUMPACK / serial-METIS only. METIS_NodeNDP returns the separator
        // tree; METIS_NodeND does not, so STRUMPACK has to rebuild a supernodal
        // tree from the elimination tree instead, and that reconstruction is
        // what produces very deep trees. Observed on tapestack3d 2026-08-12,
        // 832k magnetic dofs: "used METIS_NodeND (iso METIS_NodeNDP)",
        // "supernodal tree was built from etree", 56 levels — and STRUMPACK's
        // own diagnostic then warned that it "does not handle this safely,
        // which could lead to segmentation faults due to stack overflows",
        // recommending precisely this flag. The same deep tree is the likely
        // cause of the factorisation memory landing overwhelmingly on one rank
        // ( 23 GiB against 2.8 GiB on its siblings ).
        //
        // METIS_NodeNDP is UNDOCUMENTED in METIS and declared by STRUMPACK
        // itself, so it is kept switchable: ( metis nodendp : false ) restores
        // METIS_NodeND if a METIS build does not provide it. Irrelevant to
        // PARMETIS and SCOTCH, which ignore the flag.
        bool mUseMetisNodeNDP = true ;

        // Only used for PETSc, ignored otherwise.
        //
        // ASM ( additive Schwarz, PETSc default sub-solver ILU(0) ) rather than
        // the former GAMG-in-parallel / JACOBI-in-serial pair, for two reasons.
        //
        // 1. The consumer is not the operator GAMG is built for. Smoothed
        //    aggregation assumes a symmetric elliptic system; that describes
        //    the PICARD thermal operator ( M/dt + K, SPD ) but not the NEWTON
        //    tangent, which carries the non-symmetric mixed-operator term
        //    ( B'*grad T ) (x) ( dlambda/dT * N ) from mt_thermal_h. Measured
        //    on tapestack3d 2026-08-12, 337k thermal dofs, rtol 1e-8:
        //    GAMG needed ~74 Krylov iterations on the two Picard iterates and
        //    ~700 on every Newton iterate of the same step -- a 9x penalty
        //    that appears exactly when the algorithm flips. ASM is slightly
        //    worse on the Picard operator ( ~107 ) and does not degrade on the
        //    Newton one, which is the trade this default takes.
        //
        // 2. The old pair changed the preconditioner WITH THE RANK COUNT, so a
        //    serial and a parallel run of the same deck solved with different
        //    numerics. Same defect class as the thermal Picard freeze, where
        //    partition-order roundoff selected the nonlinear algorithm. ASM with one block
        //    degenerates to ILU, so the family is now continuous in comm_size.
        //
        // A deck that knows its operator is symmetric can still ask for GAMG
        // through ( preconditioner : gamg ) -- and should, since GAMG wins on
        // the Picard operator.
        Preconditioner    mPreconditioner     = Preconditioner::ASM ;

        // only used for PETSc and STRUMPACK. AUTO resolves to PREONLY behind a
        // direct LU preconditioner and GMRES otherwise ( cl_SolverPETSC.cpp ),
        // which is the correct choice for the non-symmetric Newton tangent
        // above -- a symmetric-only method such as CG must not become the
        // default while that term is in the operator.
        KrylovMethod      mKrylovMethod       = KrylovMethod::AUTO ;

        // Relative tolerance: PETSc KSP rtol AND STRUMPACK outer-GMRES
        // rel_tol, always applied ( STRUMPACK since 2026-08-18; before
        // that only when stated, and an unstated deck inherited the
        // library's 1e-6 — the tapestack3d A/B showed that loose exit
        // test WAS the printed nonlinear residual, see strumpacktools ).
        //
        // 1e-10 rather than 1e-8 because a transient solved for an ABSOLUTE
        // field accepts an absolute error of ~rtol * |field| per step. A
        // thermal run in kelvin therefore drifts by ~rtol * T every step, and
        // in an adiabatic or weakly forced problem nothing pulls it back:
        // measured on tapestack3d 2026-08-13, 337k thermal dofs, a uniform
        // 77 K start drifted 8.6e-5 K in 23 steps at rtol 1e-8, and 2.6e-7 K
        // at rtol 1e-10. The failure is silent -- no residual, no message --
        // so the default carries the safety margin rather than the deck.
        // ( The old reassurance "STRUMPACK pays nothing, the factorization
        // already delivers ~2e-16" was refuted by the A/B: on an
        // ill-conditioned matrix the raw factor gives ~4 digits and the
        // OUTER GMRES is what reaches deep residuals. )
        //
        // This is a mitigation, not the cure. The cure is to solve in
        // increment form, where the criterion becomes scale-free.
        real mRelativeTolerance = 1e-10;

        // Set when the deck or set_relative_tolerance() stated a value.
        // NOT a consumer gate: since 2026-08-18 STRUMPACK always applies
        // mRelativeTolerance ( same as PETSc ). Do not re-introduce the
        // have_relative_tolerance() guard in strumpacktools — that was
        // the REFINE-era split, retired after the tapestack3d
        // 1e-8 vs 1e-10 A/B; stall protection is now PREC_GMRES +
        // maxit 50 + always-set abs_tol. The flag still rides
        // synchronize() so the positional payload keeps its slot ( current
        // width 13, see the history note in synchronize() ).
        bool mHaveRelativeTolerance = false ;

        // Absolute tolerance of the ITERATIVE part of a linear solve.
        //
        // STRUMPACK: exit test of the outer GMRES next to rel_tol. The
        // library default is 1e-10 ( StrumpackOptions.hpp ), and BELFEM
        // shipped with it unset for years — on a SMALL right-hand side
        // ( transient startup, ||b|| ~ 1e-2 ) that floor sits exactly at
        // a 1e-11 RELATIVE nonlinear target and Newton converges only by
        // GMRES overshoot lottery ( three-voice jury 2026-08-17 ).
        // 1e-14 restores >= 3 decades of headroom while keeping GMRES:
        // the factorization-preconditioned iteration reaches it in a few
        // extra Arnoldi steps, and the maxit cap bounds the cost.
        //
        // PETSc: passed as KSP atol when stated in the deck; PETSC_DEFAULT
        // otherwise ( the historical behavior ).
        real mAbsoluteTolerance = 1e-14 ;

        // true once the deck ( or a programmatic caller ) stated an
        // absolute tolerance. PETSc applies the value only then;
        // STRUMPACK always applies mAbsoluteTolerance, because its
        // library default is the overshoot-lottery defect above.
        bool mHaveAbsoluteTolerance = false ;

        bool mUseInitialGuess = false ;

        // deck-stated Krylov iteration budget for ITERATIVE solvers
        // ( "max iterations" in a linear section — distinct from the
        // nonlinear key of the same name, which the controller parses ).
        // 0 = not stated, keep the library default ( PETSc: 10000 ).
        // Consumed by PETSc only; STRUMPACK's refinement cap is a fixed
        // safety net against the refinement stall, deliberately not deck-tunable
        uint mMaxNumIterations = 0 ;

    public:

        SolverParameters( const SolverType aType ) ;

        SolverParameters( const input::Section * aInput );

        SolverParameters( const SolverParameters & aOther );

        SolverParameters(SolverParameters&& aOther) = delete ;
        SolverParameters& operator=(SolverParameters&& aOther) = delete ;

        ~SolverParameters() = default ;

        void
        set_distributed_matrix_type( const DistributedMatrixType aType );

        void
        set_preconditioner( const Preconditioner aPreconditioner );

        void
        set_krylov_method( const KrylovMethod aKrylovMethod );

        void
        set_reordering_method( const ReorderingMethod aReorderingMethod );

        void
        set_compression_method( const CompressionMethod aCompressionMethod );

        void
        set_relative_tolerance( const real aEpsilon );

        void
        set_absolute_tolerance( const real aEpsilon );

        void
        set_use_initial_guess( const bool aUse );

        void
        set_matrix_matching( const bool aUse );

        SolverType
        type() const ;

        DistributedMatrixType
        distributed_matrix_type() const ;

        Preconditioner
        preconditioner() const;

        KrylovMethod
        krylov_method() const;

        ReorderingMethod
        reordering_method() const;

        CompressionMethod
        compression_method() const;

        real
        relative_tolerance() const ;

        real
        absolute_tolerance() const ;

        //! true if the tolerance was stated explicitly rather than
        //! inherited from the class default
        bool
        have_relative_tolerance() const ;

        bool
        have_absolute_tolerance() const ;

        //! deck-stated Krylov iteration budget; 0 = library default
        uint
        max_iterations() const ;

        //! validates > 0 and finite — the single check for the deck
        //! parser and programmatic callers alike
        void
        set_compression_cutoff( const real aCutoff );

        //! true if the cutoff was stated rather than inherited
        bool
        have_compression_cutoff() const ;

        real
        compression_cutoff() const ;

        //! validates > 0 — the single check for the deck parser and
        //! programmatic callers alike. Megabytes per process
        void
        set_memory_budget( const uint aMegaBytes );

        //! true if the deck ( or a caller ) stated a budget
        bool
        have_memory_budget() const ;

        //! MB per process; 0 when not stated
        uint
        memory_budget() const ;


        bool
        use_initial_guess() const ;

        bool
        use_matrix_matching() const ;

        bool
        use_metis_nodendp() const ;

        void
        synchronize();

    private:

        SolverType
        get_solver_type_from_input( const input::Section * aInput ) const ;

    };
}
#endif //BELFEM_CL_SOLVERPARAMETERS_HPP