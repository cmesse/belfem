/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * ASCII Art (c) 1996-2001   Joan G. Stark, used with permission.
 */

#include <iostream>
#include <string>
#ifdef OMP
#include <omp.h>
#include <cstdlib>   // getenv, to tell an explicit OMP_NUM_THREADS from a default
#include <fstream>   // /proc/cpuinfo, for the physical core count
#include <set>       // unique ( package, core ) pairs -- once-per-run setup only
#include <utility>
#if defined( __linux__ )
#include <sched.h>   // sched_getaffinity: what this process may actually run on
#elif defined( __APPLE__ )
#include <sys/sysctl.h>
#endif
#endif

#include "commtools.hpp"

#include "assert.hpp"
#include "fn_create_graph_from_matrix.hpp"
#include "cl_SolverWrapper.hpp"
#include "cl_Cell.hpp"
#include "fn_Graph_symrcm.hpp"
#include "fn_compute_permutation.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        Wrapper::Wrapper( const string & aLabel, const bool aUsesMPI ) :
            mCommRank(   gComm.rank() ),
            mCommSize( gComm.size() ),
            mUsesMPI( aUsesMPI ),
            mLabel( aLabel )
        {

        }
//------------------------------------------------------------------------------

        Wrapper::~Wrapper()
        {
            Wrapper::free();
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
                    const int_t aNumRhsColumns )
        {
            BELFEM_ERROR( false, "initialize() is not implemented for %s",
                         mLabel.c_str() );
        }

//------------------------------------------------------------------------------

        void
        Wrapper::free()
        {
            CommunicationObject::free();

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
        Wrapper::mat2vec(
            const Matrix< real > & aM,
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

        real
        Wrapper::get_cond1() const
        {
            BELFEM_ERROR( false,
                          "get_cond1() is not implemented for solver %s",
                          mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        real
        Wrapper::get_cond2() const
        {
            BELFEM_ERROR( false,
                          "get_cond2() is not implemented for solver %s",
                          mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        real
        Wrapper::get_forward_error() const
        {
            BELFEM_ERROR( false,
                          "get_forward_error() is not implemented for solver %s",
                          mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        real
        Wrapper::get_backward_error() const
        {
            BELFEM_ERROR( false,
                          "get_backward_error() is not implemented for solver %s",
                          mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        real
        Wrapper::get_omega2() const
        {
            BELFEM_ERROR( false,
                          "get_omega2() is not implemented for solver %s",
                          mLabel.c_str() );

            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        bool
        Wrapper::supports_factorization_reuse() const
        {
            // opt-in per wrapper. Saying "no" is always safe: the caller then
            // simply re-factorizes, which is what every path did before the
            // frozen scope existed
            return false ;
        }

//------------------------------------------------------------------------------

        void
        Wrapper::freeze_factorization( const SpMatrix & aMatrix )
        {
            // deliberately NOT a silent no-op. A caller that freezes and then
            // solves believes it is reusing a factorization ; if the wrapper
            // quietly ignored that, the caller would still get right answers
            // but at many times the cost, and nothing would ever say so
            BELFEM_ERROR( false,
                          "freeze_factorization() is not implemented for solver %s - "
                          "ask supports_factorization_reuse() before arming the scope",
                          mLabel.c_str() );
        }

//------------------------------------------------------------------------------

        void
        Wrapper::unfreeze_factorization()
        {
            // safe when nothing is frozen, so a scope guard may call it
            // unconditionally on the way out -- including on the path where
            // freeze_factorization() itself threw
        }

//------------------------------------------------------------------------------

        bool
        Wrapper::factorization_is_frozen() const
        {
            return false ;
        }


//------------------------------------------------------------------------------

#ifdef OMP
        namespace
        {
            /**
             * physical cores this process may actually run on, or 0 when that
             * cannot be established.
             *
             * The affinity mask is the right basis rather than the machine's core
             * count, because every mechanism that narrows an allocation -- a SLURM
             * cgroup, a cpuset, taskset, an MPI binding policy -- narrows the mask.
             * A turtle that counted a 24-core node inside an 8-core allocation
             * would under-warn.
             *
             * Deliberately NOT reading SLURM_CPUS_ON_NODE / SLURM_CPUS_PER_TASK:
             * sites differ on whether those count logical or physical CPUs under
             * SMT, so they cannot be compared against a physical-core budget
             * without knowing the site's convention.
             *
             * Setup-path code, called once per solver: the std::set and the
             * ifstream are the one-off-allocation case doc/coding_philosophy.md
             * allows, not a hot-path pattern.
             */
            unsigned int
            physical_cores_available( const bool aRespectAffinity )
            {
#if defined( __linux__ )
                cpu_set_t tMask ;
                CPU_ZERO( &tMask );
                if ( aRespectAffinity )
                {
                    if ( sched_getaffinity( 0, sizeof( tMask ), &tMask ) != 0 )
                    {
                        return 0 ;
                    }
                }
                else
                {
                    // count the whole machine: every cpu is "in" the mask
                    for ( int k = 0; k < CPU_SETSIZE; ++k ) CPU_SET( k, &tMask );
                }

                std::ifstream tInfo( "/proc/cpuinfo" );
                if ( ! tInfo.is_open() )
                {
                    return 0 ;
                }

                // a physical core is a unique ( package, core ) pair; the two
                // hyperthreads of one core report the same pair
                std::set< std::pair< int, int > > tCores ;

                int tCpu  = -1 ;
                int tPkg  = -1 ;
                int tCore = -1 ;

                std::string tLine ;
                while ( std::getline( tInfo, tLine ) )
                {
                    const std::size_t tColon = tLine.find( ':' );
                    if ( tColon == std::string::npos ) continue ;

                    const std::string tKey = tLine.substr( 0, tColon );
                    const int tValue = std::atoi( tLine.c_str() + tColon + 1 );

                    if ( tKey.rfind( "processor", 0 ) == 0 )
                    {
                        // a new block begins: the previous one is complete
                        if ( tCpu >= 0 && tPkg >= 0 && tCore >= 0
                             && CPU_ISSET( tCpu, &tMask ) )
                        {
                            tCores.insert( { tPkg, tCore } );
                        }
                        tCpu  = tValue ;
                        tPkg  = -1 ;
                        tCore = -1 ;
                    }
                    else if ( tKey.rfind( "physical id", 0 ) == 0 ) tPkg  = tValue ;
                    else if ( tKey.rfind( "core id",     0 ) == 0 ) tCore = tValue ;
                }

                // the last block has no successor to flush it
                if ( tCpu >= 0 && tPkg >= 0 && tCore >= 0 && CPU_ISSET( tCpu, &tMask ) )
                {
                    tCores.insert( { tPkg, tCore } );
                }

                return ( unsigned int ) tCores.size() ;

#elif defined( __APPLE__ )
                int    tCores = 0 ;
                size_t tSize  = sizeof( tCores );
                if ( sysctlbyname( "hw.physicalcpu", &tCores, &tSize, nullptr, 0 ) != 0 )
                {
                    return 0 ;
                }
                return tCores > 0 ? ( unsigned int ) tCores : 0 ;
#else
                // unknown platform: no budget, therefore no warning and no advice
                return 0 ;
#endif
            }
        }
#endif

//------------------------------------------------------------------------------

        void
        Wrapper::hatch_turtle()
        {
// NOT a mistake: this is OMP, not BELFEM_OMP.
// BELFEM_OMP gates BELFEM's OWN Fortran kernels and is OFF by default.
// hatch_turtle has nothing to do with those kernels -- it reports on the THIRD-PARTY
// thread budget ( STRUMPACK, MKL, threaded BLAS ), which is exactly the behaviour we
// cannot influence and therefore most need to warn about. Moving it to BELFEM_OMP
// would silently delete the oversubscription warning on every default build.
#ifdef OMP
            // The budget is PHYSICAL cores, not logical CPUs. Until 2026-08-21 this
            // used std::thread::hardware_concurrency() -- the logical count -- and
            // then RECOMMENDED that budget. On the 10-core / 20-thread reference
            // workstation with ten ranks that recommends OMP_NUM_THREADS = 2, i.e.
            // 20 threads on 10 cores: exactly the configuration
            // examples/scripts/Allrun refuses to launch and doc/parallel_execution.md
            // forbids. Two components of this repository disagreed about what a core
            // is; this one was wrong.
            const unsigned int tMine       = physical_cores_available( true ) ;
            const unsigned int tMaxThreads = omp_get_max_threads() ;
            const unsigned int tRanks      = gComm.node_size() ;
            const unsigned int tMachine    = tMine > 0 ? physical_cores_available( false ) : 0 ;

            // An unknown budget must never become a RECOMMENDATION -- that is the
            // defect being fixed, so there is no logical-CPU fallback. It is not a
            // reason for silence, though: an unset OMP_NUM_THREADS with several
            // ranks on a node IS the affinity-mask default, and diagnosing it needs
            // no core count at all ( Grok, 2026-08-21 ). On a platform this parse
            // cannot budget -- ARM, some VMs, a mask wider than CPU_SETSIZE -- that
            // would otherwise silence the one case the turtle exists to catch.
            const bool tKnown = tMine > 0 ;

            // Are the ranks SHARING one mask, or does each hold its own slice?
            // Without a node-local communicator this rank cannot see its siblings'
            // masks, but the two cases separate arithmetically: disjoint slices
            // must fit in the machine, a shared mask cannot.
            //
            //   unbound, 10 ranks x 10 cores  = 100 > 10  -> shared  -> product test
            //   bound,    4 ranks x  2 cores  =   8 <= 10 -> slices  -> per-rank test
            //
            // Getting this wrong in the SHARED direction is the expensive mistake:
            // warning on a correctly bound run teaches the user to ignore the box.
            const bool tSliced = tKnown && tMachine > 0 && tMine * tRanks <= tMachine ;

            // Test the PRODUCT when the mask is shared. The old form divided first
            // and returned early when the quotient was zero, so 16 ranks x 1 thread
            // on 10 cores -- genuinely oversubscribed -- reported nothing at all.
            const unsigned int tCores     = tMine ;
            const unsigned int tRequested = tSliced ? tMaxThreads : tRanks * tMaxThreads ;

            // Christian's question, answerable only here: a user who never set the
            // variable is looking at a runtime default ( the affinity mask popcount,
            // hence every hyperthread ), not at a choice they made.
            const bool tExplicit = std::getenv( "OMP_NUM_THREADS" ) != nullptr ;

            // two reasons to speak, and only one of them can name a number
            const bool tOversubscribed = tKnown && tRequested > tCores ;
            const bool tBlindDefault   = ! tKnown && ! tExplicit && tRanks > 1 ;

            if ( ! tOversubscribed && ! tBlindDefault ) return ;

            // only the root proc prints
            if ( this->rank() != 0 ) return ;

            // Threads per rank that WOULD fit. Never a rank count: the even-rank
            // policy that STRUMPACK's proportional mapping needs belongs to
            // examples/scripts/Allrun, which owns it.
            const unsigned int tFits = tSliced ? tCores
                                             : ( tCores / tRanks > 0 ? tCores / tRanks : 1 ) ;

            Cell< string > tMsg ;

            tMsg.push( "" );
            tMsg.push( "    ╔════════════════════════════════════════════════════════════════════════╗" ) ;
            tMsg.push( "    ║                                                                        ║" ) ;
            tMsg.push( "    ║                      __     WARNING:                                   ║" ) ;
            tMsg.push( "    ║           .,-;-;-,. /'_\\                                               ║" ) ;
            tMsg.push( "    ║         _/_/_/_|_\\_\\) /     You are oversubscribing threads or cores.  ║" ) ;
            tMsg.push( "    ║       '-<_><_><_><_>=/\\     This costs time and memory; it does        ║" ) ;
            tMsg.push( "    ║  jgs    `/_/====/_/-'\\_\\    not speed the calculation.                 ║" ) ;
            tMsg.push( "    ║          \"\"     \"\"    \"\"                                               ║" ) ;
            tMsg.push( "    ╠════════════════════════════════════════════════════════════════════════╣" ) ;
            tMsg.push( "    ║                                                                        ║" ) ;
            tMsg.push( tExplicit ?
                sprint("    ║   OMP_NUM_THREADS                        : %10lu ( set )          ║", ( long unsigned int ) tMaxThreads ) :
                sprint("    ║   OMP_NUM_THREADS                        : %10lu ( UNSET )        ║", ( long unsigned int ) tMaxThreads ) );
            if ( ! tExplicit )
            {
                tMsg.push( "    ║     the runtime defaulted to the affinity mask, which counts every     ║" );
                tMsg.push( "    ║     hyperthread. It is not a recommendation.                           ║" );
            }
            tMsg.push( sprint("    ║   MPI processes on this node             : %10u                  ║", ( unsigned int ) tRanks ) );
            if ( ! tKnown )
            {
                tMsg.push( "    ║   PHYSICAL cores available here          :    UNKNOWN                  ║" );
                tMsg.push( "    ║     this platform's core topology could not be read, so no value       ║" );
                tMsg.push( "    ║     is recommended here -- see doc/parallel_execution.md               ║" );
            }
            if ( tKnown ) tMsg.push( sprint( tSliced ?
                "    ║   threads requested by this rank         : %10lu                  ║" :
                "    ║   ranks x threads requested              : %10lu                  ║",
                ( long unsigned int ) tRequested ) );
            if ( tKnown ) tMsg.push( sprint("    ║   PHYSICAL cores available here          : %10lu                  ║", ( long unsigned int ) tCores ) );
            tMsg.push( "    ║                                                                        ║");
            tMsg.push( "    ╠════════════════════════════════════════════════════════════════════════╣" ) ;
            tMsg.push( "    ║                                                                        ║");
            tMsg.push( "    ║  Measured on this code ( doc/parallel_execution.md ):                  ║" );
            tMsg.push( "    ║    4 ranks x 4 threads on 10 cores ran the factorization ~15 %         ║" );
            tMsg.push( "    ║    slower than 4 x 2, and used ~13 GiB more memory.                    ║" );
            tMsg.push( "    ║                                                                        ║");
            tMsg.push( "    ║  Recommendation :                                                      ║" );
            if ( tKnown ) tMsg.push( sprint("    ║    * export OMP_NUM_THREADS=%-3lu and relaunch                           ║", ( long unsigned int ) tFits ) );
            tMsg.push( "    ║    * for binding flags see doc/parallel_execution.md                   ║" );
            tMsg.push("    ║                                                                        ║" );
            tMsg.push( "    ╚════════════════════════════════════════════════════════════════════════╝" );
            tMsg.push( "" );

            for ( const string & tLine : tMsg )
            {
                std::cerr << tLine << std::endl;
            }
#endif
        }

//------------------------------------------------------------------------------
    }
}
