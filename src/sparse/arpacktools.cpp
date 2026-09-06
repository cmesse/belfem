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

#include "typedefs.hpp"

#include "assert.hpp"
#include "cl_Logger.hpp"
#include "arpacktools.hpp"

namespace belfem
{
    namespace arpack
    {
//------------------------------------------------------------------------------

        void
        check_naupd( const int_t aInfo )
        {
            // the non-fatal exits: Ritz values have been computed, they are
            // just not all converged. dneupd still extracts what there is,
            // so the run continues -- but the caller must not trust the
            // result without also looking at NCONV
            switch ( aInfo )
            {
                case 0 :
                {
                    return ;
                }
                case 1 :
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dnaupd hit the restart iteration limit before converging.\n"
                        "                 Check the converged count before using the result.\n" );
                    return ;
                }
                case 3 :
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dnaupd could not apply any shifts - the Krylov subspace is too\n"
                        "                 narrow for this spectrum. Increase ncv relative to nev.\n" );
                    return ;
                }
                case -9999 :
                {
                    // the algorithm failed, the caller did not. Same tier as
                    // the two above: report and let the converged count decide
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dnaupd could not build an Arnoldi factorization.\n"
                        "                 No eigenvalue is available for this step.\n" );
                    return ;
                }
                default :
                {
                    break ;
                }
            }

            string tMessage ;

            switch ( aInfo )
            {
                case -1 :
                    tMessage = "N must be positive." ;
                break ;
                case -2 :
                    tMessage = "NEV must be positive." ;
                break ;
                case -3 :
                    tMessage = "NCV-NEV >= 2 and less than or equal to N." ;
                break ;
                case -4 :
                    tMessage = "The maximum number of Arnoldi update iterations must be greater than zero." ;
                break ;
                case -5 :
                    tMessage = "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'." ;
                break ;
                case -6 :
                    tMessage = "BMAT must be one of 'I' or 'G'." ;
                break ;
                case -7 :
                    tMessage = "Length of private work array is not sufficient." ;
                break ;
                case -8 :
                    tMessage = "Error return from LAPACK eigenvalue calculation." ;
                break ;
                case -9 :
                    tMessage = "Starting vector is zero." ;
                break ;
                case -10 :
                    tMessage = "IPARAM(7) must be 1,2,3,4." ;
                break ;
                case -11 :
                    tMessage = "IPARAM(7) = 1 and BMAT = 'G' are incompatible." ;
                break ;
                case -12 :
                    tMessage = "IPARAM(1) must be equal to 0 or 1." ;
                break ;
                case -9999 :
                    tMessage = "Could not build an Arnoldi factorization." ;
                break ;

                // 100 and up are raised by the BELFEM drivers, not by ARPACK.
                // 101-105 can only come from the distributed driver, where
                // every one of them is decided collectively so that all ranks
                // report the same code
                case 100 :
                    tMessage = "The driver received an unexpected reverse communication request." ;
                break ;
                case 101 :
                    tMessage = "The local row blocks do not sum to the global row count, "
                               "or the global row count is not the same on every rank." ;
                break ;
                case 102 :
                    tMessage = "At least one rank owns no rows. PARPACK needs a positive "
                               "local dimension on every rank of the communicator." ;
                break ;
                case 103 :
                    tMessage = "The local number of nonzeros contradicts the row pointers "
                               "( pointers( nloc + 1 ) - 1 must equal nnz )." ;
                break ;
                case 104 :
                    tMessage = "The job flag must be 0 ( smallest magnitude ) or 1 ( largest magnitude )." ;
                break ;
                case 105 :
                    tMessage = "An MPI call failed while the row map was being built." ;
                break ;
                default:
                    tMessage = "Unknown Error" ;
            }

            BELFEM_ERROR( false,
                "ARPACK-ng has thrown error %i in dnaupd : %s",
                    ( int ) aInfo,
                    tMessage.c_str() );
        }

//------------------------------------------------------------------------------

        void
        check_neupd( const int_t aInfo )
        {
            if( aInfo == 0 )
            {
                return ;
            }

            // -14 is the extraction side of a non-converged dnaupd: the
            // algorithm did not reach the requested accuracy, which is an
            // expected outcome rather than a fault. Aborting here would kill
            // the run before the caller's converged count can turn it into a
            // missing diagnostic
            if( aInfo == -14 )
            {
                message( InfoLevel::Verbose,
                    "ARPACK-ng dneupd found no eigenvalue of sufficient accuracy.\n"
                    "                 No eigenvalue is available for this step.\n" );
                return ;
            }

            string tMessage ;

            switch ( aInfo )
            {
                case 1 :
                    tMessage = "The Schur form computed by dlahqr could not be reordered by dtrsen." ;
                break ;
                case -1 :
                    tMessage = "N must be positive." ;
                break ;
                case -2 :
                    tMessage = "NEV must be positive." ;
                break ;
                case -3 :
                    tMessage = "NCV-NEV >= 2 and less than or equal to N." ;
                break ;
                case -5 :
                    tMessage = "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'." ;
                break ;
                case -6 :
                    tMessage = "BMAT must be one of 'I' or 'G'." ;
                break ;
                case -7 :
                    tMessage = "Length of private work WORKL array is not sufficient." ;
                break ;
                case -8 :
                    tMessage = "Error return from calculation of a real Schur form." ;
                break;
                case -9 :
                    tMessage = "Error return from calculation of eigenvectors." ;
                break;
                case -10 :
                    tMessage = "IPARAM(7) must be 1,2,3,4." ;
                break ;
                case -11 :
                    tMessage = "IPARAM(7) = 1 and BMAT = 'G' are incompatible." ;
                break ;
                case -12 :
                    tMessage = "HOWMNY = 'S' not yet implemented" ;
                break ;
                case -13 :
                    tMessage = "HOWMNY must be one of 'A' or 'P' if RVEC = .true." ;
                break ;
                case -14 :
                    tMessage = "DNAUPD did not find any eigenvalues to sufficient accuracy." ;
                break;
                case -15 :
                    tMessage = "DNEUPD got a different count of the number of converged Ritz values than DNAUPD got." ;
                break;
                default:
                    tMessage = "Unknown Error" ;
            }

            BELFEM_ERROR( false,
                "ARPACK-ng has thrown error %i in dneupd : %s",
                    ( int ) aInfo,
                    tMessage.c_str() );
        }

//------------------------------------------------------------------------------

        void
        check_saupd( const int_t aInfo )
        {
            // same non-fatal exits as the nonsymmetric driver, info = 3
            // included: dsaupd documents it ( "no shifts could be applied",
            // raise NCV ) and the Fortran driver extracts through it
            switch ( aInfo )
            {
                case 0 :
                {
                    return ;
                }
                case 1 :
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dsaupd hit the restart iteration limit before converging.\n"
                        "                 Check the converged count before using the result.\n" );
                    return ;
                }
                case 3 :
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dsaupd could not apply any shifts in a restart cycle.\n"
                        "                 Increase ncv relative to nev.\n" );
                    return ;
                }
                case -9999 :
                {
                    message( InfoLevel::Verbose,
                        "ARPACK-ng dsaupd could not build a Lanczos factorization.\n"
                        "                 No eigenvalue is available for this step.\n" );
                    return ;
                }
                default :
                {
                    break ;
                }
            }

            string tMessage ;

            switch ( aInfo )
            {
                case -1 :
                    tMessage = "N must be positive." ;
                break ;
                case -2 :
                    tMessage = "NEV must be positive." ;
                break ;
                case -3 :
                    tMessage = "NCV-NEV >= 2 and less than or equal to N." ;
                break ;
                case -4 :
                    tMessage = "The maximum number of Arnoldi update iterations must be greater than zero." ;
                break ;

                // the symmetric table admits three modes the nonsymmetric one
                // rejects, which is exactly why this decoder is separate
                case -5 :
                    tMessage = "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'." ;
                break ;
                case -6 :
                    tMessage = "BMAT must be one of 'I' or 'G'." ;
                break ;
                case -7 :
                    tMessage = "Length of private work array WORKL is not sufficient." ;
                break ;
                case -8 :
                    tMessage = "Error return from trid. eigenvalue calculation; informational error from LAPACK routine dsteqr." ;
                break ;
                case -9 :
                    tMessage = "Starting vector is zero." ;
                break ;
                case -10 :
                    tMessage = "IPARAM(7) must be 1,2,3,4,5." ;
                break ;
                case -11 :
                    tMessage = "IPARAM(7) = 1 and BMAT = 'G' are incompatible." ;
                break ;
                case -12 :
                    tMessage = "IPARAM(1) must be equal to 0 or 1." ;
                break ;

                // NOT the same as dnaupd's -13, which is about HOWMNY
                case -13 :
                    tMessage = "NEV and WHICH = 'BE' are incompatible." ;
                break ;
                case -9999 :
                    tMessage = "Could not build a Lanczos factorization." ;
                break ;

                // 100 and up are raised by the BELFEM drivers, not by ARPACK,
                // and carry the same meanings as in check_naupd
                case 100 :
                    tMessage = "The driver received an unexpected reverse communication request." ;
                break ;
                case 101 :
                    tMessage = "The local row blocks do not sum to the global row count, "
                               "or the global row count is not the same on every rank." ;
                break ;
                case 102 :
                    tMessage = "At least one rank owns no rows. PARPACK needs a positive "
                               "local dimension on every rank of the communicator." ;
                break ;
                case 103 :
                    tMessage = "The local number of nonzeros contradicts the row pointers "
                               "( pointers( nloc + 1 ) - 1 must equal nnz )." ;
                break ;
                case 104 :
                    tMessage = "The job flag must be 0 ( smallest magnitude ) or 1 ( largest magnitude )." ;
                break ;
                case 105 :
                    tMessage = "An MPI call failed while the row map was being built." ;
                break ;
                default:
                    tMessage = "Unknown Error" ;
            }

            BELFEM_ERROR( false,
                "ARPACK-ng has thrown error %i in dsaupd : %s",
                    ( int ) aInfo,
                    tMessage.c_str() );
        }

//------------------------------------------------------------------------------

        void
        check_seupd( const int_t aInfo )
        {
            if( aInfo == 0 )
            {
                return ;
            }

            // the extraction side of a non-converged dsaupd, same tier as
            // check_neupd's -14: an expected outcome, not a fault
            if( aInfo == -14 )
            {
                message( InfoLevel::Verbose,
                    "ARPACK-ng dseupd found no eigenvalue of sufficient accuracy.\n"
                    "                 No eigenvalue is available for this step.\n" );
                return ;
            }

            string tMessage ;

            switch ( aInfo )
            {
                case -1 :
                    tMessage = "N must be positive." ;
                break ;
                case -2 :
                    tMessage = "NEV must be positive." ;
                break ;
                case -3 :
                    tMessage = "NCV-NEV >= 2 and less than or equal to N." ;
                break ;
                case -5 :
                    tMessage = "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'." ;
                break ;
                case -6 :
                    tMessage = "BMAT must be one of 'I' or 'G'." ;
                break ;
                case -7 :
                    tMessage = "Length of private work WORKL array is not sufficient." ;
                break ;
                case -8 :
                    tMessage = "Error return from trid. eigenvalue calculation; informational error from LAPACK routine dsteqr." ;
                break ;
                case -9 :
                    tMessage = "Starting vector is zero." ;
                break ;
                case -10 :
                    tMessage = "IPARAM(7) must be 1,2,3,4,5." ;
                break ;
                case -11 :
                    tMessage = "IPARAM(7) = 1 and BMAT = 'G' are incompatible." ;
                break ;
                case -12 :
                    tMessage = "NEV and WHICH = 'BE' are incompatible." ;
                break ;
                case -14 :
                    tMessage = "DSAUPD did not find any eigenvalues to sufficient accuracy." ;
                break ;
                case -15 :
                    tMessage = "HOWMNY must be one of 'A' or 'S' if RVEC = .true." ;
                break ;
                case -16 :
                    tMessage = "HOWMNY = 'S' not yet implemented." ;
                break ;
                case -17 :
                    tMessage = "DSEUPD got a different count of the number of converged Ritz values than DSAUPD got." ;
                break ;
                default:
                    tMessage = "Unknown Error" ;
            }

            BELFEM_ERROR( false,
                "ARPACK-ng has thrown error %i in dseupd : %s",
                    ( int ) aInfo,
                    tMessage.c_str() );
        }

//------------------------------------------------------------------------------
    }
}
