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

#include "en_SolverEnums.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    string
    to_string( const SolverType aSolverType )
    {
        switch ( aSolverType )
        {
            case SolverType::UMFPACK  :
            {
                return "UMFPACK" ;
            }
            case SolverType::SUPERLU  :
            {
                return "SUPERLU" ;
            }
            case SolverType::MUMPS  :
            {
                return "MUMPS" ;
            }
            case SolverType::STRUMPACK :
            {
                return "STRUMPACK" ;
            }
            case SolverType::PARDISO :
            {
                return "PARDISO" ;
            }
            case SolverType::PETSc :
            {
                return "PETSc" ;
            }
            default:
            {
                return "UNKNOWN" ;
            }
        }
    }

//--------------------------------------------------------------------------

    SolverType
    solver_type( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast<  uint >( SolverType::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< SolverType >( k ) ) ) == tString )
            {
                return static_cast< SolverType >( k ) ;
            }
        }
        BELFEM_ERROR( false, "unknown solver: %s", aString.c_str() );
        return SolverType::UNDEFINED ;
    }

//------------------------------------------------------------------------------

    string
    to_string( const DistributedMatrixType aDistributedMatrixType )
    {
        switch ( aDistributedMatrixType )
        {
            case DistributedMatrixType::CSR :
            {
                return "csr" ;
            }
            case DistributedMatrixType::CSC :
            {
                return "csc" ;
            }
            case DistributedMatrixType::AIJ :
            {
                return "aij" ;
            }
            default:
            {
                return "unknown" ;
            }
        }
    }

    DistributedMatrixType
    distributed_matrix_type( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k < static_cast< uint >( DistributedMatrixType::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< DistributedMatrixType >( k ) ) ) == tString )
            {
                return static_cast< DistributedMatrixType >( k );
            }
        }
        BELFEM_ERROR( false, "unknown distributed matrix format: %s", aString.c_str() );
        return DistributedMatrixType::UNDEFINED ;
    }

    DistributedMatrixType
    distributed_matrix_format( const string & aString );

//------------------------------------------------------------------------------

    string
    to_string( const EulerMethod aEulerMethod )
    {
        switch ( aEulerMethod )
        {
            case EulerMethod::ForwardExplicit :
            {
                return "ForwardExplicit" ;
            }
            case EulerMethod::CrankNicolson :
            {
                return "CrankNicolson" ;
            }
            case EulerMethod::Galerkin :
            {
                return "Galerkin" ;
            }
            case EulerMethod::BackwardDifference1 :
            {
                return "BDF1" ;
            }
            case EulerMethod::BackwardDifference2 :
            {
                return "BDF2" ;
            }
            case EulerMethod::BackwardDifference3 :
            {
                return "BDF3" ;
            }
            case EulerMethod::BackwardDifference4 :
            {
                return "BDF4" ;
            }
            case EulerMethod::BackwardDifference5 :
            {
                return "BDF5" ;
            }
            default:
            {
                return "unknown" ;
            }
        }
    }

    EulerMethod
    euler_method( const string & aString )
    {
        string tString = string_to_lower( aString );
        for ( uint k=0; k<static_cast<  uint >( EulerMethod::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< EulerMethod >( k ) ) ) == tString )
            {
                return static_cast< EulerMethod >( k ) ;
            }
        }

        BELFEM_ERROR( false, "unknown euler method: %s", aString.c_str() );
        return EulerMethod::UNDEFINED ;
    }

//------------------------------------------------------------------------------

    string
    to_string( const Preconditioner aPreconditioner )
    {
        switch ( aPreconditioner )
        {
            case( Preconditioner::NONE ) :
            {
                return "none" ;
            }
            case( Preconditioner::ASM ) :
            {
                return "asm" ;
            }
            case( Preconditioner::GASM ) :
            {
                return "asm" ;  // mapped onto PCASM; the restricted variant is selected by PCASMSetType() in the wrapper
            }
            case( Preconditioner::GAMG ) :
            {
                return "gamg" ;
            }
            case( Preconditioner::JACOBI ) :
            {
                return "jacobi" ;
            }
            case( Preconditioner::BJACOBI ) :
            {
                return "bjacobi" ;
            }
            case( Preconditioner::LU ) :
            {
                return "lu" ;
            }
            case( Preconditioner::ILU ) :
            {
                return "ilu" ;
            }
            case( Preconditioner::ICC ) :
            {
                return "icc" ;
            }
            case( Preconditioner::HMG ) :
            {
                return "hmg" ;
            }
            case( Preconditioner::SPAI ) :
            {
                return "spai" ;
            }
            default:
            {
                return "unknown";
            }
        }
    }

//------------------------------------------------------------------------------

    Preconditioner
    preconditioner( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast<  uint >( Preconditioner::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< Preconditioner >( k ) ) ) == tString )
            {
                return static_cast< Preconditioner >( k ) ;
            }
        }

        BELFEM_ERROR( false, "unknown preconditioner: %s", aString.c_str() );
        return Preconditioner::UNDEFINED;
    }

//------------------------------------------------------------------------------

    string
    to_string( const KrylovMethod aKrylovMethod )
    {
        switch ( aKrylovMethod )
        {
            case( KrylovMethod::PREONLY ) :
            {
                return "preonly" ;
            }
            case( KrylovMethod::BCGS ) :
            {
                return "bcgs" ;
            }
            case( KrylovMethod::CG ) :
            {
                return "cg" ;
            }
            case( KrylovMethod::CGS ) :
            {
                return "cgs" ;
            }
            case( KrylovMethod::IBCGS ) :
            {
                return "ibcgs" ;
            }
            case( KrylovMethod::GMRES ) :
            {
                return "gmres" ;
            }
            case( KrylovMethod::TFQMR ) :
            {
                return "tfqmr" ;
            }
            case( KrylovMethod::AUTO ) :
            {
                return "auto" ;
            }
            default:
            {
                return "unknown";
            }
        }
    }

//----------------------------------------------------------------------------

    KrylovMethod
    krylov_method( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast<  uint >( KrylovMethod::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< KrylovMethod >( k ) ) ) == tString )
            {
                return static_cast< KrylovMethod >( k ) ;
            }
        }

        BELFEM_ERROR( false, "unknown krylov method: %s", aString.c_str() );
        return KrylovMethod::UNDEFINED ;
    }

    string
    to_string( const ReorderingMethod aReorderingMethod )
    {
        switch (  aReorderingMethod )
        {
            case ReorderingMethod::NATURAL :
            {
                return "natural" ;
            }
            case ReorderingMethod::METIS :
            {
                return "metis" ;
            }
            case ReorderingMethod::SCOTCH :
            {
                return "scotch" ;
            }
            case ReorderingMethod::AUTOMATIC :
            {
                return "automatic" ;
            }
            case ReorderingMethod::PARMETIS :
            {
                return "parmetis" ;
            }
            case ReorderingMethod::PTSCOTCH :
            {
                return "ptscotch" ;
            }
            default:
            {
                return "unknown" ;
            }
        }
    }

    ReorderingMethod
    reordering_method( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast<  uint >( ReorderingMethod::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< ReorderingMethod >( k ) ) ) == tString )
            {
                return static_cast< ReorderingMethod >( k ) ;
            }
        }

        BELFEM_ERROR( false, "unknown reordering method: %s", aString.c_str() );
        return ReorderingMethod::UNDEFINED ;
    }

//------------------------------------------------------------------------------

    string
    to_string( const CompressionMethod aCompressionMethod )
    {
        switch ( aCompressionMethod )
        {
            case CompressionMethod::OFF :
            {
                return "off" ;
            }
        case CompressionMethod::BLR :
            {
                return "blr" ;
            }
            case CompressionMethod::AUTOMATIC :
            {
                return "automatic" ;
            }
            default:
            {
                return "unknown" ;
            }
        }
    }

    CompressionMethod
    compression_method( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast<  uint >( CompressionMethod::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< CompressionMethod >( k ) ) ) == tString )
            {
                return static_cast< CompressionMethod >( k ) ;
            }
        }

        BELFEM_ERROR( false, "unknown compression method: %s", aString.c_str() );
        return CompressionMethod::UNDEFINED ;
    }

//------------------------------------------------------------------------------

    string
    to_string( const MumpsSerialReodrdering aSerialReodrdering )
    {
        switch ( aSerialReodrdering )
        {
            case ( MumpsSerialReodrdering::AMD ) :
            {
                return "amd" ;
            }
            case ( MumpsSerialReodrdering::USERPIVOT ) :
            {
                return "userpivot" ;
            }
            case ( MumpsSerialReodrdering::AMF )  :
            {
                return "amf" ;
            }
            case ( MumpsSerialReodrdering::SCOTCH )  :
            {
                return "scotch" ;
            }
            case ( MumpsSerialReodrdering::PORD )  :
            {
                return "prod" ;
            }
            case ( MumpsSerialReodrdering::METIS )  :
            {
                return "metis" ;
            }
            case ( MumpsSerialReodrdering::QAMD )  :
            {
                return "qamd" ;
            }
            case( MumpsSerialReodrdering::AUTOMATIC )  :
            {
                return "automatic" ;
            }
            default :
            {
                return "unknown" ;
            }
        }
    }

//------------------------------------------------------------------------------

    MumpsSerialReodrdering
    serial_reordering( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast< uint >( MumpsSerialReodrdering::UNDEFINED ); ++k )
        {
            if ( string_to_lower(
                to_string( static_cast< MumpsSerialReodrdering >( k ) ) ) == tString )
            {
                return static_cast< MumpsSerialReodrdering >( k ) ;
            }
        }
        BELFEM_ERROR( false, "unknown serial reordering method: %s", aString.c_str() );
        return MumpsSerialReodrdering::UNDEFINED ;

    }


//------------------------------------------------------------------------------

    std::string
    to_string( const MumpsParallelReodrdering aParallelReodrdering )
    {
        switch ( aParallelReodrdering )
        {
            case ( MumpsParallelReodrdering::AUTOMATIC ) :
            {
                return "automatic" ;
            }
            case( MumpsParallelReodrdering::PTSCOTCH ) :
            {
                return "ptscotch" ;
            }
            case( MumpsParallelReodrdering::PARMETIS ) :
            {
                return "parmetis" ;
            }
            default:
            {
                return "undefined" ;
            }
        }
    }

//------------------------------------------------------------------------------

    MumpsParallelReodrdering
    parallel_reordering( const string & aString )
    {
        string tString = string_to_lower( aString );

        if( tString == "automatic" || tString == "auto" )
        {
            return MumpsParallelReodrdering::AUTOMATIC ;
        }
        else if ( tString == "ptscotch" )
        {
            return MumpsParallelReodrdering::PTSCOTCH ;
        }
        else if ( tString == "parmetis" )
        {
            return MumpsParallelReodrdering::PARMETIS ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown parallel reordering method: %s", aString.c_str() );
            return MumpsParallelReodrdering::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------

    string
    to_string( const MumpsBlockLowRanking aBLR )
    {
        switch( aBLR )
        {
            case( MumpsBlockLowRanking::Off ) :
            {
                return "off" ;
            }
            case( MumpsBlockLowRanking::Automatic ) :
            {
                return "automatic" ;
            }
            case( MumpsBlockLowRanking::FactorizationAndSolution ) :
            {
                return "FactorizationAndSolution" ;
            }
            case( MumpsBlockLowRanking::FactorizationOnly ) :
            {
                return "FactorizationOnly" ;
            }
            default :
            {
                return "undefined" ;
            }
        }
    }
//------------------------------------------------------------------------------

    MumpsBlockLowRanking
    block_low_ranking( const string & aString )
    {
        string tString = string_to_lower( aString );

        for ( uint k=0; k<static_cast< uint >( MumpsBlockLowRanking::UNDEFINED ); ++k )
        {
            if ( string_to_lower( to_string( static_cast< MumpsBlockLowRanking >( k ) ) ) == tString )
            {
                return static_cast< MumpsBlockLowRanking >( k ) ;
            }
        }
        BELFEM_ERROR( false, "unknown block low ranking method: %s", aString.c_str() );
        return MumpsBlockLowRanking::UNDEFINED ;
    }

//------------------------------------------------------------------------------
}