//
// Created by Christian Messe on 16.07.20.
//

#include "en_SolverEnums.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    string
    to_string( const SolverType aSolverType )
    {
        switch ( aSolverType )
        {
            case( SolverType::UMFPACK ) :
            {
                return "UMFPACK" ;
            }
            case( SolverType::MUMPS ) :
            {
                return "MUMPS" ;
            }
            case( SolverType::STRUMPACK ) :
            {
                return "STRUMPACK" ;
            }
            case( SolverType::PARDISO ) :
            {
                return "PARDISO" ;
            }
            case( SolverType::PETSC ) :
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

        if( tString == "umfpack" )
        {
            return SolverType::UMFPACK ;
        }
        else if( tString == "mumps" )
        {
            return SolverType::MUMPS ;
        }
        else if( tString == "strumpack" )
        {
            return SolverType::STRUMPACK ;
        }
        else if ( tString == "pardiso" )
        {
            return SolverType::PARDISO ;
        }
        else if ( tString == "petsc" )
        {
            return SolverType::PETSC ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown solver: %s", aString.c_str() );
            return SolverType::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------

    string
    to_string( const EulerMethod aEulerMethod )
    {
        switch ( aEulerMethod )
        {
            case( EulerMethod::ForwardExplicit ) :
            {
                return "ForwardExplicit" ;
            }
            case( EulerMethod::CrankNicolson ) :
            {
                return "CrankNicolson" ;
            }
            case( EulerMethod::Galerkin ) :
            {
                return "Galerkin" ;
            }
            case( EulerMethod::BackwardImplicit ) :
            {
                return "BackwardImplicit" ;
            }
            default:
            {
                return "unknown" ;
            }
        }
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
                return "gasm" ;
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

        if( tString == "none" )
        {
            return Preconditioner::NONE ;
        }
        else if ( tString == "asm" )
        {
            return Preconditioner::ASM ;
        }
        else if ( tString ==  "gasm" )
        {
            return Preconditioner::GASM ;
        }
        else if ( tString ==  "gamg" )
        {
            return Preconditioner::GAMG ;
        }
        else if ( tString ==  "jacobi" )
        {
            return  Preconditioner::JACOBI ;
        }
        else if ( tString ==  "bjacobi" )
        {
            return  Preconditioner::BJACOBI ;
        }
        else if ( tString ==  "lu" )
        {
            return Preconditioner::LU;
        }
        else if ( tString ==  "ilu" )
        {
            return Preconditioner::ILU;
        }
        else if ( tString ==  "icc" )
        {
            return Preconditioner::ICC;
        }
        else if ( tString ==  "hmg" )
        {
            return Preconditioner::HMG;
        }
        else if ( tString ==  "spai" )
        {
            return Preconditioner::SPAI;
        }
        else
        {
            BELFEM_ERROR( false, "unknown preconditioner: %s", aString.c_str() );
            return Preconditioner::UNDEFINED;
        }
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
                return "gmres" ;
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

        if( tString ==  "preonly" )
        {
            return KrylovMethod::PREONLY ;
        }
        if( tString ==  "auto" )
        {
            return KrylovMethod::AUTO ;
        }
        else if ( tString == "bcgs" )
        {
            return KrylovMethod::BCGS ;
        }
        else if ( tString == "cg" )
        {
            return KrylovMethod::CG ;
        }
        else if ( tString == "cgs" )
        {
            return KrylovMethod::CGS ;
        }
        else if ( tString == "ibcgs" )
        {
            return KrylovMethod::IBCGS ;
        }
        else if ( tString == "gmres" )
        {
            return KrylovMethod::GMRES ;
        }
        else if ( tString == "tfqmr" )
        {
            return KrylovMethod::TFQMR ;
        }
        else if ( tString == "bcgs" )
        {
            return KrylovMethod::BCGS ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown krylov method: %s", aString.c_str() );
            return KrylovMethod::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------

    string
    to_string( const SerialReodrdering aSerialReodrdering )
    {
        switch ( aSerialReodrdering )
        {
            case ( SerialReodrdering::AMD ) :
            {
                return "amd" ;
            }
            case ( SerialReodrdering::USERPIVOT ) :
            {
                return "userpivot" ;
            }
            case ( SerialReodrdering::AMF )  :
            {
                return "amf" ;
            }
            case ( SerialReodrdering::SCOTCH )  :
            {
                return "scotch" ;
            }
            case ( SerialReodrdering::PORD )  :
            {
                return "prod" ;
            }
            case ( SerialReodrdering::METIS )  :
            {
                return "metis" ;
            }
            case ( SerialReodrdering::QAMD )  :
            {
                return "qamd" ;
            }
            case( SerialReodrdering::AUTOMATIC )  :
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

    SerialReodrdering
    serial_reordering( const string & aString )
    {
        string tString = string_to_lower( aString );

        if( tString == "amd" )
        {
            return SerialReodrdering::AMD ;
        }
        else if ( tString == "userpivot" )
        {
            return SerialReodrdering::USERPIVOT ;
        }
        else if ( tString == "amf" )
        {
            return SerialReodrdering::AMF ;
        }
        else if ( tString == "scotch" )
        {
            return SerialReodrdering::SCOTCH ;
        }
        else if ( tString == "prod" )
        {
            return SerialReodrdering::PORD ;
        }
        else if ( tString == "metis" )
        {
            return SerialReodrdering::METIS ;
        }
        else if ( tString == "qamd" )
        {
            return SerialReodrdering::QAMD ;
        }
        else if ( tString == "automatic" || tString == "auto" )
        {
            return SerialReodrdering::AUTOMATIC ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown serial reordering method: %s", aString.c_str() );
            return SerialReodrdering::UNDEFINED ;
        }
    }


//------------------------------------------------------------------------------

    std::string
    to_string( const ParallelReodrdering aParallelReodrdering )
    {
        switch ( aParallelReodrdering )
        {
            case ( ParallelReodrdering::AUTOMATIC ) :
            {
                return "automatic" ;
            }
            case( ParallelReodrdering::PTSCOTCH ) :
            {
                return "ptscotch" ;
            }
            case( ParallelReodrdering::PARMETIS ) :
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

    ParallelReodrdering
    parallel_reordering( const string & aString )
    {
        string tString = string_to_lower( aString );

        if( tString == "automatic" || tString == "auto" )
        {
            return ParallelReodrdering::AUTOMATIC ;
        }
        else if ( tString == "ptscotch" )
        {
            return ParallelReodrdering::PTSCOTCH ;
        }
        else if ( tString == "parmetis" )
        {
            return ParallelReodrdering::PARMETIS ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown parallel reordering method: %s", aString.c_str() );
            return ParallelReodrdering::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------

    string
    to_string( const BlockLowRanking aBLR )
    {
        switch( aBLR )
        {
            case( BlockLowRanking::Off ) :
            {
                return "off" ;
            }
            case( BlockLowRanking::Automatic ) :
            {
                return "automatic" ;
            }
            case( BlockLowRanking::FactorizationAndSolution ) :
            {
                return "FactorizationAndSolution" ;
            }
            case( BlockLowRanking::FactorizationOnly ) :
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

    BlockLowRanking
    block_low_ranking( const string & aString )
    {
        string tString = string_to_lower( aString );

        if( tString == "off" )
        {
            return BlockLowRanking::Off ;
        }
        else if( tString == "automatic" || tString == "auto" )
        {
            return BlockLowRanking::Automatic ;
        }
        else if ( tString == "factorizationandsolution" )
        {
            return BlockLowRanking::FactorizationAndSolution ;
        }
        else if ( tString == "factorizationonly" )
        {
            return BlockLowRanking::FactorizationOnly ;
        }
        else
        {
            BELFEM_ERROR( false, "unknown block low ranking method: %s", aString.c_str() );
            return BlockLowRanking::UNDEFINED ;
        }
    }

//------------------------------------------------------------------------------
}