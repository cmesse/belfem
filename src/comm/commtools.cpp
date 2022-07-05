//
// Created by Christian Messe on 2019-01-18.
//

#include "commtools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    proc_t
    comm_size()
    {
        return gComm.size();
    }

//------------------------------------------------------------------------------

    proc_t
    comm_rank()
    {
        return gComm.rank();
    }

//------------------------------------------------------------------------------

    void
    comm_barrier()
    {
#ifdef BELFEM_MPI
        MPI_Barrier( gComm.world() );
#endif
    }

//------------------------------------------------------------------------------

    void
    create_master_commlist( Vector< proc_t > & aCommList )
    {
        // number of procs
        proc_t tNumProcs = comm_size() ;
        proc_t tMyRank = comm_rank() ;

        if( tNumProcs > 1 )
        {
            // initialize counter
            uint tCount = 0 ;

            // allocate memory
            aCommList.set_size( tNumProcs-1 );

            // loop over all procs
            for( proc_t p = 0; p<tNumProcs; ++p )
            {
                // we do not communicate with ourselves
                if( p != tMyRank )
                {
                    aCommList( tCount++ ) = p ;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    int
    comm_tag( const proc_t aSource, const proc_t aTarget )
    {
        proc_t tMin = ( aSource < aTarget ) ? ( aSource ) : ( aTarget );
        proc_t tMax = ( aSource > aTarget ) ? ( aSource ) : ( aTarget );

        if ( tMin == aTarget )
        {
            return 1024 * ( tMax * gComm.size() + tMin ) + 1;
        }
        else
        {
            return 1024 * ( tMax * gComm.size() + tMin ) + 513;
        }
    }

//------------------------------------------------------------------------------

    Vector< int >
    split_message( const index_t aMessageLength )
    {
        // number of full packages
        index_t tFull = aMessageLength / gMaxArrayLength;

        Vector< int > aSteps( tFull + 1 );

        for ( index_t k = 0; k < tFull; ++k )
        {
            aSteps( k ) = gMaxArrayLength;
        }

        // length of final message
        aSteps( tFull ) = aMessageLength - tFull * gMaxArrayLength;

        return aSteps;
    }

}