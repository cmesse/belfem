//
// Created by christian on 8/16/21.
//

#include "typedefs.hpp"
#include "commtools.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_HDF5.hpp"
#include "cl_SpMatrix.hpp"
#include "petsctools.hpp"

using namespace belfem ;

Communicator gComm;
Logger       gLog( 3 );
//------------------------------------------------------------------------------


int main( int    argc,
          char * argv[] )
 {
    // create communicator
    gComm = Communicator( argc, argv );
#ifdef BELFEM_PETSC
    //------------------------------------------------------------------------------

    // matrix wrapper
    Mat tPetscMatrix ;

    Vector< PetscInt > mMyPointers ;
    Vector< PetscInt > mMyColumns ;
    Vector< PetscInt > mMyRows ;
    Vector< real >     mMyValues ;
    PetscInt mMyNumRows ;
    PetscInt mMyNumNNZ ;
    PetscInt mNumRows = 5 ;
    PetscInt mNumCols = 5 ;
    PetscInt mMyRowOffset = 0 ;

    if( gComm.rank() == 0 )
    {
        // load matrix
        // SpMatrix aMatrix( "matrix.hdf5");

        Matrix <real> tA( 5, 5, 0.0 );

        tA( 0, 0 ) =  1.0;
        tA( 0, 1 ) = -1.0;
        tA( 0, 3 ) = -3.0;
        tA( 1, 0 ) = -2.0;
        tA( 1, 1 ) =  5.0;
        tA( 2, 2 ) =  4.0;
        tA( 2, 3 ) =  6.0;
        tA( 2, 4 ) =  4.0;
        tA( 3, 0 ) = -4.0;
        tA( 3, 2 ) =  2.0;
        tA( 3, 3 ) =  7.0;
        tA( 4, 1 ) =  8.0;
        tA( 4, 4 ) = -5.0;


        SpMatrix aMatrix( tA, SpMatrixType::CSR );

        Vector <real> tY = { -13., 8., 56., 30., -9. };
        Vector <real> tX( 5, 0.1 );

        // get the number of rows
        PetscInt tNumAllRows =  aMatrix.n_cols();

        // get the number of procs
        PetscInt tNumProcs = gComm.size() ;

        // create a communication list
        Vector< proc_t > tCommList( tNumProcs );

        for( proc_t p=0; p<tNumProcs; ++p )
        {
            tCommList( p ) = p ;
        }
        tCommList.print("commlist");

        // make sure that matrix is cpp based and has both indices
        aMatrix.create_coo_indices();

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // define how many rows and cols exist per proc
        Vector< PetscInt > tNumRowsPerProc( tNumProcs, tNumAllRows / tNumProcs );
        Vector< PetscInt > tNumNnzProc( tNumProcs );

        // the last entry gets one more
        tNumRowsPerProc( tNumProcs - 1 ) += tNumAllRows % tNumProcs ;
        tNumRowsPerProc.print("NumRowsPerProc");

        Vector< PetscInt > tRowOffsetsPerProc( tNumProcs + 1, 0 );

        // now we must create the row offsets
        PetscInt tOff = 0 ;
        for( PetscInt p=0; p<tNumProcs; ++p )
        {
            tRowOffsetsPerProc( p ) += tOff ;
            tOff += tNumRowsPerProc( p );
        }
        tRowOffsetsPerProc( tNumProcs ) = tOff ;

        tRowOffsetsPerProc.print("RowOffsets");

        // split values
        Cell< Vector< PetscInt > > tAllPointers( tNumProcs, Vector< PetscInt >() );
        Cell< Vector< PetscInt > > tAllRows( tNumProcs, Vector< PetscInt >() );
        Cell< Vector< PetscInt > > tAllColumns( tNumProcs, Vector< PetscInt >() );
        Cell< Vector< PetscReal > > tAllValues( tNumProcs, Vector< PetscReal >() );
        // tOff = 0 ;

        int * tMatPointers =  aMatrix.pointers() ;
        //int * tRows =  aMatrix.rows() ;
        int * tMatCols =  aMatrix.cols() ;
        int * tMatRows =  aMatrix.rows() ;
        real * tMatValues = aMatrix.data() ;

        PetscInt tPointerOffset = 0 ;
        PetscInt tColumnOffset = 0 ;

        for( PetscInt p=0; p<tNumProcs; ++p )
        {
            // grab vector with pointers
            Vector< PetscInt > & tPointers = p==0 ? mMyPointers : tAllPointers( p );
            Vector< PetscInt > & tColumns = p==0 ? mMyColumns : tAllColumns( p );
            Vector< PetscInt > & tRows = p==0 ? mMyRows : tAllRows( p );
            Vector< real > & tValues = p==0 ? mMyValues : tAllValues( p );

            // number of nonzeros per proc
            PetscInt & tNNZ = p==0 ? mMyNumNNZ : tNumNnzProc( p );

            tNNZ = tMatPointers[ tRowOffsetsPerProc( p + 1 ) ] - tMatPointers[ tRowOffsetsPerProc( p ) ] ;

            PetscInt & tNumRows = p==0 ? mMyNumRows : tNumRowsPerProc( p );
            tNumRows = tNumRowsPerProc( p  ) ;

            tPointers.set_size( tNumRows + 1 );

            for( PetscInt i = 0; i<=tNumRows; ++i )
            {
                tPointers( i ) = tMatPointers[ tRowOffsetsPerProc( p ) + i ] - tPointerOffset ;
            }

            tPointerOffset += tPointers( tNumRows ) ;

            // grab vector with columns

            tColumns.set_size( tNNZ );
            tRows.set_size( tNNZ );

            tValues.set_size( tNNZ );

            for( PetscInt k=0; k<tNNZ; ++k )
            {
                tColumns( k ) = tMatCols[ tColumnOffset + k ];
                tRows( k ) = tMatRows[ tColumnOffset + k ];
                tValues( k ) = tMatValues[ tColumnOffset + k ];
            }

            tColumnOffset += tNNZ ;
            tPointers.print( "tPointers" );
            tColumns.print("tColumns");

        }


        // send list to other procs
        send( tCommList, tNumRowsPerProc );
        send( tCommList, tNumNnzProc );

        // length of tRowOffsetsPerProc does not match, do this to fix
        Vector< PetscInt > tRowOff( tNumProcs );
        for( PetscInt p=0; p<tNumProcs; ++p )
        {
            tRowOff( p ) = tRowOffsetsPerProc( p );
        }
        send( tCommList, tRowOff );

        send( tCommList, tAllPointers ) ;
        send( tCommList, tAllRows );
        send( tCommList, tAllColumns );
        send( tCommList, tAllValues );
    }
    else
    {
        receive( 0, mMyNumRows );
        receive( 0, mMyNumNNZ );
        receive( 0, mMyRowOffset );
        receive( 0, mMyPointers );
        receive( 0, mMyRows );
        receive( 0, mMyColumns );
        receive( 0, mMyValues );
    }

     mMyPointers.print("p");
     mMyColumns.print("c");
     mMyValues.print("v");
     MatCreate( gComm.world(), &tPetscMatrix );
     MatSetSizes( tPetscMatrix, mMyNumRows, mNumCols, mNumRows, mNumCols );
     MatSetType( tPetscMatrix, MATMPIAIJ );

     MatMPIAIJSetPreallocationCSR( tPetscMatrix, mMyPointers.data(), mMyColumns.data(), NULL );

     PetscInt tOffset = 0 ;
     PetscInt * tColumns = mMyColumns.data() ;
     PetscReal * tData = mMyValues.data() ;


     for( PetscInt r=0; r<mMyNumRows; ++r )
     {
         // compute num cols
         PetscInt tNumCols = mMyPointers( r+1 ) - mMyPointers( r );

         // get index of this row
         PetscInt tThisRow = r + mMyRowOffset ;

         // insert columns
         MatSetValues(
                 tPetscMatrix,
                 1,
                 &tThisRow,
                 tNumCols,
                 &tColumns[ tOffset ],
                 &tData[ tOffset ],
                 INSERT_VALUES );

         // increment offset
         tOffset += tNumCols ;
     }
     MatAssemblyBegin(tPetscMatrix,MAT_FINAL_ASSEMBLY);
     MatAssemblyEnd(tPetscMatrix,MAT_FINAL_ASSEMBLY);

     MatView( tPetscMatrix, PETSC_VIEWER_STDOUT_( gComm.world() ) );
    //------------------------------------------------------------------------------
#endif
    // close communicator
    return gComm.finalize();
}