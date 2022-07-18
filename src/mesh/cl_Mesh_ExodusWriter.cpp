//
// Created by Christian Messe on 2019-07-28.
//
// https://gsjaardema.github.io/seacas/html/deprecated.html
#include "cl_Mesh_ExodusWriter.hpp"

#include <cstring>

#ifdef BELFEM_EXODUS
#include "netcdf.h"
#include "exodusII.h"
#endif

#include "cl_Cell.hpp"
#include "cl_StringList.hpp"
#include "cl_Mesh_Field.hpp"
#include "commtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        ExodusWriter::ExodusWriter( Mesh * aMesh )
#ifdef BELFEM_EXODUS
          :  mMesh( aMesh )
#endif
        {

        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::save( const string & aPath )
        {
#ifdef BELFEM_EXODUS
            mPath = aPath;

            this->populate_header();

            this->populate_node_coords();
            this->populate_element_ids();
            this->populate_blocks();

            this->populate_element_connectivity();
            if( mNumSideSets > 0 )
            {
                this->populate_sidesets();
            }

            this->populate_time();
            this->populate_global_variables();
            this->populate_fields();

            this->close_file();


#else
            BELFEM_ERROR( false, "Exodus is not linked in this runtime." );
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::check( const string & aRoutine )
        {
#ifdef BELFEM_EXODUS
            BELFEM_ERROR( mError == 0,
                    "Exodus call %s has thorwn an error: %d",
                    aRoutine.c_str(),
                    mError );
#endif
        }
//------------------------------------------------------------------------------
        void
        ExodusWriter::close_file()
        {
#ifdef BELFEM_EXODUS
            ex_close( mHandle );
#endif
        }

//------------------------------------------------------------------------------
        void
        ExodusWriter::populate_header()
        {
#ifdef BELFEM_EXODUS
            mCpuWordSize = sizeof( double );
            mIoWordSize = 8;

            // create a file
            mHandle = ex_create(
                    mPath.c_str(), /* filename path */
                    EX_CLOBBER,  /* create mode */
                    & mCpuWordSize, /* CPU float word size in bytes */
                    & mIoWordSize /* I/O float word size in bytes */
                    );

            mNumDim   = mMesh->number_of_dimensions();
            mNumNodes = mMesh->number_of_nodes();

            mNumElements = mMesh->number_of_elements();
            mNumBlocks   = mMesh->number_of_blocks();


            // count visible sidesets
            mNumSideSets = 0 ;
            for( SideSet * tSideSet : mMesh->sidesets() )
            {
                if( ! tSideSet->is_hidden() )
                {
                    ++mNumSideSets ;
                }
            }

            mError = ex_put_init( mHandle,
                    "Go Buffs!",
                    mNumDim,
                    mNumNodes,
                    mNumElements,
                    mNumBlocks,
                    mNumNodeSets,
                    mNumSideSets );

            // check error from exodus
            this->check( "ex_put_init");
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_node_coords()
        {
#ifdef BELFEM_EXODUS
            BELFEM_ASSERT( mMesh->number_of_dimensions() == 2 ||
                          mMesh->number_of_dimensions() == 3,
                          "Number of dimensions of mesh must be 2 or 3 ( is %u ).",
                                  ( unsigned int ) mMesh->number_of_dimensions() );


            double * tCoords = ( double * ) malloc( mNumNodes * sizeof( double ) );

            index_t tCount=0;

            for( mesh::Node * tNode : mMesh->nodes() )
            {
                tCoords[ tCount++ ] = tNode->x();
            }

            mError = ex_put_coord( mHandle, tCoords, NULL, NULL );

            // check error from exodus
            this->check( "ex_put_coord (x)");

            tCount = 0;

            for( mesh::Node * tNode : mMesh->nodes() )
            {
                tCoords[ tCount++ ] = tNode->y();
            }

            mError = ex_put_coord( mHandle, NULL, tCoords, NULL );

            this->check( "ex_put_coord (y)");

            if( mNumDim == 3 )
            {
                tCount=0;

                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    tCoords[ tCount++ ] = tNode->z();
                }

                // send coordinates to exodus
                mError = ex_put_coord( mHandle, NULL, NULL, tCoords );

                this->check( "ex_put_coord (z)");

                StringList tLabels( 3 );
                tLabels.push( "x" );
                tLabels.push( "y" );
                tLabels.push( "z" );

                mError = ex_put_coord_names( mHandle, tLabels.data() );

                this->check( "ex_put_coord_names (xyz)");
            }
            else
            {
                StringList tLabels( 2 );
                tLabels.push( "x" );
                tLabels.push( "y" );

                mError = ex_put_coord_names( mHandle, tLabels.data() );
                this->check( "ex_put_coord_names (xy)");
            }

            free( tCoords );
#endif
        }
//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_element_ids()
        {
#ifdef BELFEM_EXODUS
            int * tElementIDs = ( int * ) malloc( mNumElements * sizeof( int ) );
            int64_t tCount = 0;

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for( int b=0; b<mNumBlocks; ++b )
            {
                Block * tBlock = tBlocks( b );
                int64_t tNumElements = tBlock->number_of_elements();

                for( int64_t e=0; e<tNumElements; ++e )
                {
                    tElementIDs[ tCount++ ] = tBlock->element( e )->id();
                }
            }
            mError = ex_put_map ( mHandle, tElementIDs );
            this->check( "ex_put_map");
            free( tElementIDs );
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_blocks()
        {
#ifdef BELFEM_EXODUS
            StringList tLabels( mNumBlocks );

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for( int64_t b=0; b<mNumBlocks; ++b )
            {
                Block * tBlock = tBlocks( b );

                string tElementType= this->element_string_from_type(
                        tBlock->element( 0 )->type() );

                int tNumNodes = this->fix_num_nodes( tBlock->element_type() ) ;

                mError = ex_put_block(
                        mHandle,
                        EX_ELEM_BLOCK,
                        tBlock->id(),
                        tElementType.c_str(),
                        tBlock->number_of_elements(),
                        tNumNodes,
                        0,
                        0,
                        1);

                this->check( "ex_put_block");


                mError = ex_put_prop( mHandle,
                                      EX_ELEM_BLOCK,
                                      tBlock->id(),
                                      tBlock->label().c_str(),
                                      1);

                this->check( "ex_put_prop");

                tLabels.push( tBlock->label() );
            }

            mError = ex_put_names(
                    mHandle,
                    EX_ELEM_BLOCK,
                    tLabels.data() );

            this->check( "ex_put_names");
#endif
        }

//------------------------------------------------------------------------------

        string
        ExodusWriter::element_string_from_type( const ElementType & aType )
        {
            switch( geometry_type( ( aType ) ) )
            {
                case( GeometryType::VERTEX ):
                {
                    return "sphere";
                }
                case( GeometryType::LINE ):
                {
                    return "beam";
                }
                case( GeometryType::TRI ):
                {
                    return "tri";
                }
                case( GeometryType::QUAD ):
                {
                    return "quad";
                }
                case( GeometryType::TET ):
                {
                    return "tetra";
                }
                case( GeometryType::HEX ):
                {
                    return "hex";
                }
                case( GeometryType::PENTA ):
                {
                    return "wedge";
                }
                case( GeometryType::PYRA ):
                {
                    return "pyramid";
                }
                default:
                {
                    BELFEM_ERROR( false, "unknown element type" );
                    return "Unknown";
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_element_connectivity()
        {
#ifdef BELFEM_EXODUS

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for( int64_t b=0; b<mNumBlocks; ++b )
            {
                Block * tBlock = tBlocks( b );

                uint tNumNodes = this->fix_num_nodes( tBlock->element_type() ) ;

                int64_t tNumElements = tBlock->number_of_elements();

                int * tConnectivity = ( int * ) malloc( tNumNodes * tNumElements * sizeof( int ) );

                int64_t tCount = 0;

                for( int64_t e=0; e<tNumElements; ++e )
                {
                    Element * tElement = tBlock->element( e );

                    for( uint k=0; k<tNumNodes; ++k )
                    {
                        tConnectivity[ tCount++ ] = tElement->node( k )->index() + 1;
                    }
                }

                mError = ex_put_conn(
                        mHandle,
                        EX_ELEM_BLOCK,
                        tBlock->id(),
                        tConnectivity, 0, 0 );

                this->check( "ex_put_conn");

                free( tConnectivity );

            }
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_sidesets()
        {
#ifdef BELFEM_EXODUS
            Cell< SideSet * > & tSideSets = mMesh->sidesets() ;

            StringList tLabels( mNumSideSets );


            for( SideSet * tSideSet : tSideSets )
            {
                // check if sideset exists on exodus
                if ( !tSideSet->is_hidden() )
                {
                    // number of facets
                    int64_t tNumFacets = tSideSet->number_of_facets();

                    ex_put_set_param(
                            mHandle,
                            EX_SIDE_SET,
                            tSideSet->id(),
                            tNumFacets,
                            0 );

                    this->check( "ex_put_set_param" );

                    // allocate containers
                    int * tElementIDs = ( int * ) malloc( tNumFacets * sizeof( int ));
                    int * tSideSetIDs = ( int * ) malloc( tNumFacets * sizeof( int ));

                    for ( int64_t i = 0; i < tNumFacets; ++i )
                    {
                        tElementIDs[ i ] = tSideSet->facet_by_index( i )->master_id();
                    }

                    for ( int64_t i = 0; i < tNumFacets; ++i )
                    {
                        tSideSetIDs[ i ] = tSideSet->facet_by_index( i )->master_index() + 1;
                    }

                    mError = ex_put_set(
                            mHandle,
                            EX_SIDE_SET,
                            tSideSet->id(),
                            tElementIDs,
                            tSideSetIDs );

                    this->check( "ex_put_set" );

                    free( tElementIDs );
                    free( tSideSetIDs );

                    tLabels.push( tSideSet->label() );
                }
            }

            mError = ex_put_names(
                    mHandle,
                    EX_SIDE_SET,
                    tLabels.data());

            this->check( "ex_put_names");
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_fields()
        {
#ifdef BELFEM_EXODUS

            // count fields
            uint tNodeFieldCount = 0;
            uint tElementFieldCount = 0;

            // get number of fields from mesh
            uint tNumFields = mMesh->number_of_fields();

            for( uint k=0; k<tNumFields; ++k )
            {
                mesh::Field * tField = mMesh->field( k );

                if( tField->write_field_to_file() )
                {
                    switch ( tField->entity_type())
                    {
                        case ( EntityType::NODE ) :
                        {
                            ++tNodeFieldCount;
                            break;
                        }
                        case ( EntityType::ELEMENT ):
                        {
                            ++tElementFieldCount;
                            break;
                        }
                        default :
                        {
                            break;
                        }
                    }
                }
            }

            // allocate matrices
            Cell< mesh::Field * > tNodeFields( tNodeFieldCount, nullptr );
            Cell< mesh::Field * > tElementFields( tElementFieldCount, nullptr );

            // reset counters
            tNodeFieldCount = 0;
            tElementFieldCount = 0;

            // collect fields
            for( uint k=0; k<tNumFields; ++k )
            {
                mesh::Field * tField = mMesh->field( k );
                if( tField->write_field_to_file() )
                {
                    switch ( tField->entity_type())
                    {
                        case ( EntityType::NODE ) :
                        {
                            tNodeFields( tNodeFieldCount++ ) = tField;
                            break;
                        }
                        case ( EntityType::ELEMENT ):
                        {
                            tElementFields( tElementFieldCount++ ) = tField;
                            break;
                        }
                        default:
                        {
                            // do nothing
                            break;
                        }
                    }
                }
            }
            this->populate_node_fields( tNodeFields );
            this->populate_element_fields( tElementFields );


            tNodeFields.clear();
            tElementFields.clear();

            mError = ex_update( mHandle );
            this->check( "ex_update");
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter:: populate_time()
        {
#ifdef BELFEM_EXODUS

            mTimeValue = mMesh->time_stamp();

            /* write time value */
            mError = ex_put_time(
                    mHandle,
                    mTimeStep,
                    &mTimeValue );

            this->check( "ex_put_time" );
#endif
        }
//------------------------------------------------------------------------------
        void
        ExodusWriter::populate_global_variables()
        {
#ifdef BELFEM_EXODUS
            // get number of global variables from mesh
            uint tNumGlobalVariables = mMesh->number_of_global_variables();

            if ( tNumGlobalVariables > 0 )
            {
                // get field titles
                StringList tFieldLabels( tNumGlobalVariables );

                for( uint k=0; k<tNumGlobalVariables; ++k )
                {
                    tFieldLabels.push( mMesh->global_variable( k )->label() );
                }

                /* initialize field data container */
                mError = ex_put_variable_param(
                        mHandle,
                        EX_GLOBAL,
                        tNumGlobalVariables );

                this->check( "ex_put_variable_param (global)");

                /* write names */
                mError = ex_put_variable_names(
                        mHandle,
                        EX_GLOBAL,
                        tNumGlobalVariables,
                        tFieldLabels.data() );

                this->check( "ex_put_variable_names (global)");

                /* assemble values of global variables */
                real * tGlobalVars = ( real * ) malloc( tNumGlobalVariables * sizeof( real ) );

                for( uint k=0; k<tNumGlobalVariables; ++k )
                {
                    tGlobalVars[ k ] = mMesh->global_variable( k )->value();
                }

                /* push variables */
                ex_put_var(
                        mHandle,
                        mTimeStep,
                        EX_GLOBAL,
                        1,
                        0,
                        tNumGlobalVariables,
                        tGlobalVars );

                this->check( "ex_put_var (global)");

                // tidy up memory
                free( tGlobalVars );
            }
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_node_fields( Cell< mesh::Field * > & aFields )
        {
#ifdef BELFEM_EXODUS
            uint tNumNodes = mMesh->number_of_nodes();
            uint tNumNodeFields = 2;

            // count fields that are to be written
            for ( mesh::Field * tField : aFields )
            {
                if( tField->write_field_to_file() )
                {
                    ++tNumNodeFields ;
                }
            }

            int tID = 1;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // collect field names
            StringList tFieldLabels( tNumNodeFields );

            tFieldLabels.push("NodeID");
            tFieldLabels.push("NodeOwner");

            for ( mesh::Field * tField : aFields )
            {
                if( tField->write_field_to_file() )
                {
                    tFieldLabels.push( tField->label() );
                }
            }

            /*  initialize data container */
            mError = ex_put_variable_param(
                    mHandle,
                    EX_NODAL,
                    tNumNodeFields );

            this->check( "ex_put_variable_param (node)");

            /* write names */
            mError = ex_put_variable_names(
                    mHandle,
                    EX_NODAL,
                    tNumNodeFields,
                    tFieldLabels.data() );

            this->check( "ex_put_variable_param (node)");

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // create field with node ids
            real * tNodeData = ( real * ) malloc( tNumNodes * sizeof( real ) );

            uint k=0;

            for( mesh::Node * tNode : mMesh->nodes() )
            {
                tNodeData[ k++ ] = ( real )  tNode->id();
            }

            // write node IDs
            /* write nodal variables */
            mError = ex_put_var(
                    mHandle,
                    mTimeStep,
                    EX_NODAL,
                    tID++,
                    1,
                    tNumNodes,
                    tNodeData );

            this->check( "ex_put_var (node IDs)");

            k = 0;
            for( mesh::Node * tNode : mMesh->nodes() )
            {
                tNodeData[ k++ ] = ( real )  tNode->owner();
            }

            // write node IDs
            /* write nodal variables */
            mError = ex_put_var(
                    mHandle,
                    mTimeStep,
                    EX_NODAL,
                    tID++,
                    1,
                    tNumNodes,
                    tNodeData );

            this->check( "ex_put_var (node owners)");

            free( tNodeData );

            // write other fields
            for ( mesh::Field * tField : aFields )
            {
                // check if field is to be written
                if( tField->write_field_to_file() )
                {
                    // write nodal variables
                    mError = ex_put_var(
                            mHandle,
                            mTimeStep,
                            EX_NODAL,
                            tID++,
                            1,
                            mMesh->number_of_nodes(),
                            tField->data().data() );

                    this->check( "ex_put_var (node data)" );
                }
            }
#endif
        }

//------------------------------------------------------------------------------

        void
        ExodusWriter::populate_element_fields( Cell< mesh::Field * > & aFields )
        {
#ifdef BELFEM_EXODUS

            uint tNumBlocks = mMesh->number_of_blocks();

            int tID = 0;

            // count fields with write flag
            uint tNumElementFields = 4;

            for ( mesh::Field * tField : aFields )
            {
                if( tField->write_field_to_file() )
                {
                    ++ tNumElementFields ;
                }
            }

            // get field titles
            StringList tFieldLabels( tNumElementFields );

            tFieldLabels.push( "ElementID");
            tFieldLabels.push( "ElementOwner");
            tFieldLabels.push( "GeometryTag");
            tFieldLabels.push( "PhysicalTag");

            for ( mesh::Field * tField : aFields )
            {
                if( tField->write_field_to_file() )
                {
                    tFieldLabels.push( tField->label() );
                }
            }

            /*  initialize data container */
            mError = ex_put_variable_param(
                    mHandle,
                    EX_ELEM_BLOCK,
                    tNumElementFields );

            this->check( "ex_put_variable_param (element)");

            /* write names */
            mError = ex_put_variable_names(
                    mHandle,
                    EX_ELEM_BLOCK,
                    tNumElementFields,
                    tFieldLabels.data() );

            this->check( "ex_put_variable_names (element)");

            // create the truth table for the element fields
            uint tSize = tNumElementFields * mMesh->number_of_elements();

            int * tTruthTable = ( int * ) malloc( tSize * sizeof( int  ));

            for ( uint k = 0; k < tSize; ++k )
            {
                tTruthTable[ k ] = 1;
            }

            mError = ex_put_truth_table(
                    mHandle,
                    EX_ELEM_BLOCK,
                    mMesh->number_of_blocks(),
                    tNumElementFields,
                    tTruthTable );

            this->check( "ex_put_truth_table (element)");

            // tidy up mempry
            free( tTruthTable );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Intrinsic Fields
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for( uint k=0; k<4; ++k )
            {
                // increment id counter
                ++tID;

                // loop over all blocks
                for ( uint b = 0; b < tNumBlocks; ++b )
                {
                    // get block
                    mesh::Block * tBlock = tBlocks( b );

                    uint tNumElements = tBlock->number_of_elements();

                    // allocate data object
                    real * tBlockData = ( real * ) malloc( tNumElements * sizeof( real ) );

                    switch( k )
                    {
                        case ( 0 ):
                        {
                            for ( uint e = 0; e < tNumElements; ++e )
                            {
                                tBlockData[ e ] = ( real ) tBlock->element( e )->id();
                            }
                            break;
                        }
                        case ( 1 ):
                        {
                            for ( uint e = 0; e < tNumElements; ++e )
                            {
                                tBlockData[ e ] = ( real ) tBlock->element( e )->owner();
                            }
                            break;
                        }
                        case ( 2 ):
                        {
                            for ( uint e = 0; e < tNumElements; ++e )
                            {
                                tBlockData[ e ] = ( real ) tBlock->element( e )->geometry_tag();
                            }
                            break;
                        }
                        case ( 3 ):
                        {
                            for ( uint e = 0; e < tNumElements; ++e )
                            {
                                tBlockData[ e ] = ( real ) tBlock->element( e )->physical_tag();
                            }
                            break;
                        }
                        default:
                        {
                            break;
                        }

                    }


                    // push data
                    mError = ex_put_var(
                            mHandle,
                            mTimeStep,
                            EX_ELEM_BLOCK,
                            tID,
                            tBlock->id(),
                            tNumElements,
                            tBlockData );

                    this->check( "ex_put_var (element, intrinsic)");

                    free( tBlockData );
                }
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Other fields
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for( mesh::Field * tField : aFields )
            {
                if( tField->write_field_to_file() )
                {
                    // increment id counter
                    ++tID;

                    // get field data
                    Vector< real > tData = tField->data();

                    // loop over all blocks
                    for ( uint b = 0; b < tNumBlocks; ++b )
                    {
                        // get block
                        mesh::Block * tBlock = tBlocks( b );

                        uint tNumElements = tBlock->number_of_elements();

                        // allocate data object
                        real * tBlockData = ( real * ) malloc( tNumElements * sizeof( real ));

                        // assemble data blockwise
                        for ( uint e = 0; e < tNumElements; ++e )
                        {
                            tBlockData[ e ] = tData( tBlock->element( e )->index() );
                        }

                        // push data
                        mError = ex_put_var(
                                mHandle,
                                mTimeStep,
                                EX_ELEM_BLOCK,
                                tID,
                                tBlock->id(),
                                tNumElements,
                                tBlockData );

                        this->check( "ex_put_var (element)");
                        // tidy up memory
                        free( tBlockData );
                    }
                }
            }

#endif
        }

//------------------------------------------------------------------------------

        int
        ExodusWriter::fix_num_nodes( const ElementType & aType )
        {
            switch( aType )
            {
                case( ElementType::QUAD16 ) :
                {
                    return 4 ;
                }
                case( ElementType::TRI15 ) :
                case( ElementType::TRI21 ) :
                {
                    return 3 ;
                }
                case( ElementType::TET20 ) :
                {
                    return 4 ;
                }
                case ( ElementType::PENTA15 ) :
                case ( ElementType::PENTA18 ) :
                {
                    return 6 ;
                }
                case ( ElementType::HEX64 ) :
                {
                    return 8 ;
                }
                default:
                {
                    return number_of_nodes( aType );
                }
            }
        }

    }
}