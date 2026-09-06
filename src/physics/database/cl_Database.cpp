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

#include "commtools.hpp"
#include "cl_Database.hpp"
#include "cl_HDF5.hpp"
#include "assert.hpp"
#include "cl_Solver.hpp"
#include "cl_DatabaseProjector.hpp"

namespace belfem
{
    Database::Database( const string & aFilePath, const string & aMaterial ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aMaterial )
    {
        if ( comm_rank() == 0 )
        {
            HDF5 tFile( aFilePath, FileMode::OPEN_RDONLY );
            tFile.select_group( aMaterial );

            this->load( tFile.active_group() );
            tFile.close_active_group();
            tFile.close();
        }
        else
        {
            this->load( 0 );
        }

        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );

    }

    Database::Database( const hid_t aHDF5, const string & aLabel ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aLabel )
    {
        this->load( aHDF5 );
        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );
    }

    Database::Database( Mesh * aMesh, const string & aField, bool aProject, const string aMaterial ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mLabel( aField )
    {
        uint tDimension ;
        uint tOrder ;
        Vector< real > tOrigin ;
        Vector< real > tNodeSteps ;
        Vector< index_t > tNumNodes ;

        if ( aProject )
        {
            database::Projector tProjector( aMesh );
            if ( aMaterial != "" ) tProjector.set_label( aMaterial );
            tProjector.project( aField, mValues );
        }
        else if ( mCommRank == 0 )
        {
            mValues = aMesh->field( aField )->data();
            share( mValues );
        }
        else
        {
            receive( mValues );
        }

        if ( mCommRank == 0 )
        {
            BELFEM_ERROR( aMesh->is_tensormesh(), "Database::Database( Mesh * aMesh, const string & aField ) expects a tensor mesh" );

            tDimension = aMesh->number_of_dimensions();
            tOrder     = aMesh->tensorconf()->order();
            tOrigin    = aMesh->tensorconf()->origin();
            tNodeSteps = aMesh->tensorconf()->step();
            tNumNodes  = aMesh->tensorconf()->num_nodes_vector();

            if ( mCommSize > 1 )
            {
                broadcast( tOrder );
                broadcast( tDimension );
                broadcast( tOrigin );
                broadcast( tNumNodes );
                broadcast( tNodeSteps );
            }
        }
        else
        {
            broadcast( tOrder );
            broadcast( tDimension );
            broadcast( tOrigin );
            broadcast( tNumNodes );
            broadcast( tNodeSteps );
        }

        mConfig = new TensorMeshConfig( tOrder, tNumNodes, tNodeSteps, tOrigin  );
        this->set_element_type();
        TensorMeshFactory tFactory ;
        tFactory.create_topology( mConfig, mTopology );
    }

    Database::~Database()
    {
        if ( mN != nullptr )
        {
            free( mN );
        }
        if ( mConfig != nullptr )
        {
            delete mConfig ;
        }
    }

    void Database::save( const hid_t aHDF5 )
    {
        if ( mCommRank != 0 ) return;

        hid_t tID = aHDF5 ;
        herr_t tStatus ;
        hdf5::save_scalar_to_file( tID, "order", mConfig->order(), tStatus );
        hdf5::save_vector_to_file( tID, "origin", mConfig->origin(), tStatus );
        hdf5::save_vector_to_file( tID, "points", mConfig->num_nodes_vector(), tStatus );
        hdf5::save_vector_to_file( tID, "step", mConfig->step(), tStatus );
        hdf5::save_vector_to_file( tID, "values", mValues, tStatus );
    }

    void
    Database::load( const hid_t aHDF5 )
    {
        uint tDimension ;
        uint32_t tOrder ;

        Vector< uint32_t > tNumPoints ;
        Vector< real >     tNodeSteps ;
        Vector< real >     tOrigin ;

        if ( mCommRank == 0 )
        {
            hid_t tID = aHDF5 ;
            herr_t tStatus ;
            hdf5::load_scalar_from_file( tID, "order", tOrder, tStatus );

            hdf5::load_vector_from_file( tID, "origin", tOrigin, tStatus );
            tDimension = tOrigin.length() ;

            hdf5::load_vector_from_file( tID, "points", tNumPoints, tStatus );
            hdf5::load_vector_from_file( tID, "step", tNodeSteps, tStatus );
            hdf5::load_vector_from_file( tID, "values", mValues, tStatus );

            if ( mCommSize > 1 )
            {
                broadcast( tOrder );
                broadcast( tDimension );
                broadcast( tOrigin );
                broadcast( tNumPoints );
                broadcast( tNodeSteps );
                share( mValues );
            }
        }
        else
        {
            broadcast( tOrder );
            broadcast( tDimension );
            broadcast( tOrigin );
            broadcast( tNumPoints );
            broadcast( tNodeSteps );
            receive( mValues );
        }

        mConfig = new TensorMeshConfig( tOrder, tNumPoints, tNodeSteps, tOrigin  );
    }

    void
    Database::set_element_type()
    {
        switch ( mConfig->element_type() )
        {
            case ElementType::QUAD4:
            {
                mFunction2D      = & Database::eval_quad4 ;
                mdFunction2Ddxi  = & Database::deval_quad4dxi ;
                mdFunction2Ddeta = & Database::deval_quad4deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::QUAD9:
            {
                mFunction2D      = & Database::eval_quad9 ;
                mdFunction2Ddxi  = & Database::deval_quad9dxi ;
                mdFunction2Ddeta = & Database::deval_quad9deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::QUAD16:
            {
                mFunction2D      = & Database::eval_quad16 ;
                mdFunction2Ddxi  = & Database::deval_quad16dxi ;
                mdFunction2Ddeta = & Database::deval_quad16deta ;
                mFunction3D      = nullptr ;
                mdFunction3Ddxi  = nullptr ;
                break ;
            }
            case ElementType::HEX8:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex8 ;
                mdFunction3Ddxi   = & Database::deval_hex8dxi ;
                mdFunction3Ddeta  = & Database::deval_hex8deta ;
                mdFunction3Ddzeta = & Database::deval_hex8dzeta ;
                break ;
            }
            case ElementType::HEX27:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex27 ;
                mdFunction3Ddxi   = & Database::deval_hex27dxi ;
                mdFunction3Ddeta  = & Database::deval_hex27deta ;
                mdFunction3Ddzeta = & Database::deval_hex27dzeta ;
                break ;
            }
            case ElementType::HEX64:
            {
                mFunction2D       = nullptr ;
                mdFunction2Ddxi   = nullptr ;
                mFunction3D       = & Database::eval_hex64 ;
                mdFunction3Ddxi   = & Database::deval_hex64dxi ;
                mdFunction3Ddeta  = & Database::deval_hex64deta ;
                mdFunction3Ddzeta = & Database::deval_hex64dzeta ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "unsupported element type in database %s", mLabel.c_str() );
            }
        }

        mNumNodesPerElement = mesh::number_of_nodes( mConfig->element_type() );

        if ( mN != nullptr ) free( mN );
        mN = ( real * ) malloc( mNumNodesPerElement * sizeof( real ) );
    }

    real
    Database::min( const uint aDimension ) const
    {
        return mConfig->min( aDimension );
    }

    real
    Database::max( const uint aDimension ) const
    {
        return mConfig->max( aDimension );
    }
}
