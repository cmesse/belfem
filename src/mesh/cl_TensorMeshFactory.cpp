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

#include "cl_TensorMeshFactory.hpp"
#include "assert.hpp"
#include "fn_linspace.hpp"
#include "cl_Element_Factory.hpp"
#include "commtools.hpp"

namespace belfem
{
    void
    TensorMeshFactory::populate_tensor_mesh(
        const TensorMeshConfig * aConfig,
        Mesh * aMesh,
        const proc_t aMasterProc )
    {


        // make sure that mesh is empty
        BELFEM_ERROR( aMesh->number_of_nodes() == 0, "Mesh is not empty" );
        BELFEM_ERROR( aMesh->number_of_elements() == 0, "Mesh is not empty" );
        BELFEM_ERROR( aMesh->number_of_blocks() == 0, "Mesh is not empty" );
        BELFEM_ERROR( aMesh->number_of_sidesets() == 0, "Mesh is not empty" );
        BELFEM_ERROR( aMesh->number_of_dimensions() == aConfig->dimension(),
            "Mesh dimension does not match config" );

        if( comm_rank() == aMasterProc )
        {
            mConfig = aConfig ;

            Matrix< real > tGrid ;
            this->create_grid(
                aConfig->num_nodes_vector(),
                aConfig->origin(),
                aConfig->limit(), tGrid );

            Matrix< index_t > tAdjacency ;
            this->create_topology( aConfig,
                tAdjacency );

            Cell< mesh::Node * > & tNodes = aMesh->nodes();
            index_t tNumNodes = aConfig->num_nodes();
            tNodes.set_size( tNumNodes, nullptr );

            if ( aConfig->dimension() == 2 )
            {
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    tNodes( k ) = new mesh::Node( k+1, tGrid( 0, k ), tGrid( 1, k ) );
                }
            }
            else if ( aConfig->dimension() == 3 )
            {
                for ( index_t k=0; k<tNumNodes; ++k )
                {
                    tNodes( k ) = new mesh::Node( k+1, tGrid( 0, k ), tGrid( 1, k ), tGrid( 2, k ) );
                }
            }

            index_t tNumElems = aConfig->num_elements();

            mesh::Block * tBlock = new mesh::Block( 1, tNumElems ) ;
            aMesh->blocks().set_size( 1, tBlock );
            mesh::ElementFactory tFactory ;

            index_t tNumNodesPerElem = tAdjacency.n_rows();
            for ( index_t e=0; e<tNumElems; ++e )
            {
                mesh::Element * tElement = tFactory.create_element( aConfig->element_type(), e+1 );
                for ( index_t k=0; k<tNumNodesPerElem; ++k )
                {
                    tElement->insert_node( tNodes( tAdjacency( k, e ) ), k );
                }
                tBlock->insert_element( tElement );
            }
            aMesh->finalize() ;

            mConfig = nullptr ;
        }
    }


//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_topology( const TensorMeshConfig * aConfig, Matrix< index_t > &aTopology )
    {
        mConfig = aConfig ;
        Vector< uint > aNumElems = mConfig->num_elements_vector() ;
        switch ( mConfig->element_type() )
        {
        case ElementType::QUAD4 :
        {
            this->create_topology_quad4( aNumElems, aTopology );
            break ;
        }
        case ElementType::QUAD9 :
        {
            this->create_topology_quad9( aNumElems, aTopology );
            break ;
        }
        case ElementType::QUAD16 :
        {
            this->create_topology_quad16( aNumElems, aTopology );
            break ;
        }
        case ElementType::HEX8 :
        {
            this->create_topology_hex8( aNumElems, aTopology );
            break ;
        }
        case ElementType::HEX27 :
        {
            this->create_topology_hex27( aNumElems, aTopology );
            break ;
        }
        case ElementType::HEX64 :
        {
            this->create_topology_hex64( aNumElems, aTopology );
            break;
        }
        default:
        {
            BELFEM_ERROR( false, "Invalid element type" );
        }
        }
    }

    void
    TensorMeshFactory::create_grid(
            const Vector< index_t > & aNumNodes,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
            Matrix< real > & aGrid ) const
    {
        uint tNumDim = aNumNodes.length();

        BELFEM_ASSERT( aMinPoint.length() == tNumDim,
            "dimension of min point does not match");
        BELFEM_ASSERT( aMaxPoint.length() == tNumDim,
            "dimension of max point does not match");

        BELFEM_ASSERT( tNumDim == 2 || tNumDim == 3,
                      "dimension must be 2 or 3");

        Vector< real > tX ;
        Vector< real > tY ;


        linspace( aMinPoint( 0 ), aMaxPoint( 0 ), aNumNodes( 0 ), tX );

        linspace( aMinPoint( 1 ), aMaxPoint( 1 ), aNumNodes( 1 ), tY );

        index_t nx = aNumNodes( 0 );
        index_t ny = aNumNodes( 1 );

        index_t tCount = 0 ;

        if ( tNumDim == 3 )
        {
            Vector< real > tZ ;
            linspace( aMinPoint( 2 ), aMaxPoint( 2 ), aNumNodes( 2 ), tZ );
            index_t nz = aNumNodes( 2 );

            aGrid.set_size( 3, nx * ny * nz );
            for ( index_t i=0; i<nx; ++i )
            {
                for ( index_t j=0; j<ny; ++j )
                {
                    for ( index_t k=0; k<nz; ++k )
                    {
                        aGrid( 0, tCount ) = tX( i );
                        aGrid( 1, tCount ) = tY( j );
                        aGrid( 2, tCount ) = tZ( k );
                        ++tCount ;
                    }
                }
            }
        }
        else
        {
            aGrid.set_size( 2, nx * ny );
            for ( index_t i=0; i<nx; ++i )
            {
                for ( index_t j=0; j<ny; ++j )
                {
                    aGrid( 0, tCount ) = tX( i );
                    aGrid( 1, tCount ) = tY( j );
                    ++tCount ;
                }
            }
        }

    }


    void
    TensorMeshFactory::create_bsplines( Mesh * aMesh )
    {
        BELFEM_ERROR( aMesh->is_tensormesh(), "Mesh must be a tensor mesh" );

        // only the master rank populates a tensor mesh ( see
        // populate_tensor_mesh ) -- a rank holding an empty copy has nothing
        // to link, and the loops below iterate the config-sized element grid,
        // not the actual container
        if ( aMesh->number_of_elements() == 0 )
        {
            return;
        }

        mConfig = aMesh->tensorconf();

        switch (  aMesh->number_of_dimensions() )
        {
            case 2 :
            {
                switch ( aMesh->tensorconf()->order() )
                {
                    case 2 :
                    {
                        this->create_bsplines_quad9(
                            aMesh->tensorconf(), aMesh->control_points() );
                        this->link_bsplines_quad9(
                            aMesh->tensorconf(), aMesh->control_points(), aMesh->elements() );
                        this->compute_tmatrix_quad9();
                        break ;
                    }
                    case 3 :
                    {
                        this->create_bsplines_quad16(
                            aMesh->tensorconf(), aMesh->control_points() );
                        this->link_bsplines_quad16(
                            aMesh->tensorconf(), aMesh->control_points(), aMesh->elements() );
                        this->compute_tmatrix_quad16();
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "unsupported order" );
                    }
                }
                break ;
            }
            case 3 :
            {
                switch ( aMesh->tensorconf()->order() )
                {
                    case 2 :
                    {
                        this->create_bsplines_hex27(
                            aMesh->tensorconf(), aMesh->control_points() );
                        this->link_bsplines_hex27(
                            aMesh->tensorconf(), aMesh->control_points(), aMesh->elements() );
                        this->compute_tmatrix_hex27();
                        break ;
                    }
                    case 3 :
                    {
                        this->create_bsplines_hex64(
                            aMesh->tensorconf(), aMesh->control_points() );
                        this->link_bsplines_hex64(
                            aMesh->tensorconf(), aMesh->control_points(), aMesh->elements() );
                        this->compute_tmatrix_hex64();
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "unsupported order" );
                    }
                }

                break ;
            }
            default :
            {
                BELFEM_ERROR( false, "unsupported dimension" );
            }
        }
        this->connect_nodes_to_bsplines( aMesh );
        mConfig = nullptr ;
    }

    void
    TensorMeshFactory::create_topology_quad4( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1 ) ;

        aTopology.set_size( 4, ni * nj );

        index_t c = 0 ;
        for ( index_t i=0; i<ni; ++i )
        {
            for ( index_t j=0; j<nj; ++j )
            {
                aTopology( 0, c ) = this->nidx( i, j );
                aTopology( 1, c ) = this->nidx( i+1, j );
                aTopology( 2, c ) = this->nidx( i+1, j+1 );
                aTopology( 3, c ) = this->nidx( i, j+1 );
                ++c ;
            }
        }
    }

    void
    TensorMeshFactory::create_topology_quad9( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1 ) ;

        aTopology.set_size( 9, ni * nj );

        index_t c = 0 ;
        index_t p = 1 ;
        index_t q = 1 ;

        for ( index_t i=0; i<ni; ++i )
        {
            q = 1 ;
            for ( index_t j=0; j<nj; ++j )
            {
                aTopology( 0, c ) = this->nidx( p-1,q-1);
                aTopology( 1, c ) = this->nidx( p+1,q-1);
                aTopology( 2, c ) = this->nidx( p+1,q+1);
                aTopology( 3, c ) = this->nidx( p-1,q+1);
                aTopology( 4, c ) = this->nidx( p,q-1);
                aTopology( 5, c ) = this->nidx( p+1,q);
                aTopology( 6, c ) = this->nidx( p,q+1);
                aTopology( 7, c ) = this->nidx( p-1,q);
                aTopology( 8, c ) = this->nidx( p,q);
                ++c ;
                q += 2 ;
            }
            p += 2 ;
        }
    }

    void
    TensorMeshFactory::create_topology_quad16( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1 ) ;

        aTopology.set_size( 16, ni * nj );

        index_t c = 0 ;

        index_t p = 0 ;
        index_t q = 0 ;

        for ( index_t i=0; i<ni; ++i )
        {
            q = 0 ;
            for ( index_t j=0; j<nj; ++j )
            {
                aTopology(  0, c ) = this->nidx( p, q );
                aTopology(  1, c ) = this->nidx( p+3, q );
                aTopology(  2, c ) = this->nidx( p+3, q+3 );
                aTopology(  3, c ) = this->nidx( p, q+3 );

                aTopology(  4, c ) = this->nidx( p+1, q );
                aTopology(  5, c ) = this->nidx( p+2, q );
                aTopology(  6, c ) = this->nidx( p+3, q+1 );
                aTopology(  7, c ) = this->nidx( p+3, q+2 );
                aTopology(  8, c ) = this->nidx( p+2, q+3 );

                aTopology(  9, c ) = this->nidx( p+1, q+3 );
                aTopology( 10, c ) = this->nidx( p, q+2 );
                aTopology( 11, c ) = this->nidx( p, q+1 );
                aTopology( 12, c ) = this->nidx( p+1, q+1 );
                aTopology( 13, c ) = this->nidx( p+2, q+1 );

                aTopology( 14, c ) = this->nidx( p+2, q+2 );
                aTopology( 15, c ) = this->nidx( p+1, q+2 );

                q += 3 ;
                ++c ;
            }

            p += 3 ;
        }
    }

    void
    TensorMeshFactory::create_topology_hex8( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1 ) ;
        index_t nk = aNumElems( 2 );

        aTopology.set_size( 8, ni*nj*nk );

        index_t c = 0 ;

        for ( index_t i=0; i<ni; ++i )
        {
            for ( index_t j=0; j<nj; ++j )
            {
                for ( index_t k=0; k<nk; ++k )
                {
                    aTopology( 0, c ) = this->nidx( i, j, k );
                    aTopology( 1, c ) = this->nidx( i+1, j, k );
                    aTopology( 2, c ) = this->nidx( i+1, j+1, k );
                    aTopology( 3, c ) = this->nidx( i, j+1, k );
                    aTopology( 4, c ) = this->nidx( i, j, k+1 );
                    aTopology( 5, c ) = this->nidx( i+1, j, k+1 );
                    aTopology( 6, c ) = this->nidx( i+1, j+1, k+1 );
                    aTopology( 7, c ) = this->nidx( i, j+1, k+1 );
                    ++c ;
                }
            }
        }
    }

    void
    TensorMeshFactory::create_topology_hex27( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1) ;
        index_t nk = aNumElems( 2 );

        aTopology.set_size( 27, ni * nj * nk );

        index_t c = 0 ;

        index_t p = 1 ;
        index_t q = 1 ;
        index_t r = 1 ;

        for ( index_t i=0; i<ni; ++i )
        {
            q = 1 ;
            for ( index_t j=0; j<nj; ++j )
            {
                r = 1 ;
                for ( index_t k=0; k<nk; ++k )
                {
                    aTopology(0, c ) = this->nidx(p - 1, q - 1, r - 1);
                    aTopology(1, c ) = this->nidx(p + 1, q - 1, r - 1);
                    aTopology(2, c ) = this->nidx(p + 1, q + 1, r - 1);
                    aTopology(3, c ) = this->nidx(p - 1, q + 1, r - 1);

                    aTopology(4, c ) = this->nidx(p - 1, q - 1, r + 1);
                    aTopology(5, c ) = this->nidx(p + 1, q - 1, r + 1);
                    aTopology(6, c ) = this->nidx(p + 1, q + 1, r + 1);
                    aTopology(7, c ) = this->nidx(p - 1, q + 1, r + 1);

                    aTopology(8, c ) = this->nidx(p, q - 1, r - 1);
                    aTopology(9, c ) = this->nidx(p + 1, q, r - 1);
                    aTopology(10, c) = this->nidx(p, q + 1, r - 1);
                    aTopology(11, c) = this->nidx(p - 1, q, r - 1);

                    aTopology(12, c) = this->nidx(p - 1, q - 1, r);
                    aTopology(13, c) = this->nidx(p + 1, q - 1, r);
                    aTopology(14, c) = this->nidx(p + 1, q + 1, r);
                    aTopology(15, c) = this->nidx(p - 1, q + 1, r);

                    aTopology(16, c) = this->nidx(p,q-1,r+1);
                    aTopology(17, c) = this->nidx(p+1,q,r+1);
                    aTopology(18, c) = this->nidx(p, q + 1, r + 1);
                    aTopology(19, c) = this->nidx(p-1, q, r + 1);

                    aTopology(20, c) = this->nidx(p, q, r);

                    aTopology(21, c) = this->nidx(p, q, r - 1);
                    aTopology(22, c) = this->nidx(p, q, r + 1);

                    aTopology(23, c) = this->nidx(p-1, q, r);
                    aTopology(24, c) = this->nidx(p + 1, q, r);

                    aTopology(25, c) = this->nidx(p, q-1, r);
                    aTopology(26, c) = this->nidx(p, q+1, r);

                    ++c ;
                    r+=2 ;
                }
                q+=2 ;
            }
            p += 2 ;
        }
    }

    void
    TensorMeshFactory::create_topology_hex64( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const
    {
        index_t ni = aNumElems( 0 );
        index_t nj = aNumElems( 1 ) ;
        index_t nk =  aNumElems( 2 );

        aTopology.set_size( 64, ni * nj * nk );

        index_t c = 0 ;

        index_t p = 0 ;
        index_t q = 0 ;
        index_t r = 0 ;
        for ( index_t i=0; i<ni; ++i )
        {
            q = 0 ;
            for ( index_t j=0; j<nj; ++j )
            {
                r = 0 ;
                for ( index_t k=0; k<nk; ++k )
                {
                    aTopology(0, c) = this->nidx(p, q, r);
                    aTopology(1, c) = this->nidx(p + 3, q, r);
                    aTopology(2, c) = this->nidx(p + 3, q + 3, r);
                    aTopology(3, c) = this->nidx(p, q + 3, r);
                    aTopology(4, c) = this->nidx(p, q, r + 3);
                    aTopology(5, c) = this->nidx(p + 3, q, r + 3);
                    aTopology(6, c) = this->nidx(p + 3, q + 3, r + 3);
                    aTopology(7, c) = this->nidx(p, q + 3, r + 3);
                    aTopology(8, c) = this->nidx(p + 1, q, r);
                    aTopology(9, c) = this->nidx(p + 2, q, r);
                    aTopology(10, c) = this->nidx(p, q + 1, r);
                    aTopology(11, c) = this->nidx(p, q + 2, r);
                    aTopology(12, c) = this->nidx(p, q, r + 1);
                    aTopology(13, c) = this->nidx(p, q, r + 2);
                    aTopology(14, c) = this->nidx(p + 3, q + 1, r);
                    aTopology(15, c) = this->nidx(p + 3, q + 2, r);
                    aTopology(16, c) = this->nidx(p + 3, q, r + 1);
                    aTopology(17, c) = this->nidx(p + 3, q, r + 2);
                    aTopology(18, c) = this->nidx(p + 2, q + 3, r);
                    aTopology(19, c) = this->nidx(p + 1, q + 3, r);
                    aTopology(20, c) = this->nidx(p + 3, q + 3, r + 1);
                    aTopology(21, c) = this->nidx(p + 3, q + 3, r + 2);
                    aTopology(22, c) = this->nidx(p, q + 3, r + 1);
                    aTopology(23, c) = this->nidx(p, q + 3, r + 2);
                    aTopology(24, c) = this->nidx(p + 1, q, r + 3);
                    aTopology(25, c) = this->nidx(p + 2, q, r + 3);
                    aTopology(26, c) = this->nidx(p, q + 1, r + 3);
                    aTopology(27, c) = this->nidx(p, q + 2, r + 3);
                    aTopology(28, c) = this->nidx(p + 3, q + 1, r + 3);
                    aTopology(29, c) = this->nidx(p + 3, q + 2, r + 3);
                    aTopology(30, c) = this->nidx(p + 2, q + 3, r + 3);
                    aTopology(31, c) = this->nidx(p + 1, q + 3, r + 3);
                    aTopology(32, c) = this->nidx(p + 1, q + 1, r);
                    aTopology(33, c) = this->nidx(p + 1, q + 2, r);
                    aTopology(34, c) = this->nidx(p + 2, q + 2, r);
                    aTopology(35, c) = this->nidx(p + 2, q + 1, r);
                    aTopology(36, c) = this->nidx(p + 1, q, r + 1);
                    aTopology(37, c) = this->nidx(p + 2, q, r + 1);
                    aTopology(38, c) = this->nidx(p + 2, q, r + 2);
                    aTopology(39, c) = this->nidx(p + 1, q, r + 2);
                    aTopology(40, c) = this->nidx(p, q + 1, r + 1);
                    aTopology(41, c) = this->nidx(p, q + 1, r + 2);
                    aTopology(42, c) = this->nidx(p, q + 2, r + 2);
                    aTopology(43, c) = this->nidx(p, q + 2, r + 1);
                    aTopology(44, c) = this->nidx(p + 3, q + 1, r + 1);
                    aTopology(45, c) = this->nidx(p + 3, q + 2, r + 1);
                    aTopology(46, c) = this->nidx(p + 3, q + 2, r + 2);
                    aTopology(47, c) = this->nidx(p + 3, q + 1, r + 2);
                    aTopology(48, c) = this->nidx(p + 2, q + 3, r + 1);
                    aTopology(49, c) = this->nidx(p + 1, q + 3, r + 1);
                    aTopology(50, c) = this->nidx(p + 1, q + 3, r + 2);
                    aTopology(51, c) = this->nidx(p + 2, q + 3, r + 2);
                    aTopology(52, c) = this->nidx(p + 1, q + 1, r + 3);
                    aTopology(53, c) = this->nidx(p + 2, q + 1, r + 3);
                    aTopology(54, c) = this->nidx(p + 2, q + 2, r + 3);
                    aTopology(55, c) = this->nidx(p + 1, q + 2, r + 3);
                    aTopology(56, c) = this->nidx(p + 1, q + 1, r + 1);
                    aTopology(57, c) = this->nidx(p + 2, q + 1, r + 1);
                    aTopology(58, c) = this->nidx(p + 2, q + 2, r + 1);
                    aTopology(59, c) = this->nidx(p + 1, q + 2, r + 1);
                    aTopology(60, c) = this->nidx(p + 1, q + 1, r + 2);
                    aTopology(61, c) = this->nidx(p + 2, q + 1, r + 2);
                    aTopology(62, c) = this->nidx(p + 2, q + 2, r + 2);
                    aTopology(63, c) = this->nidx(p + 1, q + 2, r + 2);

                    ++c ;

                    r+=3 ;
                }
                q +=3 ;
            }
            p+= 3 ;
        }
    }

//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_bsplines_quad9(
        const TensorMeshConfig * aConfig,
        Cell< mesh::ControlPoint * > & aPoints )
    {
        aPoints.set_size( aConfig->num_control_points(), nullptr );

        index_t tCount = 0 ;

        index_t ni = aConfig->num_control_points( 0 ) ;
        index_t nj = aConfig->num_control_points( 1 ) ;

        real dx = aConfig->step()(0);
        real dy = aConfig->step()(1);

        real x0 = aConfig->origin()(0);
        real y0 = aConfig->origin()(1);


        real x = x0- 0.5 * dx ;
        for ( index_t i=0; i<ni; ++i, x += dx )
        {
            real y = y0 - 0.5 * dy ;
            for ( index_t j=0; j<nj; ++j, y += dy )
            {
                aPoints( tCount ) = new mesh::ControlPoint( tCount+1, x, y );
                ++tCount ;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_bsplines_quad16(
        const TensorMeshConfig * aConfig,
        Cell< mesh::ControlPoint * > & aPoints )
    {
        aPoints.set_size( aConfig->num_control_points(), nullptr );

        index_t ni = aConfig->num_control_points( 0 ) ;
        index_t nj = aConfig->num_control_points( 1 ) ;

        real dx = aConfig->step()(0);
        real dy = aConfig->step()(1);

        real x0 = aConfig->origin()(0);
        real y0 = aConfig->origin()(1);

        index_t tCount = 0 ;

        real x = x0 - dx ;
        for ( index_t i=0; i<ni; ++i, x += dx )
        {
            real y = y0 - dy ;
            for ( index_t j=0; j<nj; ++j, y += dy )
            {
                aPoints( tCount ) = new mesh::ControlPoint( tCount+1, x, y );
                ++tCount ;
            }
        }
    }

    void
    TensorMeshFactory::create_bsplines_hex27(
        const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints )
    {
        aPoints.set_size( aConfig->num_control_points(), nullptr );

        index_t ni = aConfig->num_control_points( 0 ) ;
        index_t nj = aConfig->num_control_points( 1 ) ;
        index_t nk = aConfig->num_control_points( 2 ) ;

        real dx = aConfig->step()(0);
        real dy = aConfig->step()(1);
        real dz = aConfig->step()(2);

        real x0 = aConfig->origin()(0);
        real y0 = aConfig->origin()(1);
        real z0 = aConfig->origin()(2);

        index_t tCount = 0 ;

        real x = x0 - 0.5 * dx ;
        for ( index_t i=0; i<ni; ++i, x += dx )
        {
            real y = y0 - 0.5 * dy ;
            for ( index_t j=0; j<nj; ++j, y+= dy )
            {
                real z = z0 - 0.5 * dz ;
                for ( index_t k=0; k<nk; ++k, z+=dz )
                {
                    aPoints( tCount ) = new mesh::ControlPoint( tCount+1, x, y, z );
                    ++tCount ;
                }
            }

        }
    }

    void
    TensorMeshFactory::create_bsplines_hex64(
        const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints )
    {
        aPoints.set_size( aConfig->num_control_points(), nullptr );

        index_t ni = aConfig->num_control_points( 0 ) ;
        index_t nj = aConfig->num_control_points( 1 ) ;
        index_t nk = aConfig->num_control_points( 2 ) ;

        real dx = aConfig->step()(0);
        real dy = aConfig->step()(1);
        real dz = aConfig->step()(2);

        real x0 = aConfig->origin()(0);
        real y0 = aConfig->origin()(1);
        real z0 = aConfig->origin()(2);

        index_t tCount = 0 ;

        real x = x0 - dx ;
        for ( index_t i=0; i<ni; ++i, x+= dx )
        {
            real y = y0 - dy ;
            for ( index_t j=0; j<nj; ++j, y+= dy )
            {
                real z = z0 - dz ;
                for ( index_t k=0; k<nk; ++k, z+=dz )
                {
                    aPoints( tCount ) = new mesh::ControlPoint( tCount+1, x, y, z );
                    ++tCount ;
                }
            }
        }
    }

    void
    TensorMeshFactory::link_bsplines_quad9(
        const TensorMeshConfig * aConfig,
        Cell< mesh::ControlPoint *  > & aPoints,
        Cell< mesh::Element * >       & aElements )
    {
        index_t ni = aConfig->num_elements( 0 ) ;
        index_t nj = aConfig->num_elements( 1 ) ;

        index_t tCount = 0 ;

        for ( index_t p=1; p<=ni; ++p )
        {
            for ( index_t q=1; q<=nj; ++q )
            {
                mesh::Element * tElement = aElements( tCount++ );

                tElement->allocate_control_points_container( 9 );

                tElement->insert_control_point( aPoints( bidx(  p-1, q-1 ) ), 0 );
                tElement->insert_control_point( aPoints( bidx(  p+1, q-1 ) ), 1 );
                tElement->insert_control_point( aPoints( bidx(  p+1, q+1 ) ), 2 );
                tElement->insert_control_point( aPoints( bidx(  p-1, q+1 ) ), 3 );

                tElement->insert_control_point( aPoints( bidx(  p  , q-1 ) ), 4 );
                tElement->insert_control_point( aPoints( bidx(  p+1, q   ) ), 5 );
                tElement->insert_control_point( aPoints( bidx(  p  , q+1 ) ), 6 );
                tElement->insert_control_point( aPoints( bidx(  p-1, q   ) ), 7 );
                tElement->insert_control_point( aPoints( bidx(  p  , q   ) ), 8 );
            }
        }
    }

    void
    TensorMeshFactory::link_bsplines_quad16(
       const TensorMeshConfig * aConfig,
       Cell< mesh::ControlPoint *  > & aPoints,
       Cell< mesh::Element * >       & aElements )
    {
        index_t ni = aConfig->num_elements( 0 );
        index_t nj = aConfig->num_elements( 1 );

        index_t tCount = 0 ;

        for ( index_t p=0; p<ni; ++p )
        {
            for ( index_t q=0; q<nj; ++q )
            {
                mesh::Element * tElement = aElements( tCount++ );

                tElement->allocate_control_points_container( 16 );

                tElement->insert_control_point( aPoints( bidx(  p  , q   ) ), 0 );
                tElement->insert_control_point( aPoints( bidx(  p+3, q   ) ), 1 );
                tElement->insert_control_point( aPoints( bidx(  p+3, q+3 ) ), 2 );
                tElement->insert_control_point( aPoints( bidx(  p  , q+3 ) ), 3 );

                tElement->insert_control_point( aPoints( bidx(  p  , q+1 ) ), 4 );
                tElement->insert_control_point( aPoints( bidx(  p  , q+2 ) ), 5 );

                tElement->insert_control_point( aPoints( bidx(  p+3, q+1 ) ), 6 );
                tElement->insert_control_point( aPoints( bidx(  p+3, q+2 ) ), 7 );

                tElement->insert_control_point( aPoints( bidx(  p+2, q+3 ) ), 8 );
                tElement->insert_control_point( aPoints( bidx(  p+1, q+3 ) ), 9 );

                tElement->insert_control_point( aPoints( bidx(  p  , q+2 ) ), 10 );
                tElement->insert_control_point( aPoints( bidx(  p  , q+1 ) ), 11 );

                tElement->insert_control_point( aPoints( bidx(  p+1, q+1 ) ), 12 );
                tElement->insert_control_point( aPoints( bidx(  p+2, q+1 ) ), 13 );

                tElement->insert_control_point( aPoints( bidx(  p+2, q+2 ) ), 14 );
                tElement->insert_control_point( aPoints( bidx(  p+1, q+2 ) ), 15 );
            }
        }
    }

    void
    TensorMeshFactory::link_bsplines_hex27(
        const TensorMeshConfig * aConfig,
        Cell< mesh::ControlPoint *  > & aPoints,
        Cell< mesh::Element * > & aElements )
    {
        index_t ni = aConfig->num_elements( 0 );
        index_t nj = aConfig->num_elements( 1 );
        index_t nk = aConfig->num_elements( 2 );
        index_t tCount = 0 ;

        for ( index_t p=1; p<=ni; ++p )
        {
            for ( index_t q=1; q<=nj; ++q )
            {
                for ( index_t r=1; r<=nk; ++r )
                {
                    mesh::Element * tElement = aElements( tCount++ );

                    tElement->allocate_control_points_container( 27 );

                    tElement->insert_control_point( aPoints( bidx(p-1,q-1,r-1) ), 0);
                    tElement->insert_control_point( aPoints( bidx(p+1,q-1,r-1) ), 1);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r-1) ), 2);
                    tElement->insert_control_point( aPoints( bidx(p-1,q+1,r-1) ), 3);
                    tElement->insert_control_point( aPoints( bidx(p-1,q-1,r+1) ), 4);
                    tElement->insert_control_point( aPoints( bidx(p+1,q-1,r+1) ), 5);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r+1) ), 6);
                    tElement->insert_control_point( aPoints( bidx(p-1,q+1,r+1) ), 7);
                    tElement->insert_control_point( aPoints( bidx(p,q-1,r-1) ), 8);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r-1) ), 9);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r-1) ), 10);
                    tElement->insert_control_point( aPoints( bidx(p-1,q,r-1) ), 11);
                    tElement->insert_control_point( aPoints( bidx(p-1,q-1,r) ), 12);
                    tElement->insert_control_point( aPoints( bidx(p+1,q-1,r) ), 13);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r) ), 14);
                    tElement->insert_control_point( aPoints( bidx(p-1,q+1,r) ), 15);
                    tElement->insert_control_point( aPoints( bidx(p,q-1,r+1) ), 16);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r+1) ), 17);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r+1) ), 18);
                    tElement->insert_control_point( aPoints( bidx(p-1,q,r+1) ), 19);
                    tElement->insert_control_point( aPoints( bidx(p,q,r) ), 20);
                    tElement->insert_control_point( aPoints( bidx(p,q,r-1) ), 21);
                    tElement->insert_control_point( aPoints( bidx(p,q,r+1) ), 22);
                    tElement->insert_control_point( aPoints( bidx(p-1,q,r) ), 23);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r) ), 24);
                    tElement->insert_control_point( aPoints( bidx(p,q-1,r) ), 25);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r) ), 26);

                }
            }
        }
    }


    void
    TensorMeshFactory::link_bsplines_hex64(
        const TensorMeshConfig * aConfig,
        Cell< mesh::ControlPoint *  > & aPoints,
        Cell< mesh::Element * > & aElements )
    {
        index_t ni = aConfig->num_elements( 0 );
        index_t nj = aConfig->num_elements( 1 );
        index_t nk = aConfig->num_elements( 2 );
        index_t tCount = 0 ;

        for ( index_t p=0; p<ni; ++p )
        {
            for ( index_t q=0; q<nj; ++q )
            {
                for ( index_t r=0; r<nk; ++r )
                {
                    mesh::Element * tElement = aElements( tCount++ );
                    tElement->allocate_control_points_container( 64 );

                    tElement->insert_control_point( aPoints( bidx(p,q,r) ),0);
                    tElement->insert_control_point( aPoints( bidx(p+3,q,r) ),1);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+3,r) ),2);
                    tElement->insert_control_point( aPoints( bidx(p,q+3,r) ),3);
                    tElement->insert_control_point( aPoints( bidx(p,q,r+3) ),4);
                    tElement->insert_control_point( aPoints( bidx(p+3,q,r+3) ),5);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+3,r+3) ),6);
                    tElement->insert_control_point( aPoints( bidx(p,q+3,r+3) ),7);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r) ),8);
                    tElement->insert_control_point( aPoints( bidx(p+2,q,r) ),9);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r) ),10);
                    tElement->insert_control_point( aPoints( bidx(p,q+2,r) ),11);
                    tElement->insert_control_point( aPoints( bidx(p,q,r+1) ),12);
                    tElement->insert_control_point( aPoints( bidx(p,q,r+2) ),13);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+1,r) ),14);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+2,r) ),15);
                    tElement->insert_control_point( aPoints( bidx(p+3,q,r+1) ),16);
                    tElement->insert_control_point( aPoints( bidx(p+3,q,r+2) ),17);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+3,r) ),18);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+3,r) ),19);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+3,r+1) ),20);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+3,r+2) ),21);
                    tElement->insert_control_point( aPoints( bidx(p,q+3,r+1) ),22);
                    tElement->insert_control_point( aPoints( bidx(p,q+3,r+2) ),23);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r+3) ),24);
                    tElement->insert_control_point( aPoints( bidx(p+2,q,r+3) ),25);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r+3) ),26);
                    tElement->insert_control_point( aPoints( bidx(p,q+2,r+3) ),27);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+1,r+3) ),28);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+2,r+3) ),29);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+3,r+3) ),30);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+3,r+3) ),31);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r) ),32);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+2,r) ),33);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+2,r) ),34);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+1,r) ),35);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r+1) ),36);
                    tElement->insert_control_point( aPoints( bidx(p+2,q,r+1) ),37);
                    tElement->insert_control_point( aPoints( bidx(p+2,q,r+2) ),38);
                    tElement->insert_control_point( aPoints( bidx(p+1,q,r+2) ),39);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r+1) ),40);
                    tElement->insert_control_point( aPoints( bidx(p,q+1,r+2) ),41);
                    tElement->insert_control_point( aPoints( bidx(p,q+2,r+2) ),42);
                    tElement->insert_control_point( aPoints( bidx(p,q+2,r+1) ),43);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+1,r+1) ),44);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+2,r+1) ),45);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+2,r+2) ),46);
                    tElement->insert_control_point( aPoints( bidx(p+3,q+1,r+2) ),47);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+3,r+1) ),48);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+3,r+1) ),49);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+3,r+2) ),50);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+3,r+2) ),51);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r+3) ),52);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+1,r+3) ),53);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+2,r+3) ),54);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+2,r+3) ),55);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r+1) ),56);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+1,r+1) ),57);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+2,r+1) ),58);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+2,r+1) ),59);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+1,r+2) ),60);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+1,r+2) ),61);
                    tElement->insert_control_point( aPoints( bidx(p+2,q+2,r+2) ),62);
                    tElement->insert_control_point( aPoints( bidx(p+1,q+2,r+2) ),63);
                }
            }
        }
    }

    void
    TensorMeshFactory::compute_tmatrix_quad4()
    {
        mTMatrix.set_size( 4,4, 0.0 );
        for ( uint k=0; k<4; ++k )
        {
            mTMatrix( k, k ) = 1.0 ;
        }
    }

    void
    TensorMeshFactory::compute_tmatrix_quad9()
    {
        // populate the node coordinates
        Matrix< index_t > idx =
            {
            { 0,2,2,0,1,2,1,0,1 },
                   { 0,0,2,2,0,1,2,1,1 }
            };

        Matrix< real > & T = mTMatrix ;
        T.set_size( 9,9,0.0 );
        Vector< real > Nxi( 3 );
        Vector< real > Neta( 3 );

        for ( index_t i=0; i<9; ++i )
        {
            // compute parameter coordinates of this point
            real xi  = idx(0,i ) - 1.0 ;
            real eta = idx(1,i) - 1.0 ;

            // evaluate 1D shape functions
            this->bspline_quad( xi, Nxi );
            this->bspline_quad( eta, Neta );

            // assemble matrix
            for ( index_t j=0; j<9; ++j )
            {
                T( i, j ) =
                    Nxi(idx(0,j))
                 * Neta(idx(1,j));
            }
        }
    }

    void
    TensorMeshFactory::compute_tmatrix_quad16()
    {
        // populate the node coordinates
        Matrix< index_t > idx =
        {
            { 0,3,3,0,1,2,3,3,2,1,0,0,1,2,2,1 },
                   { 0,0,3,3,0,0,1,2,3,3,2,1,1,1,2,2 }
        };

        Matrix< real > & T = mTMatrix ;
        T.set_size( 16,16,0.0 );
        Vector< real > Nxi( 4 );
        Vector< real > Neta( 4 );

        for ( index_t i=0; i<16; ++i )
        {
            // compute parameter coordinates of this point
            real xi  = ( idx(0,i )*2. - 3. ) / 3. ;
            real eta = ( idx(1,i )*2. - 3. ) / 3. ;

            // evaluate 1D shape functions
            this->bspline_quad( xi, Nxi );
            this->bspline_quad( eta, Neta );

            // assemble matrix
            for ( index_t j=0; j<16; ++j )
            {
                T( i, j ) =
                    Nxi(idx(0,j))
                 * Neta(idx(1,j));
            }
        }
    }

    void
    TensorMeshFactory::compute_tmatrix_hex8()
    {
        mTMatrix.set_size( 8,8, 0.0 );
        for ( uint k=0; k<8; ++k )
        {
            mTMatrix( k, k ) = 1.0 ;
        }
    }

    void
    TensorMeshFactory::compute_tmatrix_hex27()
    {
        // populate the node coordinates
        Matrix< index_t > idx =
        {
       {0,2,2,0,0,2,2,0,1,
               2,1,0,0,2,2,0,1,2,
               1,0,1,1,1,0,2,1,1},
              {0,0,2,2,0,0,2,2,0,
               1,2,1,0,0,2,2,0,1,
               2,1,1,1,1,1,1,0,2},
              {0,0,0,0,2,2,2,2,0,
               0,0,0,1,1,1,1,2,2,
               2,2,1,0,2,1,1,1,1}
        };

        Matrix< real > & T = mTMatrix ;
        T.set_size( 27,27,0.0 );
        Vector< real > Nxi( 3 );
        Vector< real > Neta( 3 );
        Vector< real > Nzeta( 3 );
        for ( index_t i=0; i<27; ++i )
        {
            // compute parameter coordinates of this point
            real xi   = idx(0,i) - 1.0 ;
            real eta  = idx(1,i) - 1.0 ;
            real zeta = idx(2,i) - 1.0 ;

            // evaluate 1D shape functions
            this->bspline_quad( xi, Nxi );
            this->bspline_quad( eta, Neta );
            this->bspline_quad( zeta, Nzeta );

            // assemble matrix
            for ( index_t j=0; j<27; ++j )
            {
                T( i, j ) =
                     Nxi(idx(0,j))
                 *  Neta(idx(1,j))
                 * Nzeta(idx(2,j));
            }
        }

    }

    void
    TensorMeshFactory::compute_tmatrix_hex64()
    {
        // populate the node coordinates
        Matrix< index_t > idx =
        {{0,3,3,0,0,3,3,0,1,2,0,0,0,0,3,3,
                 3,3,2,1,3,3,0,0,1,2,0,0,3,3,2,1,
                 1,1,2,2,1,2,2,1,0,0,0,0,3,3,3,3,
                 2,1,1,2,1,2,2,1,1,2,2,1,1,2,2,1},
                {0,0,3,3,0,0,3,3,0,0,1,2,0,0,1,2,
                 0,0,3,3,3,3,3,3,0,0,1,2,1,2,3,3,
                 1,2,2,1,0,0,0,0,1,1,2,2,1,2,2,1,
                 3,3,3,3,1,1,2,2,1,1,2,2,1,1,2,2},
                {0,0,0,0,3,3,3,3,0,0,0,0,1,2,0,0,
                 1,2,0,0,1,2,1,2,3,3,3,3,3,3,3,3,
                 0,0,0,0,1,1,2,2,1,2,2,1,1,1,2,2,
                 1,1,2,2,3,3,3,3,1,1,1,1,2,2,2,2}
        };
        Matrix< real > & T = mTMatrix ;
        T.set_size( 64,64,0.0 );
        Vector< real >   Nxi( 4 );
        Vector< real >  Neta( 4 );
        Vector< real > Nzeta( 4 );

        for ( index_t i=0; i<64; ++i )
        {
            // compute parameter coordinates of this point
            real  xi  = ( idx(0,i )*2. - 3. ) / 3. ;
            real  eta = ( idx(1,i )*2. - 3. ) / 3. ;
            real zeta = ( idx(2,i )*2. - 3. ) / 3. ;

            // evaluate 1D shape functions
            this->bspline_quad( xi, Nxi );
            this->bspline_quad( eta, Neta );
            this->bspline_quad( zeta, Nzeta );

            // assemble matrix
            for ( index_t j=0; j<64; ++j )
            {
                T( i, j ) =
                     Nxi(idx(0,j))
                 *  Neta(idx(1,j))
                 * Nzeta(idx(2,j));
            }
        }

    }

    void
    TensorMeshFactory::bspline_quad( const real xi, Vector< real > & aResult )
    {
        aResult( 0 ) = 0.125*(xi-1.)*(xi-1.);
        aResult( 1 ) = 0.25*(3.-xi*xi);
        aResult( 2 ) = 0.125*(xi+1)*(xi+1.);
    }

    void
    TensorMeshFactory::bspline_cub( const real xi, Vector< real > & aResult )
    {
        aResult( 0 ) = xi*(xi*(xi + 3.) + 3.) + 1.;
        aResult( 1 ) = xi * ( 15. - xi * ( 3.*xi + 3. )) +23. ;
        aResult( 2 ) = xi*(xi*(3*xi - 3) - 15. )+23. ;
        aResult( 3 ) = 1. + xi * ( xi * ( 3. -xi)  - 3. );
        aResult /= 48. ;
    }


    void
    TensorMeshFactory::connect_nodes_to_bsplines( Mesh * aMesh )
    {
        Cell< mesh::Node * > & tNodes = aMesh->nodes();
        Cell< mesh::ControlPoint * > & tControlPoints = aMesh->control_points();

        DynamicBitset tBitset( aMesh->number_of_control_points() );


        index_t tCount = 0 ;
        // make sure that indices are correct
        // and that everything is unflagged
        for ( mesh::ControlPoint * tControlPoint : tControlPoints )
        {
            tControlPoint->set_index( tCount++ );
            tControlPoint->unflag();
        }

        Cell< index_t > tIndices ;

        Vector< real > tWeigts ;
        Cell< mesh::Basis * > tSources ;

        Cell< std::pair< mesh::ControlPoint *, real > > tPairs ;

        for ( mesh::Node * tNode : tNodes )
        {
            tBitset.reset();

            // loop over all elements connected to this node
            for ( uint e=0; e<tNode->number_of_elements(); ++e )
            {
                // flag all control points connected to this element
                mesh::Element * tElement = tNode->element( e );

                for ( uint c=0; c<tElement->number_of_control_points(); ++c )
                {
                    mesh::ControlPoint * tControlPoint = tElement->control_point( c );
                    tControlPoint->flag() ;
                    tBitset.set( tControlPoint->index() );
                }
            }

            // now we collect the control points
            tBitset.where( tIndices );

            // next, we set the local indices
            tCount = 0 ;
            for ( index_t b=0; b<tIndices.size(); ++b )
            {
                tControlPoints( tIndices(b) )->set_level( tCount++ );
            }

            // now we assemble the T-Matrix
            for ( uint e=0; e<tNode->number_of_elements(); ++e )
            {
                // now we need to know which node we are
                mesh::Element * tElement = tNode->element( e );
                index_t i = gNoIndex ;
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    if ( tElement->node( k )->id() == tNode->id() )
                    {
                        i = k ;
                        break ;
                    }
                }

                for ( uint j=0; j<tElement->number_of_control_points(); ++j )
                {
                    mesh::ControlPoint * tControlPoint = tElement->control_point( j );
                    if ( tControlPoint->is_flagged() )
                    {
                        // get the weight
                        real tWeight = mTMatrix( i, j );

                        if ( std::abs( tWeight ) > BELFEM_MESH_EPSILON )
                        {
                            tPairs.push( std::pair< mesh::ControlPoint *, real >( tControlPoint, tWeight ) );
                        }

                        // we can unflag this point because the weight between a node and a control point
                        // is unique and not related to the element
                        tControlPoint->unflag();
                    }
                }
            }

            // sort sources after IDs. This is strictly not neccessary, but cleaner
            std::sort( tPairs.begin(), tPairs.end(), [](
                const std::pair< mesh::ControlPoint *, real > & a,
                const std::pair< mesh::ControlPoint *, real > & b ) { return a.first->id() < b.first->id(); } );

            tSources.set_size( tPairs.size() , nullptr );
            tWeigts.set_size( tPairs.size());

            // assemble vectors
            tCount = 0 ;
            for ( auto & tPair : tPairs )
            {
                tSources( tCount ) = tPair.first ;
                tWeigts( tCount ) = tPair.second ;
                ++tCount ;
            }

            // set dependencies
            tNode->set_sources( tSources, tWeigts );
        }

    }

    const Matrix< real > &
    TensorMeshFactory::t_matrix( const TensorMeshConfig * aConfig )
    {
        // check if T-matrix has been populated
        if ( aConfig == nullptr )
        {
            return mTMatrix ;
        }

        switch (  aConfig->element_type() )
        {
            case ElementType::QUAD4 :
            {
                this->compute_tmatrix_quad4();
                break ;
            }
            case ElementType::QUAD9 :
            {
                this->compute_tmatrix_quad9();
                break ;
            }
            case ElementType::QUAD16 :
            {
                this->compute_tmatrix_quad16();
                break ;
            }
            case ElementType::HEX8 :
            {
                this->compute_tmatrix_hex8();
                break ;
            }
            case ElementType::HEX27 :
            {
                this->compute_tmatrix_hex27();
                break ;
            }
            case ElementType::HEX64 :
            {
                this->compute_tmatrix_hex64();
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "invalid element type" );
            }
        }
        return mTMatrix ;
    }

}
