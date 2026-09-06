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

#ifndef BELFEM_CL_TENSORMESHFACTORY_HPP
#define BELFEM_CL_TENSORMESHFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    class TensorMeshFactory
    {
        const TensorMeshConfig * mConfig = nullptr ;

        Matrix< real > mTMatrix ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        TensorMeshFactory() = default ;

//------------------------------------------------------------------------------

        ~TensorMeshFactory() = default ;

//------------------------------------------------------------------------------

        void
        populate_tensor_mesh(
            const TensorMeshConfig * aConfig,
            Mesh * aMesh,
            const proc_t aMasterProc = 0 );

        void
        create_topology(
            const TensorMeshConfig * aConfig,
            Matrix< index_t > & aTopology ) ;

        void
        create_grid(
            const Vector< index_t > & aNumNodes,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
            Matrix< real > & aGrid ) const ;

        void
        create_bsplines( Mesh * aMesh );

        const Matrix< real > &
        t_matrix( const TensorMeshConfig * aConfig = nullptr ) ;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        create_topology_quad4( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const ;

        void
        create_topology_quad9( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const ;

        void
        create_topology_quad16( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const ;

        void
        create_topology_hex8( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const ;

        void
        create_topology_hex27( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const ;

        void
        create_topology_hex64( const Vector< index_t > & aNumElems, Matrix< index_t > & aTopology ) const;

        void
        create_bsplines_quad9( const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints );

        void
        create_bsplines_quad16( const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints );

        void
        create_bsplines_hex27( const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints );

        void
        create_bsplines_hex64( const TensorMeshConfig * aConfig, Cell< mesh::ControlPoint *  > & aPoints );

        void
        link_bsplines_quad9(
            const TensorMeshConfig * aConfig,
            Cell< mesh::ControlPoint *  > & aPoints,
            Cell< mesh::Element * > & aElements );

        void
        link_bsplines_quad16(
            const TensorMeshConfig * aConfig,
            Cell< mesh::ControlPoint *  > & aPoints,
            Cell< mesh::Element * > & aElements );

        void
        link_bsplines_hex27(
            const TensorMeshConfig * aConfig,
            Cell< mesh::ControlPoint *  > & aPoints,
            Cell< mesh::Element * > & aElements );

        void
        link_bsplines_hex64(
           const TensorMeshConfig * aConfig,
           Cell< mesh::ControlPoint *  > & aPoints,
           Cell< mesh::Element * > & aElements );


        index_t
        nidx( const index_t i, const index_t j ) const;

        index_t
        nidx( const index_t i, const index_t j, const index_t k ) const;

        index_t
        bidx( const index_t p, const index_t q ) const;

        index_t
        bidx( const index_t p, const index_t q, const index_t r ) const;

        void
        compute_tmatrix_quad4();

        void
        compute_tmatrix_quad9();

        void
        compute_tmatrix_quad16();

        void
        compute_tmatrix_hex8();

        void
        compute_tmatrix_hex27();

        void
        compute_tmatrix_hex64();

        void
        bspline_quad( const real xi, Vector< real > & aResult );

        void
        bspline_cub( const real xi, Vector< real > & aResult );

        void
        connect_nodes_to_bsplines( Mesh * aMesh );

    };

    inline index_t
    TensorMeshFactory::nidx( const index_t i, const index_t j ) const
    {
        return i * mConfig->num_nodes( 1 ) + j ;
    }

    inline index_t
    TensorMeshFactory::nidx( const index_t i, const index_t j, const index_t k ) const
    {
        return mConfig->num_nodes( 2 ) *
            ( i * mConfig->num_nodes( 1 ) +  j ) + k ;
    }

    inline index_t
    TensorMeshFactory::bidx( const index_t p, const index_t q ) const
    {
        return p * mConfig->num_control_points( 1 ) + q ;
    }

    inline index_t
    TensorMeshFactory::bidx( const index_t p, const index_t q, const index_t r ) const
    {
        return mConfig->num_control_points( 2 ) *
            ( p * mConfig->num_control_points( 1 ) +  q ) + r ;
    }

}

#endif //BELFEM_CL_TENSORMESHFACTORY_HPP
