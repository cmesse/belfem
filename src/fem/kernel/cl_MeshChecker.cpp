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
#include "assert.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_MeshChecker.hpp"

#include "cl_FEM_Element.hpp"


namespace belfem
{
    MeshChecker::MeshChecker( Mesh* aMesh ) :
        mCommRank( comm_rank() ),
        mMesh( aMesh )
    {
        BELFEM_ERROR( mCommRank == 0, "The mesh checker must not be run in parallel" );
        BELFEM_ERROR( ! mMesh->edges_exist(), "The mesh checker must be run before edges have been created" );
        BELFEM_ERROR( ! mMesh->faces_exist(), "The mesh checker must be run before faces have been created" );

        // note that the mesh checker shouldn't run in parallel
        // because the mesh has not been distributed yet.
        mPipette = new mesh::Pipette();
        for ( mesh::Block * tBlock : mMesh->blocks() )
        {
            this->process_block( tBlock );
        }

        // a nonzero count means the input mesh ( usually gmsh ) contained
        // clockwise elements. Reported at Verbose only: on a clockwise mesh
        // every run reorients the same elements, and the reoriented mesh is
        // the correct one, so the count is a mesh diagnostic, not a warning.
        // 80-column layout: four blanks each side, content <= 72 chars
        if ( mElementCount > 0 )
        {
            message( InfoLevel::Verbose,
                "    MeshChecker : reoriented %lu elements with negative volume\n"
                "                  ( clockwise input mesh? )",
                ( long unsigned int ) mElementCount );
        }
        else
        {
            message( InfoLevel::Verbose,
                "    MeshChecker : all element volumes positive, nothing to reorient" );
        }

        // now we can tag this mesh as checked
        mMesh->set_mesh_checker_flag();
    }

    MeshChecker::~MeshChecker()
    {
        if ( mPipette != nullptr )
        {
            delete mPipette;
        }
    }

    void
    MeshChecker::link_to_block( mesh::Block * aBlock )
    {
        ElementType tType = aBlock->element_type() ;

        mPipette->set_element_type( tType );

        switch ( tType )
        {
            case ElementType::TRI3 :
            {
                mFunSwap = & MeshChecker::swap_tri3;
                break ;
            }
            case ElementType::TRI6 :
            {
                mFunSwap = & MeshChecker::swap_tri6;
                break ;
            }
            case ElementType::TRI10 :
            {
                mFunSwap = & MeshChecker::swap_tri10;
                break ;
            }
            case ElementType::TRI15 :
            {
                mFunSwap = & MeshChecker::swap_tri15;
                break ;
            }
            case ElementType::QUAD4 :
            case ElementType::QUAD4TS :
            {
                mFunSwap = & MeshChecker::swap_quad4;
                break ;
            }
            case ElementType::QUAD8 :
            case ElementType::QUAD9 :
            case ElementType::QUAD9TS :
            {
                mFunSwap = & MeshChecker::swap_quad9;
                break ;
            }
            case ElementType::QUAD16 :
            {
                mFunSwap = & MeshChecker::swap_quad16;
                break ;
            }
            case ElementType::TET4 :
            {
                mFunSwap = & MeshChecker::swap_tet4;
                break ;
            }
            case ElementType::TET10 :
            {
                mFunSwap = & MeshChecker::swap_tet10;
                break ;
            }
            case ElementType::TET20 :
            {
                mFunSwap = & MeshChecker::swap_tet20;
                break ;
            }
            case ElementType::TET35 :
            {
                mFunSwap = & MeshChecker::swap_tet35;
                break ;
            }
            case ElementType::PYRA5 :
            {
                mFunSwap = & MeshChecker::swap_pyra5 ;
                break ;
            }
            case ElementType::PYRA13 :
            case ElementType::PYRA14 :
            {
                mFunSwap = & MeshChecker::swap_pyra14 ;
                break ;
            }
            case ElementType::PENTA6 :
            case ElementType::PENTA6TS :
            {
                mFunSwap = & MeshChecker::swap_penta6 ;
                break ;
            }
            case ElementType::PENTA15 :
            {
                mFunSwap = & MeshChecker::swap_penta15 ;
                break ;
            }
            case ElementType::PENTA18 :
            case ElementType::PENTA18TS :
            {
                mFunSwap = & MeshChecker::swap_penta18 ;
                break ;
            }
            case ElementType::HEX8 :
            case ElementType::HEX8TS :
            {
                mFunSwap = & MeshChecker::swap_hex8 ;
                break ;
            }
            case ElementType::HEX8TB :
            {
                mFunSwap = & MeshChecker::swap_hex8tb ;
                break ;
            }
            case ElementType::HEX20 :
            {
                mFunSwap = & MeshChecker::swap_hex20 ;
                break ;
            }
            case ElementType::HEX27 :
            {
                mFunSwap = & MeshChecker::swap_hex27 ;
                break ;
            }
            case ElementType::HEX64 :
            {
                mFunSwap = & MeshChecker::swap_hex64 ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "No swapping function implemented for element type %s",
                    to_string( tType ).c_str() );
            }
        }
    }


    void
    MeshChecker::process_block( mesh::Block * aBlock )
    {
        // flag them if their volumes are negative
        this->link_to_block( aBlock );

        //std::cout << "Checking block " << aBlock->id() << " " << to_string( aBlock->element_type() ) << std::endl;
        Cell< mesh::Element * > & tElements = aBlock->elements();
        for ( mesh::Element * tElement : tElements )
        {
            if ( mPipette->measure( tElement ) < 0 )
            {
                this->swap( tElement );
                ++mElementCount ;
            }
        }
    }

}
