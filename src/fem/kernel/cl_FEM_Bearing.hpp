//
// Created by Christian Messe on 15.11.19.
//

#ifndef BELFEM_CL_FEM_BEARING_HPP
#define BELFEM_CL_FEM_BEARING_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_FEM_Dof.hpp"

namespace belfem
{
    namespace fem
    {
        class DofManagerBase;
        class Field ;
        class DofManager ;

//------------------------------------------------------------------------------

        class Bearing
        {
            DofManagerBase * mParent    = nullptr;

            const id_t mID;

            mesh::Node * mNode = nullptr;

            Cell< Dof * > mDOFs;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // constructor for empty bearing
            Bearing( DofManagerBase * aParent );

//------------------------------------------------------------------------------

            Bearing( Field * aParent, const id_t aID, mesh::Node * aNode );

//------------------------------------------------------------------------------

            Bearing( DofManager * aParent, const id_t aID, mesh::Node * aNode );

//------------------------------------------------------------------------------

            ~Bearing() = default;

//------------------------------------------------------------------------------

            id_t
            id() const;

//------------------------------------------------------------------------------

            mesh::Node *
            node();

//------------------------------------------------------------------------------

            void
            impose_dirichlet( const real aValue, const uint aDofType=0 );

//------------------------------------------------------------------------------

            void
            free();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline id_t
        Bearing::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline mesh::Node *
        Bearing::node()
        {
            return mNode;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_BEARING_HPP
