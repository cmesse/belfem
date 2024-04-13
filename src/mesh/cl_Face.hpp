//
// Created by christian on 7/12/21.
//

#ifndef BELFEM_CL_FACE_HPP
#define BELFEM_CL_FACE_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_Vertex.hpp"

namespace belfem
{
    namespace mesh
    {
        class Element ;

//-----------------------------------------------------------------------------

        class Face : public Vertex
        {

            Element * mMaster   = nullptr ;
            const uint mIndexOnMaster ;

            Element * mSlave = nullptr ;
            const uint mIndexOnSlave ;

            // orientation on master is always zero
            const uint mOrientationOnSlave ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            // 2d constructor
            Face( Element    * aParent );

            // 3d constructor
            Face( Element    * aMaster,
                  const uint aIndexOnMaster,
                  Element    * aSlave,
                  const uint aIndexOnSlave );

//-----------------------------------------------------------------------------

            ~Face();

//-----------------------------------------------------------------------------


            EntityType
            entity_type() const ;

//-----------------------------------------------------------------------------

            Element *
            master();

//-----------------------------------------------------------------------------

            Element *
            slave();

//-----------------------------------------------------------------------------

            uint
            index_on_master() const ;

//-----------------------------------------------------------------------------

            uint
            index_on_slave() const ;

//-----------------------------------------------------------------------------

            uint
            orientation_on_slave() const ;

//-----------------------------------------------------------------------------

            bool
            edge_orientation( const uint aEdgeIndex ) const ;

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            uint
            compute_orientation(
                    Element * aMaster,
                    const uint aIndexOnMaster,
                    Element * aSlave,
                    const uint aIndexOnSlave );

//-----------------------------------------------------------------------------
        };
//-----------------------------------------------------------------------------

        inline EntityType
        Face::entity_type() const
        {
            return EntityType::FACE ;
        }

//-----------------------------------------------------------------------------

        inline Element *
        Face::master()
        {
            return mMaster ;
        }

//-----------------------------------------------------------------------------

        inline Element *
        Face::slave()
        {
            return mSlave ;
        }

//-----------------------------------------------------------------------------

        inline uint
        Face::index_on_master() const
        {
            return mIndexOnMaster ;
        }

//-----------------------------------------------------------------------------

        inline uint
        Face::index_on_slave() const
        {
            return mIndexOnSlave ;
        }

//-----------------------------------------------------------------------------

        inline uint
        Face::orientation_on_slave() const
        {
            return mOrientationOnSlave ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FACE_HPP
