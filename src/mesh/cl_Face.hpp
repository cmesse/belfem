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
            uint mIndexOnMaster ;

            Element * mSlave = nullptr ;
            uint mIndexOnSlave ;

            // orientation on master is always zero
            uint mOrientationOnSlave = gNoIndex ;

            Face * mPeriodic = nullptr ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            // 2d constructor
            Face( Element    * aParent );

            // 3d constructor
            Face( Element    * aMaster,
                  const uint aIndexOnMaster,
                  Element    * aSlave,
                  const uint aIndexOnSlave,
                  const uint aOrientationOnSlave = gNoIndex );

//-----------------------------------------------------------------------------

            ~Face() override;

//-----------------------------------------------------------------------------


            EntityType
            entity_type() const override ;

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
            edge_direction( const uint aEdgeIndex ) const ;

//-----------------------------------------------------------------------------

            uint
            number_of_corner_nodes() ;

//-----------------------------------------------------------------------------

            void
            flag_corner_nodes();

//-----------------------------------------------------------------------------

            size_t
            memory() const ;

//-----------------------------------------------------------------------------

            uint
            number_of_elements() const override ;

//-----------------------------------------------------------------------------

            Element *
            element(const uint aIndex) override;

//-----------------------------------------------------------------------------

            const Element *
            element(const uint aIndex) const override;

//------------------------------------------------------------------------------

            void
            set_master( Element * aElement, const uint aIndex );

 //------------------------------------------------------------------------------

            void
            set_slave( Element * aElement, const uint aIndex, const uint aOrientation );

            void
            set_periodic( Face * aPeriodic );

            Face *
            periodic() ;

            const Face *
            periodic() const ;

            bool
            is_periodic() const ;

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

        inline uint
        Face::number_of_elements() const
        {
            uint aNumElems = 0 ;
            if ( mMaster != nullptr ) ++ aNumElems ;
            if ( mSlave != nullptr ) ++ aNumElems ;
            return aNumElems ;
        }

//-----------------------------------------------------------------------------

        inline Element *
        Face::element(const uint aIndex)
        {
            if ( aIndex == 0 )
            {
                if ( mMaster != nullptr ) return mMaster ;
                if ( mSlave != nullptr ) return mSlave ;
                BELFEM_ASSERT( false, "Element index %u out of bounds.", aIndex );
                return nullptr ;
            }
            else
            {
                BELFEM_ASSERT( mMaster != nullptr, "Element index %u out of bounds.", aIndex );
                BELFEM_ASSERT( mSlave != nullptr, "Element index %u out of bounds.", aIndex );
                return mSlave ;
            }
        }

//-----------------------------------------------------------------------------

        inline const Element *
        Face::element(const uint aIndex) const
        {
            if ( aIndex == 0 )
            {
                if ( mMaster != nullptr ) return mMaster ;
                if ( mSlave != nullptr ) return mSlave ;
                BELFEM_ASSERT( false, "Element index %u out of bounds.", aIndex );
                return nullptr ;
            }
            else
            {
                BELFEM_ASSERT( mMaster != nullptr, "Element index %u out of bounds.", aIndex );
                BELFEM_ASSERT( mSlave != nullptr, "Element index %u out of bounds.", aIndex );
                return mSlave ;
            }
        }

//-----------------------------------------------------------------------------

        inline size_t
        Face::memory() const
        {
            return sizeof( Face ) + this->array_memory() ;
        }

//-----------------------------------------------------------------------------

        inline void
        Face::set_periodic( Face * aPeriodic )
        {
            mPeriodic = aPeriodic ;
        }

        inline Face *
        Face::periodic()
        {
            return mPeriodic ;
        }

        inline const Face *
        Face::periodic() const
        {
            return mPeriodic ;
        }

        inline bool
        Face::is_periodic() const
        {
            return mPeriodic != nullptr ;
        }

    }
}

#endif //BELFEM_CL_FACE_HPP
