//
// Created by christian on 9/20/21.
//

#ifndef BELFEM_CL_FEM_DOMAINGROUP_HPP
#define BELFEM_CL_FEM_DOMAINGROUP_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "en_FEM_DomainType.hpp"

namespace belfem
{
    namespace fem
    {
        class DomainGroup
        {
            //! the type of this domain
            const DomainType mType ;

            //! the label of this domain
            const string mLabel ;

            //! IDs of associated blocks or sidesets on the mesh
            Vector< id_t > mGroupIDs ;

            //! multi-purpose flag
            //! needed by Maxwell Factory to tell if sideset
            //! has been associated with boundary conditions
            bool mFlag = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DomainGroup( const DomainType aType,
                         const string & aLabel,
                         const Vector< id_t > & aGroupIDs );

            virtual
            ~DomainGroup() = default ;

//------------------------------------------------------------------------------

            const string &
            label() const ;

            const Vector< id_t > &
            groups() const ;

//------------------------------------------------------------------------------

            /**
             * manually overwrite groups
             */
             void
             set_groups( const Vector< id_t > & aGroups );

//------------------------------------------------------------------------------

            /**
             * return the type
             */
            DomainType
            type() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            virtual const Vector< id_t > &
            plus() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            virtual const Vector< id_t > &
            minus() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            virtual const Vector< id_t > &
            master() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            virtual id_t
            master_block_id() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            virtual id_t
            slave_block_id() const ;

//------------------------------------------------------------------------------

            void
            flag();

//------------------------------------------------------------------------------

            void
            unflag();

//------------------------------------------------------------------------------

            bool
            is_flagged() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        class DomainCut : public DomainGroup
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            const Vector< id_t > mPlus ;
            const Vector< id_t > mMinus ;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DomainCut( const string & aLabel,
                       const Vector< id_t > & aGroupIDs,
                       const Vector< id_t > & aPlus,
                       const Vector< id_t > & aMinus,
                       const DomainType aType = DomainType::Cut );

//------------------------------------------------------------------------

            ~DomainCut() = default ;

//------------------------------------------------------------------------------

            // the block IDs for the plus side
            const Vector< id_t > &
            plus() const ;

//------------------------------------------------------------------------------

            // the block IDs for the minus side
            const Vector< id_t > &
            minus() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        class DomainThinShell : public DomainCut
        {
            const Vector< id_t > mMaster ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DomainThinShell( const string & aLabel,
                       const Vector< id_t > & aGroupIDs,
                       const Vector< id_t > & aPlus,
                       const Vector< id_t > & aMinus,
                       const Vector< id_t > & aMaster,
                       const DomainType aType = DomainType::ThinShell );

//------------------------------------------------------------------------

            ~DomainThinShell() = default ;

//------------------------------------------------------------------------

            // the master blocks for the orientation
            const Vector< id_t > &
            master() const ;

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        class DomainSolid : public DomainGroup
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DomainSolid(
                    const DomainType aType,
                    const string & aLabel,
                    const Vector< id_t > & aGroupIDs );

            ~DomainSolid() = default ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        class DomainInterface : public DomainGroup
        {
            const id_t mMasterBlockID ;
            const id_t mSlaveBlockID ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DomainInterface(
                    const  DomainType aType,
                    const id_t aSideSetID,
                    const id_t aMasterBlockID,
                    const id_t aSlaveBlockID );

            ~DomainInterface() = default ;


//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            id_t
            master_block_id() const ;

//------------------------------------------------------------------------------

            // these functions do nothing for the base class
            id_t
            slave_block_id() const ;
//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline DomainType
        DomainGroup::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        inline void
        DomainGroup::flag()
        {
            mFlag = true ;
        }

//------------------------------------------------------------------------------

        inline void
        DomainGroup::unflag()
        {
            mFlag = false ;
        }

//------------------------------------------------------------------------------

        inline bool
        DomainGroup::is_flagged() const
        {
            return mFlag ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        DomainCut::plus() const
        {
            return mPlus ;
        }

//------------------------------------------------------------------------------

        // the block IDs for the minus side
        inline const Vector< id_t > &
        DomainCut::minus() const
        {
            return mMinus ;
        }

//------------------------------------------------------------------------------

        inline id_t
        DomainInterface::master_block_id() const
        {
            return mMasterBlockID ;
        }

//------------------------------------------------------------------------------

        inline id_t
        DomainInterface::slave_block_id() const
        {
            return mSlaveBlockID ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        DomainThinShell::master() const
        {
            return mMaster ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_DOMAINGROUP_HPP
