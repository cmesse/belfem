//
// Created by christian on 9/20/21.
//

#include "cl_FEM_DomainGroup.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        DomainGroup::DomainGroup(
                const DomainType aType,
                const string & aLabel,
                const Vector< id_t > & aGroupIDs  ) :
                mType( aType ),
                mLabel( aLabel ),
                mGroupIDs( aGroupIDs )
        {

        }

//------------------------------------------------------------------------------

        const string &
        DomainGroup::label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        const Vector< id_t > &
        DomainGroup::groups() const
        {
            return mGroupIDs ;
        }

//------------------------------------------------------------------------------

        void
        DomainGroup::set_groups( const Vector< id_t > & aGroups )
        {
            mGroupIDs = aGroups ;
        }

//------------------------------------------------------------------------------

        id_t
        DomainGroup::master_block_id() const
        {
            BELFEM_ERROR( false,
                         "Invalid call to base class function DomainGroup::master_block_id()") ;
            return gNoID ;
        }

//------------------------------------------------------------------------------

        id_t
        DomainGroup::slave_block_id() const
        {
            BELFEM_ERROR( false,
                         "Invalid call to base class function DomainGroup::slave_block_id()") ;
            return gNoID ;
        }

//------------------------------------------------------------------------------

        const Vector< id_t > &
        DomainGroup::plus() const
        {
            BELFEM_ERROR( false,
                         "Invalid call to base class function DomainGroup::plus()") ;
            return mGroupIDs ;
        }


//------------------------------------------------------------------------------

        const Vector< id_t > &
        DomainGroup::minus() const
        {
            BELFEM_ERROR( false,
                         "Invalid call to base class function DomainGroup::minus()") ;
            return mGroupIDs ;
        }

//------------------------------------------------------------------------------

        const Vector< id_t > &
        DomainGroup::master() const
        {
            BELFEM_ERROR( false,
                          "Invalid call to base class function DomainGroup::master()") ;
            return mGroupIDs ;
        }

//------------------------------------------------------------------------------

        DomainCut::DomainCut( const string & aLabel,
                   const Vector< id_t > & aGroupIDs,
                   const Vector< id_t > & aPlus,
                   const Vector< id_t > & aMinus,
                  const DomainType aType ) :
                DomainGroup( aType, aLabel, aGroupIDs ),
                mPlus( aPlus ),
                mMinus( aMinus )
        {

        }

//------------------------------------------------------------------------------

        DomainThinShell::DomainThinShell( const string & aLabel,
                         const Vector< id_t > & aGroupIDs,
                         const Vector< id_t > & aPlus,
                         const Vector< id_t > & aMinus,
                         const Vector< id_t > & aMaster,
                         const DomainType aType ) :
                DomainCut( aLabel, aGroupIDs, aPlus, aMinus, aType ) ,
                mMaster( aMaster )
        {

        }

//------------------------------------------------------------------------------

        DomainSolid::DomainSolid(
                const DomainType aType,
                const string & aLabel,
                const Vector< id_t > & aGroupIDs ) :
                DomainGroup( aType, aLabel, aGroupIDs )
        {

        }

//------------------------------------------------------------------------------

        DomainInterface::DomainInterface(
                const  DomainType aType,
                const id_t aSideSetID,
                const id_t aMasterBlockID,
                const id_t aSlaveBlockID ) :
                DomainGroup( aType,
                             "Interface",
                             Vector< id_t >( 1, aSideSetID ) ) ,
                             mMasterBlockID( aMasterBlockID ),
                             mSlaveBlockID( aSlaveBlockID )
        {
        }

//------------------------------------------------------------------------------
    }
}