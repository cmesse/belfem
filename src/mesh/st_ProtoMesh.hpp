/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_ST_PROTOMESH_HPP
#define BELFEM_ST_PROTOMESH_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Map.hpp"
#include "Mesh_Enums.hpp"
#include "en_DomainType.hpp"

namespace belfem
{
    namespace mesh
    {
        namespace proto
        {
            struct MetaData
            {
                uint         mNumberOfDimensions ;
                Cell< uint > mNumberOfEntities ;
                Cell< uint > mNumberOfGroups ;
                size_t       mChecksum = 0 ;

                //! provenance of the writing build ( empty on older files ):
                //! the checksum is blind to the construction code, so "which
                //! build wrote this file" is the first question on any
                //! reload oddity
                string       mBelfemVersion ;
                string       mGitHash ;

                //! fingerprint of the SETTINGS the mesh was enriched with, and
                //! the canonical text behind it. The checksum above identifies
                //! the base mesh only, so these are what catch an edited deck
                //! on an unchanged .msh. Zero / empty on files written before
                //! the tag existed
                uint64_t     mConfigTag = 0 ;
                string       mConfigText ;
            };

            struct NodeData
            {
                Cell< id_t >      mIDs ;
                Cell< proc_t >    mOwners ;
                Matrix< real >    mCoords ;
                Cell< id_t >      mDuplicateData ;
            };

            struct ControlPointData
            {
                Cell< id_t >      mIDs ;
                Cell< proc_t >    mOwners ;
                Matrix< real   >  mCoords ;
                Cell< id_t >      mElementTopology ;
            };

            struct TMatrixData
            {
                Cell< index_t >   mNumTargets ;
                Cell< id_t >      mTargetIDs ;
                Cell< uint >      mCounters ;
                Cell< id_t >      mSourceIDs ;
                Cell< uchar >     mTypes ;
                Cell< real >      mWeights ;
            };

            struct ElementData
            {
                Cell< id_t >    mIDs ;
                Cell< proc_t >  mOwners ;
                Cell< uchar >   mTypes ;
                Cell< uint >    mGeometryTags ;
                Cell< uint >    mPhysicalTags ;
                Cell< id_t >    mTopology ;
            };

            struct ElementExtra
            {
                Cell< id_t >   mNeighborData ;
                Cell< id_t >   mEdgeData ;
                Cell< id_t >   mFaceData ;
                Cell< uchar >  mEdgeDirections ;
                Cell< id_t >   mCurvedElementIDs ;
            };

            struct EdgeData
            {
                Cell< id_t >     mIDs ;
                Cell< proc_t >   mOwners ;
                Cell< id_t >     mTopology ;
            };

            struct FaceData
            {
                Cell< id_t >    mIDs ;
                Cell< proc_t >  mOwners ;
                Cell< id_t >    mMasterIDs ;
                Cell< id_t >    mSlaveIDs ;
                Cell< uchar >   mIndicesOnMaster ;
                Cell< uchar >   mIndicesOnSlave ;
                Cell< uchar >   mOrientationsOnSlave ;
            };

            struct FacetData : FaceData
            {
                Cell< uchar >  mTypes ;
                Cell< uint >   mGeometryTags ;
                Cell< uint >   mPhysicalTags ;
                Cell< id_t >   mTopology ;
            };

            struct FacetExtra
            {
                Cell< id_t >   mNeighborData ;
                Cell< id_t >   mCurvedFacetIDs ;
            };

            struct ThinShellData
            {
                id_t           mSideSetID ;
                id_t           mGhostSideSetID = gNoID ;
                Cell< id_t >   mBlocksIDs ;
                Cell< real >   mThicknesses ;
                Cell< string > mMaterials ;
            };

            /*struct SourceData
            {
                Cell< id_t >   mNodeData ;
                Cell< id_t >   mEdgeData ;
                Cell< id_t >   mFaceData ;
                Cell< id_t >   mElementData ;
                Cell< id_t >   mFacetData ;
            };*/

            struct GroupData
            {
                id_t         mID = gNoID ;
                index_t      mNumElements = 0 ;
                string       mLabel ;
                ElementType  mElementType = ElementType::EMPTY ;
                DomainType   mDomainType  = DomainType::Default ;
                bool         mHidden   = false ;
                bool         mHasEdges = false ;
                bool         mHasFaces = false ;

                //! physical block thickness ( thin-shell layers and side
                //! connector walls; NaN for volume blocks ). Must travel to
                //! the other procs: EF_HEX8TB reads the wall width from it
                real         mThickness = BELFEM_QUIET_NAN ;
            };

            struct PeriodicityData
            {
                Cell< id_t >   mMasterPlane ;
                Cell< id_t >   mSlavePlane ;
                Cell< id_t >   mMasterSideSets ;
                Cell< id_t >   mSlaveSideSets ;
                Cell< id_t >   mMasterNodes ;
                Cell< id_t >   mSlaveNodes ;
                Cell< id_t >   mMasterEdges ;
                Cell< id_t >   mSlaveEdges ;
                Cell< id_t >   mMasterFaces ;
                Cell< id_t >   mSlaveFaces ;
                Cell< id_t >   mMasterFacets ;
                Cell< id_t >   mSlaveFacets ;
            };

            // todo: Still missing: CurveData
            // todo: Still missing: function to provide tables to Kernel

        }
    }
}
#endif //BELFEM_ST_PROTOMESH_HPP
