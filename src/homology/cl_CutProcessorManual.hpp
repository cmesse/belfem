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

#ifndef CL_CUTPROCESSORMANUAL_HPP
#define CL_CUTPROCESSORMANUAL_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {

        class CutProcessorManual
        {
            Mesh * mMesh ;
            Vector< id_t > mAirBlocks ;
            Vector< id_t > mCutSideSets ;
            Vector< id_t > mCutMasterVolumes ;
            Cell< Node * > mAbstractNodes ;
            Cell< Node * > mOriginalNodes ;
            Cell< Node * > mDuplicateNodes ;
            Cell< Node * > mOrphanNodes ;
            Matrix< real > mIncidentMatrix ;

            struct ManualSet
            {
                 Node * mAbstractNode = nullptr ;
                 Cell< Node * > mOriginalNodes ;
                 Cell< Node * > mDuplicateNodes ;
                 Cell< Element * > mElements ;

            };
            Cell< ManualSet * > mSets ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------
            CutProcessorManual( Mesh * aMesh,
                               const Vector< id_t > & aAirBlocks,
                               const Vector< id_t > & aCutSideSets,
                               const Vector< id_t > & aCutMasterVolumes );

            ~CutProcessorManual();

            Cell< Node * > &
            abstract_nodes();

            Cell< Node * > &
            interface_originals();

            Cell< Node * > &
            interface_duplicates();

            Cell< Node * > &
            orphaned_nodes();

            const Vector< id_t > &
            air_blocks() const ;

            const Matrix< real > &
            incident_matrix() const ;

        private:

            ManualSet *
            create_set( const index_t aIndex, id_t & aMaxNodeID );


        };

        inline Cell< Node * > &
        CutProcessorManual::abstract_nodes()
        {
            return mAbstractNodes;
        }

        inline Cell< Node * > &
        CutProcessorManual::orphaned_nodes()
        {
            return mOrphanNodes;
        }

        inline const Vector<id_t>& CutProcessorManual::air_blocks() const
        {
            return mAirBlocks;
        }

        inline const Matrix< real > &
        CutProcessorManual::incident_matrix() const
        {
            return mIncidentMatrix;
        }

    }
}
#endif //CL_CUTPROCESSORMANUAL_HPP
