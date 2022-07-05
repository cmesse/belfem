//
// Created by christian on 10/28/21.
//

#ifndef BELFEM_FN_FEM_BIOT_SAVART_HPP
#define BELFEM_FN_FEM_BIOT_SAVART_HPP

#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
        namespace biotsavart
        {
            void
            identify_current_blocks(  DofManager * aField,
                                      Vector< id_t > & aAirBlocks,
                                      Vector< id_t > & aCoilBlocks,
                                      Vector< id_t > & aScBlocks );

            void
            get_node_coords(
                    DofManager * aField,
                    Vector< index_t >    & aNodeIndices,
                    Matrix< real >       & aNodeCoords );

            void
            get_sideset_node_coords(
                    DofManager * aField,
                    const Vector< id_t > & aSideSetIDs,
                    Vector< index_t >    & aNodeIndices,
                    Matrix< real >       & aNodeCoords );

            void
            compute_biot_savart_2d(
                    DofManager * aField,
                    const Cell< string >       & aBFields,
                    const Vector< id_t >       & aCoilBlocks,
                    const Vector< id_t >       & aScBlocks,
                    const Vector< index_t >    & aNodeIndices,
                    const Matrix< real >       & aNodeCoords );

            void
            compute_biot_savart_3d(
                    DofManager * aField,
                    const Cell< string >       & aBFields,
                    const Vector< id_t >       & aCoilBlocks,
                    const Vector< id_t >       & aScBlocks,
                    const Vector< index_t >    & aNodeIndices,
                    const Matrix< real >       & aNodeCoords );
        }

        void
        compute_biot_savart(  DofManager * aField );

    }
}
#endif //BELFEM_FN_BIOT_SAVART_HPP
