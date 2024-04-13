//
// Created by christian on 4/10/24.
//

#ifndef BELFEM_CL_MESHOPTIONS_HPP
#define BELFEM_CL_MESHOPTIONS_HPP

#include "typedefs.hpp"

namespace belfem
{
    class MeshOptions
    {
        bool mComputeConnectivities = true ;
        bool mCreateEdges = false ;
        bool mCreateFaces = false ;
        bool mCheckElements = true ;

        string mUnitString = "mm" ;
        real   mUnitScale   = 0.001 ;

    public:

         MeshOptions() = default ;

        ~MeshOptions() = default ;
    };
}

#endif //BELFEM_CL_MESHOPTIONS_HPP
