//
// Created by christian on 4/15/22.
//

#ifndef BELFEM_ST_LAYER_HPP
#define BELFEM_ST_LAYER_HPP

namespace belfem
{
    namespace mesh
    {
        /**
         *  this is a supporting data type for the tape generation
         *  the mesh entities are passed over to the mesh
         *  and will be deleted by the mesh during destruction
         */
        struct Layer
        {
            // cell with cloned nodes
            Cell< Node * >    Nodes ;
            Cell< Edge * >    Edges ;
            Cell< Face * >    Faces ;
            Cell< Element * > Elements ;
        };
    }

}



#endif //BELFEM_ST_LAYER_HPP
