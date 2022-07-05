//
// Created by Christian Messe on 2019-08-18.
//
#include "fn_GT_data_path.hpp"

namespace belfem
{
    namespace gastables
    {
        string
        data_path()
        {
            // get the tpls environment variabl
            string tTPLS = std::getenv("TPLS");

            // go two steps up and add the cea directory
            return tTPLS.substr( 0, tTPLS.find_last_of('/') ) +"/cea";
        }
    }
}