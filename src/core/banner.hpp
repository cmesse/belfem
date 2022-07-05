//
// Created by Christian Messe on 2018-12-30.
//

#ifndef BELFEM_BANNER_HPP
#define BELFEM_BANNER_HPP

#include <string>
namespace belfem
{
//------------------------------------------------------------------------------

    std::string
    exec( const std::string & aCommand );

//------------------------------------------------------------------------------

    /**
     * returns the Unix version ( Linux or Darwin )
     */
    std::string
    uname();

//------------------------------------------------------------------------------

    /**
     * grabs the cpu info for the banner
     */
    std::string
    cpu_info();

//------------------------------------------------------------------------------

    /**
     * prints the banner
     */
    void
    print_banner( const std::string aExecName = "" );

//------------------------------------------------------------------------------
}
#endif //BELFEM_BANNER_HPP