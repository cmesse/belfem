//
// Created by christian on 3/1/23.
//
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_HDF5.hpp"
#include "cl_FEM_SideSet.hpp"
#include "FEM_geometry.hpp"
#include "fn_norm.hpp"
#include "fn_det.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "cl_Pipette.hpp"
#include "cl_Element_Factory.hpp"

using namespace belfem ;
using namespace fem ;

Communicator gComm;
Logger       gLog( 5 );

bool
EXPECT_NEAR( const real aA, const real aB, const real aEpsilon )
{
    const bool aCheck = std::abs( aA - aB ) < aEpsilon ;
    if( aCheck )
    {
        std::cout << "    pass! " << std::endl ;
    }
    else
    {
        std::cout << "    FAIL! | " << aA << " - " << aB << " | = " << std::abs( aA - aB )  <<  " > " << aEpsilon << std::endl ;
        BELFEM_ERROR( false, "fail");
    }
    return aCheck ;
}


//------------------------------------------------------------------------------

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

    print_banner();

    // todo: only test on one core
    // todo: make path relative
    const string tPath ="/home/christian/codes/belfem/test/fem/test_geometry.hdf5" ;

    // load the file
    HDF5 tFile( tPath, FileMode::OPEN_RDONLY );

    tFile.close() ;


    // close communicator
    return gComm.finalize();
}
