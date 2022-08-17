//
// Created by Christian Messe on 17.08.20.
//

#include <iostream>


#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"

#include "fn_linspace.hpp"
#include "fn_r2.hpp"

#define private public
#define protected public
#include "cl_GM_EoS_Hydrogen.hpp"
#include "cl_GM_EoS_Oxygen.hpp"
#include "cl_GM_EoS_Methane.hpp"
#undef protected
#undef private

#include "cl_Gas.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_GM_Helmholz_DerivTest.hpp"
#include "fn_r2.hpp"
using namespace belfem;
using namespace gasmodels;

Communicator gComm;
Logger       gLog( 3 );

#define EXPECT_NEAR( A, B, C ) { std::cout << A << std::endl; };
int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( argc, argv );

//----------------------------------------------------------------------------


    Gas tGas( "O2", GasModel::HELMHOLTZ );

    std::cout << tGas.rho( 90.0, 1e5 ) << std::endl ;
    std::cout << tGas.eos()->hvap( 90.0, 1e5 ) << std::endl ;

//----------------------------------------------------------------------------

    return gComm.finalize();
}
