//
// Created by Christian Messe on 25.08.19.
//

#ifndef BELFEM_CL_GT_INPUTTRANSPORT_HPP
#define BELFEM_CL_GT_INPUTTRANSPORT_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_TransportPoly.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        class InputTransport: public Ascii
        {
            Map <string, uint> mMap;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * Default Constructor. Loads an input file from a path
             * @param aPath
             */
            InputTransport( const string & aPath );

            ~InputTransport() = default;

//------------------------------------------------------------------------------

            void
            read_data( RefGas * aRefgas );

//------------------------------------------------------------------------------

            uint
            entry_exists( RefGas * aRefgas  );

//------------------------------------------------------------------------------
            /*
             * test if interaction parameter exists
             */
            bool
            interaction_parameter_exists(
                    const string & aLabelA,
                    const string & aLabelB );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            TransportPoly *
            read_polynomial( uint & aLineCount );

//------------------------------------------------------------------------------

            real
            word_to_real( const string & aWord );

        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_INPUTTRANSPORT_HPP
