//
// Created by Christian Messe on 02.12.19.
//

#include "cl_Atmosphere.hpp"
#include "cl_AtmosphereModel_ISA1976.hpp"
#include "assert.hpp"
namespace belfem
{
//------------------------------------------------------------------------------
    Atmosphere::Atmosphere( const belfem::AtmosphereType aType ) :
            mType( aType )
    {
        switch ( aType )
        {
            case ( AtmosphereType::ISA1976 ) :
            {
                mModel = new atmosphere::AtmosphereModel_ISA1976();
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown atmosphere model " );
            }
        }
    }

//------------------------------------------------------------------------------

    Atmosphere::~Atmosphere()
    {
        delete mModel;
    }

//------------------------------------------------------------------------------

    void
    Atmosphere::compute_T_and_P( const real & aAltitude,
                                 real & aTemperature,
                                 real & aPressure ) const
    {
        mModel->compute_T_and_p( aAltitude, aTemperature, aPressure );
    }

//------------------------------------------------------------------------------
}