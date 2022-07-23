//
// Created by Christian Messe on 14.11.19.
//

#include "cl_MaterialFactory.hpp"
#include "cl_IsotropicMaterial.hpp"
#include "cl_Material_Aluminum.hpp"
#include "cl_Material_CCSIC.hpp"
#include "cl_Material_AltraMat80.hpp"
#include "cl_Material_Rohacell51.hpp"
#include "cl_Material_Inconel718.hpp"
#include "cl_Material_Copper.hpp"
#include "cl_Material_Ti6Al4V.hpp"
#include "cl_Material_CuCrZr.hpp"
#include "cl_Material_Zirconia.hpp"
#include "cl_Material_Inconel750X.hpp"
#include "cl_Material_Hastelloy.hpp"
#include "cl_Material_Silver.hpp"
#ifdef BELFEM_GASMODELS
#include "cl_Material_Air.hpp"
#endif


namespace belfem
{
//----------------------------------------------------------------------------

    Material *
    MaterialFactory::create_material( const MaterialType aMaterial )
    {
        switch( aMaterial )
        {
            case ( MaterialType::Simple ) :
            {
                return create_isotropic_material( aMaterial );
            }
#ifdef BELFEM_GASMODELS
            case ( MaterialType::Air ) :
            {
                return new material::Air ;
            }
#endif
            case ( MaterialType::Aluminum ) :
            {
                return new material::Aluminum ;
            }
            case ( MaterialType::Inconel718 ) :
            {
                return new material::Inconel718 ;
            }
            case ( MaterialType::Copper ) :
            {
                return new material::Copper ;
            }
            case ( MaterialType::TI6AL4V ) :
            {
                return new material::Ti6Al4V ;
            }
            case ( MaterialType::CCSIC ) :
            {
                return new material::CCSIC ;
            }
            case( MaterialType::Altramat80 ) :
            {
                return new material::AltraMat80 ;
            }
            case( MaterialType::Rohacell51 ) :
            {
                return new material::Rohacell51 ;
            }
            case( MaterialType::CuCrZr ) :
            {
                return new material::CuCrZr ;
            }
            case( MaterialType::Zirconia ) :
            {
                return new material::Zirconia ;
            }
            case( MaterialType::Inconel750X ) :
            {
                return new material::Inconel750X ;
            }
            case( MaterialType::Hastelloy ) :
            {
                return new material::Hastelloy ;
            }
            case( MaterialType::Silver ) :
            {
                return new material::Silver ;
            }
            default :
            {
                BELFEM_ERROR( false, "Unknown Material or material not implemented in factory." );
                return nullptr;
            }
        }
    }

//----------------------------------------------------------------------------

    Material *
    MaterialFactory::create_material( const string & aLabel )
    {
        return this->create_material( string_to_material_type( aLabel ) );
    }

//----------------------------------------------------------------------------

    Material *
    MaterialFactory::create_isotropic_material( const MaterialType aMaterial )
    {
        IsotropicMaterial * aMat = nullptr;

        switch( aMaterial )
        {
            case ( MaterialType::Simple ) :
            {
                aMat = new IsotropicMaterial( MaterialType::Simple );

                aMat->mLabel = "Simple";
                aMat->mYoungPoly.set_size( 1, 70.0e9 );
                aMat->mShearPoly.set_size( 1, 27.5e9 );
                aMat->mSpecificHeatPoly.set_size( 1, 900.0 );
                aMat->mThermalConductivityPoly.set_size( 1, 235.0 );
                aMat->mDensityPoly.set_size( 1, 2700.0 );

                break;
            }
            default :
            {
                BELFEM_ERROR( false, "Unknown Material" );
                break;
            }
        }

        return aMat;
    }

//----------------------------------------------------------------------------

}
