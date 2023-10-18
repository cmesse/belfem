//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_EN_IWGS_HPP
#define BELFEM_EN_IWGS_HPP

namespace belfem
{
//------------------------------------------------------------------------------

    enum class IwgType
    {
        SimpleDiffusion,
        StaticHeatConduction,       // thermal conduction, static
        TransientHeatConduction,    // thermal conduction, transient
        Gradient2D,                 // computes derivatives in volume
        Gradient3D,                 // computes derivatives in volume
        SurfaceGradient,            // computes derivatives on surface
        PlaneStress,                // plane stress in 2D
        LinearElasticity,           // Linear Elasticity
        TATCAD,                     // aerothermodymanics

        MAXWELL_HA_TRI3,
        MAXWELL_HA_TET4,
        MAXWELL_HA_TRI6,
        MAXWELL_HA_TET10,

        MAXWELL_HPHI_TRI3,
        MAXWELL_HPHI_TET4,
        MAXWELL_HPHI_TRI6,
        MAXWELL_HPHI_TET10,

        MAXWELL_L2J_TRI3,
        MAXWELL_L2J_TET4,
        MAXWELL_L2J_TRI6,
        MAXWELL_L2J_TET10,

        MAXWELL_HA_L2B_TRI3,
        MAXWELL_HA_L2B_TET4,
        MAXWELL_HA_L2B_TRI6,
        MAXWELL_HA_L2B_TET10,

        MAXWELL_HPHI_L2B_TRI3,
        MAXWELL_HPHI_L2B_TET4,
        MAXWELL_HPHI_L2B_TRI6,
        MAXWELL_HPHI_L2B_TET10,

        MAXWELL_L2,

        UNDEFINED // ADD NEW IWGS TO is_maxwell AS WELL !!!
    };

//------------------------------------------------------------------------------

    inline bool
    is_maxwell_ha( const IwgType aType )
    {
        switch( aType )
        {
            case( IwgType::MAXWELL_HA_TRI3 ) :
            case( IwgType::MAXWELL_HA_TRI6  ) :
            case( IwgType::MAXWELL_HA_TET4 ) :
            case( IwgType::MAXWELL_HA_TET10  ) :
            {
                return true;
            }
            default :
            {
                return false;
            }
        }
    }

//------------------------------------------------------------------------------

    inline bool
    is_maxwell_hphi( const IwgType aType )
    {
        switch( aType )
        {
            case( IwgType::MAXWELL_HPHI_TRI3 ) :
            case( IwgType::MAXWELL_HPHI_TRI6  ) :
            case( IwgType::MAXWELL_HPHI_TET4 ) :
            case( IwgType::MAXWELL_HPHI_TET10  ) :
            {
                return true;
            }
            default :
            {
                return false;
            }
        }
    }

//------------------------------------------------------------------------------

    inline bool
    is_maxwell_postproc( const IwgType aType )
    {
        switch( aType )
        {
            case( IwgType::MAXWELL_L2 ) :
            case( IwgType::MAXWELL_L2J_TRI3 ) :
            case( IwgType::MAXWELL_L2J_TET4 ) :
            case( IwgType::MAXWELL_L2J_TRI6 ) :
            case( IwgType::MAXWELL_L2J_TET10 ) :

            case( IwgType::MAXWELL_HA_L2B_TRI3 ) :
            case( IwgType::MAXWELL_HA_L2B_TET4 ) :
            case( IwgType::MAXWELL_HA_L2B_TRI6 ) :
            case( IwgType::MAXWELL_HA_L2B_TET10 ) :

            case( IwgType::MAXWELL_HPHI_L2B_TRI3 ) :
            case( IwgType::MAXWELL_HPHI_L2B_TET4 ) :
            case( IwgType::MAXWELL_HPHI_L2B_TRI6 ) :
            case( IwgType::MAXWELL_HPHI_L2B_TET10 ) :
            {
                return true;
            }
            default :
            {
                return false;
            }
        }
    }

//------------------------------------------------------------------------------

    inline bool
    is_maxwell( const IwgType aType )
    {
        return
           is_maxwell_ha( aType )
        || is_maxwell_hphi( aType )
        || is_maxwell_postproc( aType );
    }

//------------------------------------------------------------------------------

    enum class IwgMode
    {
        Direct         =1,
        Iterative      =2,
        UNDEFINED      =3
    };

//------------------------------------------------------------------------------

    enum class SolverAlgorithm
    {
        Direct        = 0,
        NewtonRaphson = 1,
        Picard        = 2,
        UNDEFINED     = 3
    };

//------------------------------------------------------------------------------
}

#endif //BELFEM_EN_IWGS_HPP
