/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_EN_IWGS_HPP
#define BELFEM_EN_IWGS_HPP

namespace belfem
{
//------------------------------------------------------------------------------

    enum class ModelDimensionality
    {
        TwoD,
        AxSymmX,
        AxSymmY,
        ThreeD,
        UNDEFINED
    };

    //
    enum class Nfunction
    {
        Scalar = 0,
        Vec2d  = 1,
        Vec3d  = 2,
        UNDEFINED = 3
    };

    enum class Bfunction
    {
        Gradient     = 0,
        Planestress  = 1,
        Voigt        = 2,
        UNDEFINED = 3
    };


    /*enum class Ntype
    {
        Scalar,
        Vector2d,
        Vector3d,
        UNDEFINED
    };

    enum class Btype
    {
        Gradient2d,
        Gradient3d,
        Planestress,
        Voigt,
        UNDEFINED
    }; */

    enum class IwgType
    {
        Poisson,
        StaticHeatConduction,       // thermal conduction, static
        TransientHeatConduction,    // thermal conduction, transient
        Gradient2D,                 // computes derivatives in volume
        Gradient3D,                 // computes derivatives in volume
        SurfaceGradient,            // computes derivatives on surface
        PlaneStress,                // plane stress in 2D
        LinearElasticity,           // Linear Elasticity
        Maxwell,
        MaxwellThermal,
        UNDEFINED // ADD NEW IWGS TO is_maxwell AS WELL !!!
    };

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

    inline bool is_maxwell( const IwgType aType )
    {
        return aType == IwgType::Maxwell || aType == IwgType::MaxwellThermal ;
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_EN_IWGS_HPP
