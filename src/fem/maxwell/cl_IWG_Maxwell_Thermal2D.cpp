//
// Created by christian on 2/9/23.
//
#include "cl_IWG_Maxwell_Thermal2D.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_Thermal2D::IWG_Maxwell_Thermal2D( const uint aNumberOfDimensions ) :
                IWG( IwgType::TransientHeatConduction,
                     IwgMode::Iterative,
                     SymmetryMode::PositiveDefiniteSymmetric )
        {
            mNumberOfDofsPerNode = 1;
            mNumberOfDerivativeDimensions = aNumberOfDimensions;
            mNumberOfRhsCols = 1;

            // set the names for the fields
            mDofFields = { "T" };
            mFluxFields = { "ej" };
            mOtherFields = { "b", "T0" };

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_Thermal2D::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->print_dofs( aElement );
            exit( 0 );
        }

//------------------------------------------------------------------------------

    }
}