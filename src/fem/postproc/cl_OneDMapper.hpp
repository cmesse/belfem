//
// Created by Christian Messe on 30.12.20.
//

#ifndef BELFEM_CL_ONEDMAPPER_HPP
#define BELFEM_CL_ONEDMAPPER_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Mesh.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_Matrix.hpp"
#include "cl_Solver.hpp"
namespace belfem
{
    class OneDMapper
    {
        const index_t mNumTargetNodes ;
        const uint    mOrder ;
        const proc_t  mMasterRank ;
        const proc_t  mMyRank ;

        // the type of these elements
        ElementType mElementType = ElementType::EMPTY ;

        Cell< mesh::Node * >    mNodes ;
        Cell< mesh::Element * > mElements ;

        SpMatrix  * mMassMatrix = nullptr ;

        fem::InterpolationFunction * mShape = nullptr ;

        // the evaluated shape function
        Cell< Matrix< real > > mN ;

        // the evaluated shape function derivative
        Cell< Matrix< real > > mdNdXi ;

        // the integraiton points
        Matrix< real > mPoints ;

        // the integration weights
        Vector< real > mWeights ;

        // the solver
        Solver * mSolver = nullptr ;

        // the right hand side
        Vector< real > mRHS ;

//------&------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        OneDMapper( const Vector< real > & aTargetCoords,
                    const uint             aOrder = 1,
                    const proc_t           aMasterRank = 0 );

//------------------------------------------------------------------------------

        ~OneDMapper();

//------------------------------------------------------------------------------

        void
        project( const Vector< real > & aSourceNodes,
                 const Vector< real > & aSourceValues,
                       Vector< real > & aTargetValues );

//------------------------------------------------------------------------------

        void
        derive( const Vector< real > & aValues,
                      Vector< real > & aDerivatives );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        create_target_nodes( const Vector< real > & aTargetCoords );

//------------------------------------------------------------------------------

        void
        create_elements();

//------------------------------------------------------------------------------

        void
        compute_adjency();

//------------------------------------------------------------------------------

        void
        create_shape_function();

//------------------------------------------------------------------------------

        void
        create_mass_matrix();

//------------------------------------------------------------------------------

        index_t
        find_interval( const Vector< real > & aSourceNodes, const real aX );

//------------------------------------------------------------------------------

        void
        compute_rhs( const Vector< real > & aSourceNodes,
                     const Vector< real > & aSourceValues,
                           Vector< real > & aRHS );

//------------------------------------------------------------------------------

        void
        compute_derivative_rhs( const Vector< real > & aValues,
                                Vector< real > & aRHS );

//------------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_ONEDMAPPER_HPP
