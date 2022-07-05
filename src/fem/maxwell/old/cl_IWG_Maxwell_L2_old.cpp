//
// Created by Christian Messe on 03.02.22.
//

#include "cl_IWG_Maxwell_L2_old.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_unique.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_L2::IWG_Maxwell_L2(
                const ElementType aElementType,
                        const Maxwell_L2_Mode aMode ) :
                IWG_Maxwell( aElementType,
                             IwgType::MAXWELL_L2,
                             IwgMode::Direct,
                             SymmetryMode::PositiveDefiniteSymmetric,
                             SideSetDofLinkMode::FacetOnly,
                             true )
        {
            switch( aElementType )
            {
                case( ElementType::TRI3 ) :
                {
                    this->set_mode_tri3( aMode );
                    break ;
                }
                case( ElementType::TRI6 ) :
                {
                    this->set_mode_tri6( aMode );
                    break ;
                }
                case( ElementType::TET4 ) :
                {
                    this->set_mode_tet4( aMode );
                    break ;
                }
                case( ElementType::TET10 ) :
                {
                    this->set_mode_tet10( aMode );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unknown Element Type");
                }
            }
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_L2::~IWG_Maxwell_L2()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::set_mode_tri3( const Maxwell_L2_Mode aMode )
        {
            switch( aMode )
            {
                case( Maxwell_L2_Mode::B2H ) :
                {
                    mFields.Superconductor = { "edge_h" };
                    mFields.MagneticFieldDensity = { "bx", "by" };

                    mEdgeDofMultiplicity = 1 ;
                    mFaceDofMultiplicity = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_edge ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_b2h ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2B ) :
                {
                    mFields.Superconductor = { "bx", "by" };
                    // ´mFields.MagneticFieldDensity = { "edge_h" };
                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tri3_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2b ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2J ) :
                {
                    mFields.Superconductor = { "jz" };
                    //mFields.MagneticFieldDensity = { "edge_h" };

                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tri3_scalar ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2jz ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unknown L2 mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::set_mode_tri6( const Maxwell_L2_Mode aMode )
        {
            switch( aMode )
            {
                case( Maxwell_L2_Mode::B2H ) :
                {
                    mFields.Superconductor = { "edge_h", "face_h" };
                    mFields.MagneticFieldDensity = { "bx", "by" };

                    mEdgeDofMultiplicity = 2 ;
                    mFaceDofMultiplicity = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_edge ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_b2h ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2B ) :
                {
                    mFields.Superconductor = { "bx", "by" };
                    //mFields.MagneticFieldDensity = { "edge_h", "face_h" };

                    mNumberOfRhsDofsPerEdge = 2 ;
                    mNumberOfRhsDofsPerFace = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tri6_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2b ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2J ) :
                {
                    mFields.Superconductor = { "jz" };
                   // mFields.MagneticFieldDensity = { "edge_h", "face_h" };

                    mNumberOfRhsDofsPerEdge = 2 ;
                    mNumberOfRhsDofsPerFace = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tri6_scalar ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2jz ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unknown L2 mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::set_mode_tet4( const Maxwell_L2_Mode aMode )
        {
            switch( aMode )
            {
                case( Maxwell_L2_Mode::B2H ) :
                {
                    mFields.Superconductor = { "edge_h" };
                    mFields.MagneticFieldDensity = { "bx", "by", "bz" };

                    mEdgeDofMultiplicity = 1 ;
                    mFaceDofMultiplicity = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_edge ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_b2h ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2B ) :
                {
                    mFields.Superconductor = { "bx", "by", "bz" };
                    //mFields.MagneticFieldDensity = { "edge_h" };

                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tet4_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2b ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2J ) :
                {
                    mFields.Superconductor = { "jx", "jy", "jz" };
                    //mFields.MagneticFieldDensity = { "edge_h" };

                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tet4_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2j ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unknown L2 mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::set_mode_tet10( const Maxwell_L2_Mode aMode )
        {
            switch( aMode )
            {
                case( Maxwell_L2_Mode::B2H ) :
                {
                    mFields.Superconductor = { "edge_h", "face_h" };
                    mFields.MagneticFieldDensity = { "bx", "by", "bz" };

                    mEdgeDofMultiplicity = 2 ;
                    mFaceDofMultiplicity = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_edge ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_b2h ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2B ) :
                {
                    mFields.Superconductor = { "bx", "by", "bz" };
                    //mFields.MagneticFieldDensity = { "edge_h", "face_h" };

                    mNumberOfRhsDofsPerEdge = 2 ;
                    mNumberOfRhsDofsPerFace = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tet10_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2b ;
                    break ;
                }
                case( Maxwell_L2_Mode::H2J ) :
                {
                    mFields.Superconductor = { "jx", "jy", "jz" };
                    //mFields.MagneticFieldDensity = { "edge_h", "face_h"};

                    mNumberOfRhsDofsPerEdge = 2 ;
                    mNumberOfRhsDofsPerFace = 2 ;

                    mFunJacobian = & IWG_Maxwell_L2::l2_node_tet10_vector ;
                    mFunRHS      = & IWG_Maxwell_L2::l2_h2j ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unknown L2 mode");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::set_block(
                const id_t aBlockID,
                const DomainType aDomainType )
        {
            Vector< id_t > tBlockIDs( mBlockIDs );
            mBlockIDs.set_size( tBlockIDs.length() + 1 );
            uint tCount = 0 ;
            for( id_t tID : tBlockIDs )
            {
                mBlockIDs( tCount++ ) = tID ;
            }
            mBlockIDs( tCount++ ) = aBlockID ;
            unique( mBlockIDs );

            mBlockTypes[ aBlockID ] = aDomainType;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::link_jacobian_function( Group * aGroup )
        {
            mEdgeFunction->link( aGroup );
            BELFEM_ERROR( aGroup->domain_type() == DomainType::SuperConductor,
                         "The L2 projector can only be connected to a superconducting group");

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::allocate_work_matrices( Group    * aGroup )
        {
            mGroup->node_coords().set_size(
                mesh::number_of_nodes( aGroup->element_type() ),
                mMesh->number_of_dimensions() );

            switch ( aGroup->element_type() )
            {
                case ( ElementType::TRI3 ) :
                {
                    // precomputed solution for ∫ N'N dxi, scalar field
                    mGroup->work_Sigma() = {
                            { 2., 1., 1. },
                            { 1., 2., 1. },
                            { 1., 1., 2. } };
                    mGroup->work_Sigma() /= 24.0;

                    // precomputed solution for ∫ N'N dxi, vector field
                    mGroup->work_Tau() = {
                            { 2., 0., 1., 0., 1., 0. },
                            { 0., 2., 0., 1., 0., 1. },
                            { 1., 0., 2., 0., 1., 0. },
                            { 0., 1., 0., 2., 0., 1. },
                            { 1., 0., 1., 0., 2., 0. },
                            { 0., 1., 0., 1., 0., 2. } };
                    mGroup->work_Tau() /= 24.0;

                    mGroup->work_J().set_size( 2, 2 );

                    // contain jz as node data
                    mGroup->work_phi().set_size( 3 );
                    mGroup->work_Phi().set_size( 3, 2 );

                    // vector field
                    mGroup->work_chi().set_size( 2 );

                    // nedelec data
                    mGroup->work_psi().set_size( 3 );
                    break;
                }
                case ( ElementType::TRI6 ) :
                {
                    // precomputed solution for ∫ N'N dxi, straight edge only, scalar field
                    mGroup->work_Sigma() = { { 6., -1., -1., 0., -4., 0. },
                                             { -1., 6., -1., 0., 0., -4. },
                                             { -1., -1., 6., -4., 0., 0. },
                                             { 0., 0., -4., 32., 16., 16. },
                                             { -4., 0., 0., 16., 32., 16. },
                                             { 0., -4., 0., 16., 16., 32. } };
                    mGroup->work_Sigma() /= 360.0 ;

                    // precomputed solution for ∫ N'N dxi, straight edge only, vector field
                    mGroup->work_Tau() = { { 6., 0., -1., 0., -1., 0., 0., 0., -4., 0., 0., 0. },
                                             { 0., 6., 0., -1., 0., -1., 0., 0., 0., -4., 0., 0. },
                                             { -1., 0., 6., 0., -1., 0., 0., 0., 0., 0., -4., 0. },
                                             { 0., -1., 0., 6., 0., -1., 0., 0., 0., 0., 0., -4. },
                                             { -1., 0., -1., 0., 6., 0., -4., 0., 0., 0., 0., 0. },
                                             { 0., -1., 0., -1., 0., 6., 0., -4., 0., 0., 0., 0. },
                                             { 0., 0., 0., 0., -4., 0., 32., 0., 16., 0., 16., 0. },
                                             { 0., 0., 0., 0., 0., -4., 0., 32., 0., 16., 0., 16. },
                                             { -4., 0., 0., 0., 0., 0., 16., 0., 32., 0., 16., 0. },
                                             { 0., -4., 0., 0., 0., 0., 0., 16., 0., 32., 0., 16. },
                                             { 0., 0., -4., 0., 0., 0., 16., 0., 16., 0., 32., 0. },
                                             { 0., 0., 0., -4., 0., 0., 0., 16., 0., 16., 0., 32. } };
                    mGroup->work_Tau() /= 360.0;

                    // contains jz
                    mGroup->work_phi().set_size( 6 );

                    // contains b
                    mGroup->work_Phi().set_size( 6, 2 );

                    // vector field
                    mGroup->work_chi().set_size( 2 );

                    // nedelec data
                    mGroup->work_psi().set_size( 8 );

                    mGroup->work_J().set_size( 2,2 );

                    break;
                }
                case ( ElementType::TET4 ) :
                {
                    // precomputed solution for ∫ N'N dxi, straight edge only, scalar field
                    mGroup->work_Sigma() = { { 2., 1., 1., 1. },
                                             { 1., 2., 1., 1. },
                                             { 1., 1., 2., 1. },
                                             { 1., 1., 1., 2. } };
                    mGroup->work_Sigma() /= 120.0 ;

                    // precomputed solution for ∫ N'N dxi, vector field
                    mGroup->work_Tau() = {
                            { 2., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0. },
                            { 0., 2., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. },
                            { 0., 0., 2., 0., 0., 1., 0., 0., 1., 0., 0., 1. },
                            { 1., 0., 0., 2., 0., 0., 1., 0., 0., 1., 0., 0. },
                            { 0., 1., 0., 0., 2., 0., 0., 1., 0., 0., 1., 0. },
                            { 0., 0., 1., 0., 0., 2., 0., 0., 1., 0., 0., 1. },
                            { 1., 0., 0., 1., 0., 0., 2., 0., 0., 1., 0., 0. },
                            { 0., 1., 0., 0., 1., 0., 0., 2., 0., 0., 1., 0. },
                            { 0., 0., 1., 0., 0., 1., 0., 0., 2., 0., 0., 1. },
                            { 1., 0., 0., 1., 0., 0., 1., 0., 0., 2., 0., 0. },
                            { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 2., 0. },
                            { 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 2. } };
                    mGroup->work_Tau() /= 120.0;
                    mGroup->work_J().set_size( 3, 3 );

                    // contains b or j
                    mGroup->work_Phi().set_size( 4, 3 );

                    // vector field
                    mGroup->work_chi().set_size( 3 );

                    // nedelec data
                    mGroup->work_psi().set_size( 6 );

                    break;
                }
                case ( ElementType::TET10 ) :
                {
                    // precomputed solution for ∫ N'N dxi, straight edge only, scalar field
                    mGroup->work_Sigma() = { { 6., 1., 1., 1., -4., -6., -4., -4., -6., -6. },
                                             { 1., 6., 1., 1., -4., -4., -6., -6., -4., -6. },
                                             { 1., 1., 6., 1., -6., -4., -4., -6., -6., -4. },
                                             { 1., 1., 1., 6., -6., -6., -6., -4., -4., -4. },
                                             { -4., -4., -6., -6., 32., 16., 16., 16., 16., 8. },
                                             { -6., -4., -4., -6., 16., 32., 16., 8., 16., 16. },
                                             { -4., -6., -4., -6., 16., 16., 32., 16., 8., 16. },
                                             { -4., -6., -6., -4., 16., 8., 16., 32., 16., 16. },
                                             { -6., -4., -6., -4., 16., 16., 8., 16., 32., 16. },
                                             { -6., -6., -4., -4., 8., 16., 16., 16., 16., 32. } };
                    mGroup->work_Sigma() /= 2520.0 ;

                    // precomputed solution for ∫ N'N dxi, straight edge only, vector field
                    mGroup->work_Tau() = { {  6.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.,
                                                0.,  0., -4.,  0.,  0., -6.,  0.,  0., -4.,  0.,
                                                0., -4.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0.  },
                                             {  0.,  6.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,
                                                1.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -4.,
                                                0.,  0., -4.,  0.,  0., -6.,  0.,  0., -6.,  0.  },
                                             {  0.,  0.,  6.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,
                                                0.,  1.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0.,
                                               -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -6.  },
                                             {  1.,  0.,  0.,  6.,  0.,  0.,  1.,  0.,  0.,  1.,
                                                0.,  0., -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,
                                                0., -6.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0.  },
                                             {  0.,  1.,  0.,  0.,  6.,  0.,  0.,  1.,  0.,  0.,
                                                1.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0., -6.,
                                                0.,  0., -6.,  0.,  0., -4.,  0.,  0., -6.,  0.  },
                                             {  0.,  0.,  1.,  0.,  0.,  6.,  0.,  0.,  1.,  0.,
                                                0.,  1.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0.,
                                               -6.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -6.  },
                                             {  1.,  0.,  0.,  1.,  0.,  0.,  6.,  0.,  0.,  1.,
                                                0.,  0., -6.,  0.,  0., -4.,  0.,  0., -4.,  0.,
                                                0., -6.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0.  },
                                             {  0.,  1.,  0.,  0.,  1.,  0.,  0.,  6.,  0.,  0.,
                                                1.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -4.,
                                                0.,  0., -6.,  0.,  0., -6.,  0.,  0., -4.,  0.  },
                                             {  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  6.,  0.,
                                                0.,  1.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0.,
                                                -4.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0., -4.  },
                                             {  1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  6.,
                                                0.,  0., -6.,  0.,  0., -6.,  0.,  0., -6.,  0.,
                                                0., -4.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0.  },
                                             {  0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,
                                                6.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0., -6.,
                                                0.,  0., -4.,  0.,  0., -4.,  0.,  0., -4.,  0.  },
                                             {  0.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  1.,  0.,
                                                0.,  6.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0.,
                                               -6.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0., -4.  },
                                             { -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -6.,
                                                0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,
                                                0.,  16.,  0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.  },
                                             {  0., -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -6.,
                                                0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,
                                                0.,  16.,  0.,  0.,  16.,  0.,  0.,  8.,  0.  },
                                             {  0.,  0., -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,
                                                0., -6.,  0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,
                                                16.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  8.  },
                                             { -6.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0., -6.,  0.,
                                                0.,  16.,  0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,
                                                8.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.  },
                                             {  0., -6.,  0.,  0., -4.,  0.,  0., -4.,  0.,  0.,
                                               -6.,  0.,  0.,  16.,  0.,  0.,  32.,  0.,  0.,  16.,
                                                0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.,  16.,  0.  },
                                             {  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -4.,  0.,
                                                0., -6.,  0.,  0.,  16.,  0.,  0.,  32.,  0.,  0.,
                                                16.,  0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.,  16.  },
                                             { -4.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -6.,
                                                0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  32.,  0.,
                                                0.,  16.,  0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.  },
                                             {  0., -4.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0.,
                                               -6.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  32.,
                                                0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.,  16.,  0.  },
                                             {  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -4.,  0.,
                                                0., -6.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,
                                                32.,  0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.,  16.  },
                                             { -4.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0., -4.,
                                                0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.,  16.,  0.,
                                                0.,  32.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.  },
                                             {  0., -4.,  0.,  0., -6.,  0.,  0., -6.,  0.,  0.,
                                                -4.,  0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.,  16.,
                                                0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,  16.,  0.  },
                                             {  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -6.,  0.,
                                                0., -4.,  0.,  0.,  16.,  0.,  0.,  8.,  0.,  0.,
                                                16.,  0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.,  16.  },
                                             { -6.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0., -4.,
                                                0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  8.,  0.,
                                                0.,  16.,  0.,  0.,  32.,  0.,  0.,  16.,  0.,  0.  },
                                             {  0., -6.,  0.,  0., -4.,  0.,  0., -6.,  0.,  0.,
                                               -4.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  8.,
                                                0.,  0.,  16.,  0.,  0.,  32.,  0.,  0.,  16.,  0.  },
                                             {  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -6.,  0.,
                                                0., -4.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,
                                                8.,  0.,  0.,  16.,  0.,  0.,  32.,  0.,  0.,  16.  },
                                             { -6.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0., -4.,
                                                0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,
                                                0.,  16.,  0.,  0.,  16.,  0.,  0.,  32.,  0.,  0.  },
                                             {  0., -6.,  0.,  0., -6.,  0.,  0., -4.,  0.,  0.,
                                               -4.,  0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.,  16.,
                                                0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  32.,  0.  },
                                             {  0.,  0., -6.,  0.,  0., -6.,  0.,  0., -4.,  0.,
                                                0., -4.,  0.,  0.,  8.,  0.,  0.,  16.,  0.,  0.,
                                                16.,  0.,  0.,  16.,  0.,  0.,  16.,  0.,  0.,  32.  } };
                    mGroup->work_Tau() /= 2520.0 ;
                    mGroup->work_J().set_size( 3, 3 );

                    // contains b or j
                    mGroup->work_Phi().set_size( 10, 3 );

                    // vector field
                    mGroup->work_chi().set_size( 3 );

                    // nedelec data
                    mGroup->work_psi().set_size( 20 );

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type");
                }
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tri3_scalar(
                 Element * aElement,
                 Matrix< real > & aK )
        {
           aK = mGroup->work_Sigma() *
                    std::abs( det( this->J_tri_straight( aElement ) ) );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tri3_vector(
                Element * aElement,
                Matrix< real > & aK )
        {
            aK = mGroup->work_Tau() *
                 std::abs( det( this->J_tri_straight( aElement ) ) );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tri6_scalar(
                Element * aElement,
                Matrix< real > & aK )
        {
            if( aElement->element()->is_curved() )
            {
                aElement->get_node_coors( mGroup->node_coords() );
                aK.fill( 0.0 );
                const Vector< real > & tW = mGroup->integration_weights() ;
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                    const Matrix< real > & tN = mGroup->N( k );

                    aK += ( tW( k ) * std::abs( det( tJ ) ) ) *
                          trans( tN ) * tN ;
                }
            }
            else
            {
                aK = mGroup->work_Sigma()
                     * std::abs( det( this->J_tri_straight( aElement ) ) );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tri6_vector(
                Element * aElement,
                Matrix< real > & aK )
        {
            if( aElement->element()->is_curved() )
            {
                aElement->get_node_coors( mGroup->node_coords() );
                aK.fill( 0.0 );
                const Vector< real > & tW = mGroup->integration_weights() ;
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                    const Matrix< real > & tN = mGroup->Nvector( k );

                    aK += ( tW( k ) * std::abs( det( tJ ) ) ) *
                            trans( tN ) * tN ;
                }
            }
            else
            {
                aK = mGroup->work_Tau()
                        * std::abs( det( this->J_tri_straight( aElement ) ) );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tet4_scalar(
                Element * aElement,
                Matrix< real > & aK )
        {
            aK = mGroup->work_Sigma()
                 * std::abs( det( this->J_tet_straight( aElement ) ) );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tet4_vector(
                Element * aElement,
                Matrix< real > & aK )
        {
            aK = mGroup->work_Tau()
                            * std::abs( det( this->J_tet_straight( aElement ) ) );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tet10_scalar(
                Element * aElement,
                Matrix< real > & aK )
        {
            if( aElement->element()->is_curved() )
            {
                aElement->get_node_coors( mGroup->node_coords() );
                aK.fill( 0.0 );
                const Vector< real > & tW = mGroup->integration_weights() ;
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                    const Matrix< real > & tN = mGroup->N( k );

                    aK += ( tW( k ) * std::abs( det( tJ ) ) )
                          * trans( tN ) * tN ;
                }
            }
            else
            {
                aK = mGroup->work_Sigma() *
                     std::abs( det( this->J_tet_straight( aElement ) ) );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_node_tet10_vector(
                Element * aElement,
                Matrix< real > & aK )
        {
            if( aElement->element()->is_curved() )
            {
                aElement->get_node_coors( mGroup->node_coords() );
                aK.fill( 0.0 );
                const Vector< real > & tW = mGroup->integration_weights() ;
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                    const Matrix< real > & tN = mGroup->Nvector( k );

                    aK += ( tW( k ) * std::abs( det( tJ ) ) )
                            * trans( tN ) * tN ;
                }
            }
            else
            {
                aK = mGroup->work_Tau() *
                        std::abs( det( this->J_tet_straight( aElement ) ) );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_edge(
                Element * aElement,
                Matrix< real > & aK )
        {
            mEdgeFunction->link( aElement, true, false, false );

            const Vector< real > & tW = mGroup->integration_weights() ;

            aK.fill( 0.0 );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tE = mEdgeFunction->E( k );
                aK += ( tW( k ) * mEdgeFunction->abs_det_J() ) *
                        trans( tE ) * tE ;
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_L2::J_tri_straight( Element * aElement )
        {
            mesh::Element * tElement = aElement->element() ;

            Matrix< real > & aJ = mGroup->work_J() ;
            aJ( 0, 0 ) =  tElement->node( 0 )->x() - tElement->node( 2 )->x() ;
            aJ( 1, 0 ) =  tElement->node( 1 )->x() - tElement->node( 2 )->x() ;
            aJ( 0, 1 ) =  tElement->node( 0 )->y() - tElement->node( 2 )->y() ;
            aJ( 1, 1 ) =  tElement->node( 1 )->y() - tElement->node( 2 )->y() ;

            return aJ ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_L2::J_tet_straight( Element * aElement )
        {
            mesh::Element * tElement = aElement->element() ;

            Matrix< real > & aJ = mGroup->work_J() ;
            aJ( 0, 0 ) =  tElement->node( 0 )->x() - tElement->node( 3 )->x() ;
            aJ( 1, 0 ) =  tElement->node( 1 )->x() - tElement->node( 3 )->x() ;
            aJ( 2, 0 ) =  tElement->node( 2 )->x() - tElement->node( 3 )->x() ;

            aJ( 0, 1 ) =  tElement->node( 0 )->y() - tElement->node( 3 )->y() ;
            aJ( 1, 1 ) =  tElement->node( 1 )->y() - tElement->node( 3 )->y() ;
            aJ( 2, 1 ) =  tElement->node( 2 )->y() - tElement->node( 3 )->y() ;

            aJ( 0, 2 ) =  tElement->node( 0 )->z() - tElement->node( 3 )->z() ;
            aJ( 1, 2 ) =  tElement->node( 1 )->z() - tElement->node( 3 )->z() ;
            aJ( 2, 2 ) =  tElement->node( 2 )->z() - tElement->node( 3 )->z() ;

            return aJ ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_b2h( Element * aElement, Vector< real > & aRHS )
        {
            mEdgeFunction->link( aElement, true,
                                 false, false );

            const Vector< real > & tW = mGroup->integration_weights() ;
            Matrix< real > & tB = mGroup->work_Phi() ;

            this->collect_node_data( aElement, mFields.MagneticFieldDensity,tB );
            aRHS.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tE = mEdgeFunction->E( k );

                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->N( k ) * tB  * tE );

            }
            // todo: uncomment  this
            //aRHS *= mMaterial->nu_s();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_h2b( Element * aElement, Vector< real > & aRHS )
        {
            mEdgeFunction->link( aElement, true,
                                 false, false );

            const Vector< real > & tW = mGroup->integration_weights() ;
            Vector< real > & tH = mGroup->work_psi();

            this->collect_nedelec_data( aElement, NedelecField::H, tH );
            aRHS.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tE = mEdgeFunction->E( k );


                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->Nvector( k ) ) * tE * tH ;


            }
            // todo: uncomment  this
            //aRHS *= constant::mu0 * mMaterial->mu_r();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_h2jz( Element * aElement, Vector< real > & aRHS )
        {
            mEdgeFunction->link( aElement, true,
                                 false, false );

            const Vector< real > & tW = mGroup->integration_weights() ;
            Vector< real > & tH = mGroup->work_psi();

            this->collect_nedelec_data( aElement, NedelecField::H, tH );
            aRHS.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tC = mEdgeFunction->C( k );

                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->N( k ) ) * tC * tH ;

            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2::l2_h2j( Element * aElement, Vector< real > & aRHS )
        {
            mEdgeFunction->link( aElement, true,
                                 false, false );

            const Vector< real > & tW = mGroup->integration_weights() ;
            Vector< real > & tH = mGroup->work_psi();

            this->collect_nedelec_data( aElement, NedelecField::H, tH );
            aRHS.fill( 0.0 );
            Vector< real > & tJ = mGroup->work_chi() ;

            Vector< real > & tX = mMesh->field_data( "ElJx");
            Vector< real > & tY = mMesh->field_data( "ElJy");
            Vector< real > & tZ = mMesh->field_data( "ElJz");
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tC = mEdgeFunction->C( k );

                // todo: remove this maybe
                tJ = tC * tH ;
                if( k == 0 )
                {
                    tX( aElement->element()->index()) = tJ( 0 );
                    tY( aElement->element()->index()) = tJ( 1 );
                    tZ( aElement->element()->index()) = tJ( 2 );
                }


                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->Nvector( k ) ) * tJ ;

            }
        }

//------------------------------------------------------------------------------
    }
}