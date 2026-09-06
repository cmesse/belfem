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

#include "globals.hpp"
#include "commtools.hpp"
#include "cl_MaxwellPostprocessor.hpp"
#include "cl_FEM_Controller.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_posv.hpp"
#include "fn_combine.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace fem
    {
        MaxwellPostprocessor::MaxwellPostprocessor(
            Kernel * aKernel,
            const Map< id_t, DomainType > & aBlockTypes,
            const Map< id_t, string >     & aMaterialMap,
            const MaxwellPostprocessorType aType,
            const bool aCreateElementFields ) :
            Postprocessor( aKernel ),
            mType( aType ),
            mCreateElementFields( aCreateElementFields )
        {
            mDomainType = this->select_domain_type( aType );
            this->set_type( aType );
            this->select_blocks_and_materials( aBlockTypes, aMaterialMap );

        }

        MaxwellPostprocessor::~MaxwellPostprocessor()
        {
            mMaterialMap.clear() ;

        }

        void
        MaxwellPostprocessor::set_type( const MaxwellPostprocessorType aType )
        {
            switch ( aType )
            {
                case MaxwellPostprocessorType::Air :
                {
                    mUpdateFunction = & MaxwellPostprocessor::update_dofs_lagrange ;
                    mComputeFunction = & MaxwellPostprocessor::compute_air ;

                    this->set_source_fields({ "phi" } );

                    if ( mNumDimensions == 2 )
                    {
                        this->set_target_fields( { "Hx", "Hy" } );
                    }
                    else
                    {
                        this->set_target_fields( { "Hx", "Hy", "Hz" } );
                    }
                    break ;
                }
                case MaxwellPostprocessorType::Ferro :
                {
                    mUpdateFunction  = & MaxwellPostprocessor::update_dofs_lagrange ;
                    mComputeFunction = & MaxwellPostprocessor::compute_ferro ;

                    this->set_source_fields({ "phi" } );
                    if ( mNumDimensions == 2 )
                    {
                        this->set_target_fields( { "Hx", "Hy", "Bx", "By" } );
                    }
                    else
                    {
                        this->set_target_fields( { "Hx", "Hy", "Hz", "Bx", "By", "Bz" } );
                    }
                    break ;
                }
                case MaxwellPostprocessorType::Conductor :
                case MaxwellPostprocessorType::ThinShellConductor :
                case MaxwellPostprocessorType::SideConnector :
                {
                    mUpdateFunction = & MaxwellPostprocessor::update_dofs_nedelec ;

                    // the wall recovers the full field ( ht + hb + hn ) via
                    // MaxwellData; E*q alone would miss the normal and
                    // binomial components. Thin-shell layers recover the
                    // normal component from the air phi
                    mComputeFunction =
                          aType == MaxwellPostprocessorType::SideConnector ?
                            & MaxwellPostprocessor::compute_side_connector
                        : aType == MaxwellPostprocessorType::ThinShellConductor ?
                            & MaxwellPostprocessor::compute_conductor_ts
                        : & MaxwellPostprocessor::compute_conductor ;

                    if ( mMesh->field_exists( "face_h" ) )
                    {
                        this->set_source_fields( { "edge_h", "face_h" } );
                    }
                    else
                    {
                        this->set_source_fields({ "edge_h" } );
                    }
                    if ( mNumDimensions == 2 )
                    {
                        this->set_target_fields( { "Hx", "Hy", "Bx", "By", "Jz"} );
                    }
                    else
                    {
                       this->set_target_fields( { "Hx", "Hy", "Hz", "Bx", "By", "Bz", "Jx", "Jy", "Jz" } );
                    }

                    break ;
                }
                case MaxwellPostprocessorType::SuperConductor :
                case MaxwellPostprocessorType::ThinShellSuperConductor :
                {
                    mUpdateFunction  = & MaxwellPostprocessor::update_dofs_nedelec ;
                    if ( aType == MaxwellPostprocessorType::SuperConductor )
                    {
                        mComputeFunction = & MaxwellPostprocessor::compute_superconductor ; // <-- needs to be replaced through special function for TS
                    }
                    else if ( aType == MaxwellPostprocessorType::ThinShellSuperConductor )
                    {
                        mComputeFunction = & MaxwellPostprocessor::compute_superconductor_ts ; // <-- needs to be replaced through special function for TS
                    }

                    if ( mMesh->field_exists( "face_h" ) )
                    {
                       this->set_source_fields( { "edge_h", "face_h" } );
                    }
                    else
                    {
                        this->set_source_fields( { "edge_h" } );
                    }

                    if ( mNumDimensions == 2 )
                    {
                        this->set_target_fields( { "Hx", "Hy", "Bx", "By", "Jz", "JJCz" } );
                    }
                    else
                    {
                        this->set_target_fields({ "Hx", "Hy", "Hz", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "JJCx", "JJCy", "JJCz" } );
                    }

                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Undefined postprocessor type" );
                }
            }

        }

        void
        MaxwellPostprocessor::select_blocks_and_materials(
            const Map< id_t, DomainType > & aBlockTypes,
            const Map< id_t, string >     & aMaterialMap )
        {
            Cell< id_t > tBlockIDs;

            switch ( mType )
            {
                case MaxwellPostprocessorType::Air :
                {
                    index_t tCount = 0 ;
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if ( tPair.second == DomainType::Air || tPair.second == DomainType::Buffer )
                        {
                            ++tCount ;
                        }
                    }

                    tBlockIDs.set_size( tCount );
                    tCount = 0 ;
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if ( tPair.second == DomainType::Air || tPair.second == DomainType::Buffer )
                        {
                            tBlockIDs( tCount++ ) = tPair.first ;
                        }
                    }
                    break ;
                }
                case MaxwellPostprocessorType::Ferro :
                {
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if ( tPair.second == DomainType::Ferro )
                        {
                            mMaterialMap[ tPair.first ] = mKernel->material( aMaterialMap( tPair.first) );
                        }
                    }

                    tBlockIDs.set_size( mMaterialMap.size() );
                    index_t tCount = 0 ;
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if ( tPair.second == DomainType::Ferro )
                        {
                            tBlockIDs( tCount++ ) = tPair.first ;
                        }
                    }

                    break ;
                }
                case MaxwellPostprocessorType::SideConnector :
                {
                    // edge coating walls: pure metal by the factory gate, so
                    // no jc filter is needed
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if (   tPair.second == DomainType::LeftCoating
                            || tPair.second == DomainType::RightCoating )
                        {
                            mMaterialMap[ tPair.first ] = mKernel->material( aMaterialMap( tPair.first) );
                        }
                    }
                    tBlockIDs.set_size( mMaterialMap.size() );
                    index_t tCount = 0 ;
                    for ( auto & tPair : mMaterialMap )
                    {
                        tBlockIDs( tCount++ ) = tPair.first ;
                    }
                    break ;
                }
                case MaxwellPostprocessorType::Conductor :
                case MaxwellPostprocessorType::ThinShellConductor :
                {
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if (   tPair.second == DomainType::Conductor
                            || tPair.second == DomainType::ThinShell )
                        {
                            // grab material
                            Material * tMaterial = mKernel->material( aMaterialMap( tPair.first) );
                            if ( ! tMaterial->have(MaterialProperty::jc) )
                            {
                                mMaterialMap[ tPair.first ] = mKernel->material( aMaterialMap( tPair.first) );
                            }
                        }
                    }
                    tBlockIDs.set_size( mMaterialMap.size() );
                    index_t tCount = 0 ;
                    for ( auto & tPair : mMaterialMap )
                    {
                        tBlockIDs( tCount++ ) = tPair.first ;
                    }
                    break ;
                }
                case MaxwellPostprocessorType::SuperConductor :
                case MaxwellPostprocessorType::ThinShellSuperConductor :
                {
                    for ( auto & tPair : aBlockTypes )
                    {
                        if ( ! mKernel->mesh()->block_exists( tPair.first ) ) continue ;
                        if ( tPair.second == DomainType::Conductor || tPair.second == DomainType::ThinShell )
                        {
                            // grab material
                            Material * tMaterial = mKernel->material( aMaterialMap( tPair.first ) );
                            if ( tMaterial->have(MaterialProperty::jc) )
                            {
                                mMaterialMap[ tPair.first ] = mKernel->material( aMaterialMap( tPair.first) );
                            }
                        }
                    }
                    tBlockIDs.set_size( mMaterialMap.size() );
                    index_t tCount = 0 ;
                    for ( auto & tPair : mMaterialMap )
                    {
                        tBlockIDs( tCount++ ) = tPair.first ;
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Undefined postprocessor type" );
                }
            }

            this->set_block_ids( tBlockIDs );
        }

        void
        MaxwellPostprocessor::initialize()
        {
            Postprocessor::initialize();

            this->create_other_fields();

            mX.set_size( mNumDimensions );
            mX0.set_size( mNumDimensions );

            mH.set_size( mNumDimensions );
            mB.set_size( mNumDimensions );
            mJ.set_size( mNumDimensions ==  2 ? 1 : 3 );
            mJJc.set_size( mNumDimensions ==  2 ? 1 : 3 );

            mY.set_size( mNumTargetFields );
            mZ.set_size( mNumTargetFields );


            switch ( mType )
            {
                case MaxwellPostprocessorType::Ferro :
                case MaxwellPostprocessorType::Conductor :
                case MaxwellPostprocessorType::ThinShellConductor :
                case MaxwellPostprocessorType::SuperConductor :
                case MaxwellPostprocessorType::ThinShellSuperConductor :
                {
                    if ( mCreateElementFields )
                    {
                        mComputeElementFields = true ;
                        // allocate mElementIndices
                        this->select_owned_elements();
                        this->create_element_fields() ;
                    }

                    break ;
                }
                default:
                {
                    break ;
                }
            }
        }

        void
        MaxwellPostprocessor::create_other_fields()
        {

            // note: even for 2D, we create also the Hz and Bz to make ParaView understand
            Cell< string > tOtherFields = { "Bx", "By", "Hz", "Bz" };
            for ( const string & tField : tOtherFields )
            {
                if ( ! mMesh->field_exists( tField ) )
                {
                    mMesh->create_field( tField );
                }
            }
        }

        void
        MaxwellPostprocessor::select_owned_elements()
        {
            const Cell< index_t > & tOwnedIndices = this->my_element_indices();

            mMyOwnedElementIndices.reserve( tOwnedIndices.size() );

            Cell< mesh::Element * > & tElements = mMesh->elements();
            for ( index_t e : tOwnedIndices )
            {
                if ( tElements( e )->owner() == mCommRank )
                {
                    mMyOwnedElementIndices.push( e );
                }
            }

            mMyOwnedElementIndices.shrink_to_fit() ;

            comm_barrier() ;
            if ( mCommRank == 0 )
            {
                collect( mAllOwnedElementIndices );
            }
            else
            {
                send( mMyOwnedElementIndices );
            }
        }

        void
        MaxwellPostprocessor::run()
        {
            Postprocessor::run();

            switch ( mType )
            {
                case MaxwellPostprocessorType::Air :
                {
                    Vector< real > & tHx = mMesh->field( "Hx" )->data();
                    Vector< real > & tHy = mMesh->field( "Hy" )->data();
                    Vector< real > & tBx = mMesh->field( "Bx" )->data();
                    Vector< real > & tBy = mMesh->field( "By" )->data();

                    tBx = constant::mu0 * tHx;
                    tBy = constant::mu0 * tHy;
                    if ( mMesh->number_of_dimensions() == 3 )
                    {
                        Vector< real > & tHz = mMesh->field( "Hz" )->data();
                        Vector< real > & tBz = mMesh->field( "Bz" )->data();
                        tBz = constant::mu0 * tHz;
                    }
                    break;
                }
                case MaxwellPostprocessorType::Conductor :
                case MaxwellPostprocessorType::ThinShellConductor :
                case MaxwellPostprocessorType::SuperConductor :
                case MaxwellPostprocessorType::ThinShellSuperConductor :
                case MaxwellPostprocessorType::Ferro :
                {
                    if ( mComputeElementFields )
                    {
                        this->compute_element_data() ;
                        this->collect_element_data() ;
                    }
                    break ;
                }
                case MaxwellPostprocessorType::SideConnector :
                {
                    this->copy_seam_fields() ;
                    break ;
                }
                default:
                {
                    break ;
                }
            }
        }

        void
        MaxwellPostprocessor::create_element_fields()
        {
            mElementTargetFields.set_size( mNumTargetFields, "");
            index_t tCount = 0 ;
            for ( string & tField : mTargetFields )
            {
                string tElementField = "element" + tField ;
                mElementTargetFields( tCount++ ) = tElementField ;
                if ( ! mMesh->field_exists( tElementField ) )
                {
                    mMesh->create_field( tElementField, EntityType::ELEMENT );
                }
            }
        }

        void
        MaxwellPostprocessor::update_dofs_lagrange()
        {
            mDOFs = mCalculator->node_data( "phi" );
        }

        void
        MaxwellPostprocessor::update_dofs_nedelec()
        {
            mDOFs = mCalculator->nedelec_data_h() ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute( const uint aK )
        {
            return ( this->*mComputeFunction )( aK );
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_air( const uint aK )
        {
            // compute the magnetic field
            mH = mCalculator->B( aK ) * mDOFs ;

            mH *= -1. ;

            return mH ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_ferro( const uint aK )
        {
            // compute the magnetic field
            mH = mCalculator->B( aK ) * mDOFs ;
            mH *= -1.;

            // compute the  permeability and the magnetic flux density
            mB = mMaterial->mu( norm( mH ) ) * mH ;

            combine( mH, mB, mY );

            return mY ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_conductor( const uint aK )
        {
            // compute the magnetic field
            mH = mCalculator->E( aK ) * mDOFs ;

            // compute the  permeability and the magnetic flux density
            mB = mMaterial->mu(  norm( mH ) ) * mH ;

            // compute the current
            mJ = mCalculator->C( aK ) * mDOFs ;

            combine( mH, mB, mJ, mY );

            return mY ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_conductor_ts( const uint aK )
        {
            // compute the magnetic field ( in-plane part, E has no
            // through-thickness row )
            mH = mCalculator->E( aK ) * mDOFs ;

            // compute the  permeability and the magnetic flux density
            const real tMu = mMaterial->mu( norm( mH ) );
            mB = tMu * mH ;

            // compute the current
            mJ = mCalculator->C( aK ) * mDOFs ;

            // recover the normal flux from the volume traces like the
            // superconductor path: per-side dispatch ( phi-region or
            // h-conductor, compute_h_trace — the solve's compute_hn uses
            // the same helper ), average of both sides, projected onto
            // the shell normal

            // dofs from master element
            Vector< real > & phi_m = mCalculator->vector("phi_m");

            // dofs from slave element
            Vector< real > & phi_s = mCalculator->vector("phi_s");

            // normal field component
            Vector< real > & bn = mCalculator->vector("bn");

            // other calculator for facet, and the kind of each volume side
            bool tMasterIsConductor ;
            bool tSlaveIsConductor ;

            Calculator * tCalc = mCalculator->get_normal_calculator(
                    phi_m, phi_s, tMasterIsConductor, tSlaveIsConductor );

            Vector< real > & hk = mCalculator->vector("hk");
            Vector< real > & hm = mCalculator->vector("hm");
            Vector< real > & hs = mCalculator->vector("hs");

            // since we are always having linear elements, we only need to
            // compute the normal field component once
            compute_h_trace( tCalc, true,  tMasterIsConductor, phi_m, hk, hm );
            compute_h_trace( tCalc, false, tSlaveIsConductor,  phi_s, hk, hs );

            bn = hm + hs ;
            bn *= 0.5*constant::mu0 ;

            //Normal vector ;
            const Vector< real > & n = tCalc->normal();
            bn = dot( bn, n ) * n ;

            // the flux keeps the air-side normal component ( [n.b] = 0 ),
            // the field divides it by the permeability
            mB += bn ;
            mH += ( 1.0 / tMu ) * bn ;

            combine( mH, mB, mJ, mY );

            return mY ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_side_connector( const uint aK )
        {
            // full wall field h = ht + hb + hn, recovered the same way as
            // the assembly kernel: MaxwellData links the master layer
            // element and adds the normal and binomial components that the
            // wall's own E * q cannot see
            mH = mCalculator->maxwell()->compute_h( aK );

            // constant mu, enforced at setup by the MaxwellData ctor gate
            mB = mMaterial->mu( norm( mH ) ) * mH ;

            // the current lives in the wall's own edge dofs
            mJ = mCalculator->C( aK ) * mDOFs ;

            combine( mH, mB, mJ, mY );

            return mY ;
        }

        void
        MaxwellPostprocessor::copy_seam_fields()
        {
            // the master holds the full mesh and the gathered fields; the
            // wall nodes are decoupled duplicates whose original() is the
            // live tape rim node, so the copy must run here and not in the
            // rank-local assembly. The source H/B entries are written by the
            // tape postprocessor instances, which run before this one
            // ( factory creation order )
            if ( mCommRank != 0 )
            {
                return ;
            }

            // T only exists on thermal runs; the copy overwrites whatever
            // the patch recovery produced on the wall nodes — J is not in
            // this list on purpose, it stays the wall's own C*q
            const string tCandidates[ 7 ] =
                { "T", "Hx", "Hy", "Hz", "Bx", "By", "Bz" };

            Cell< Vector< real > * > tFields ;
            tFields.reserve( 7 );

            for ( uint f=0; f<7; ++f )
            {
                if ( mMesh->field_exists( tCandidates[ f ] ) )
                {
                    tFields.push( & mMesh->field_data( tCandidates[ f ] ) );
                }
            }

            // station pairing of the lateral faces, see
            // MaxwellData::prepare_side_connector_frame
            const uint tLeftFace[ 4 ]  = { 3, 2, 6, 7 };
            const uint tRightFace[ 4 ] = { 0, 1, 5, 4 };

            for ( mesh::Block * tBlock : mMesh->blocks() )
            {
                const bool tIsLeft =
                    tBlock->domain_type() == DomainType::LeftCoating ;

                if ( ! tIsLeft &&
                     tBlock->domain_type() != DomainType::RightCoating )
                {
                    continue ;
                }

                const uint * tTapeFace  = tIsLeft ? tRightFace : tLeftFace ;
                const uint * tOuterFace = tIsLeft ? tLeftFace  : tRightFace ;

                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    for ( uint k=0; k<4; ++k )
                    {
                        const index_t tSrc = tElement->node(
                            tTapeFace[ k ] )->original()->index();
                        const index_t tIn  = tElement->node(
                            tTapeFace[ k ] )->index();
                        const index_t tOut = tElement->node(
                            tOuterFace[ k ] )->index();

                        // both wall faces show the tape seam value
                        for ( Vector< real > * tField : tFields )
                        {
                            ( *tField )( tIn )  = ( *tField )( tSrc );
                            ( *tField )( tOut ) = ( *tField )( tSrc );
                        }
                    }
                }
            }
        }

        void
        MaxwellPostprocessor::compute_element_data()
        {
            mElementData.fill( 0.0 );

            index_t tCount = 0 ;

            id_t tLastBlockID = gNoID ;

            Cell< mesh::Element * > & tElements = mMesh->elements();

            for ( index_t e : mMyOwnedElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                if ( tElement->owner() == mCommRank )
                {
                    id_t tBlockID = tElement->block_id() ;

                    if ( tBlockID != tLastBlockID )
                    {
                        // get the material
                        mBlock = mField->block( tElement->block_id() );

                        mMaterial = mMaterialMap[ tBlockID ];

                        // get the calculator
                        mCalculator = mBlock->calculator();

                        tLastBlockID = tBlockID ;
                    }

                    mElement = mBlock->element( tElement->id() );

                    // link calculator to element
                    mCalculator->link( mElement );

                    this->update_element_dofs();

                    mZ.fill( 0.0 );
                    const Vector< real > & tW = mCalculator->integration()->weights() ;

                    for ( uint k=0; k<mCalculator->num_intpoints(); ++k )
                    {
                        mZ += tW( k ) * this->compute( k );
                    }

                    mZ /= sum( tW );

                    mElementData.set_col( tCount++, mZ );
                }
            }
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_superconductor( const uint aK )
        {
            // compute the magnetic field
            mH = mCalculator->E( aK ) * mDOFs ;

            // compute the  permeability and the magnetic flux density
            mB = mMaterial->mu(  norm( mH ) ) * mH ;

            // compute the current
            mJ = mCalculator->C( aK ) * mDOFs ;

            // divide the current through the critical current
            mJJc = mJ ;

            Kernel * tThermalKernel = mCalculator->group()->parent()->parent()->controller()->thermal_kernel() ;

            // aura-only blocks may be absent from the thermal dof manager
            bool tIsThermal = tThermalKernel != nullptr
                && tThermalKernel->dofmgr()->block_exists( mCalculator->group()->id() ) ;
            real Temp = gTbulk ;
            if ( tIsThermal )
            {
                Calculator * tCalculator =
                    tThermalKernel->dofmgr()->block( mCalculator->group()->id() )->calculator() ;

                // the thermal peer keeps its previous element when this one is
                // missing from the thermal group ( see
                // link_element_maxwell_thermal ), so its q() would be another
                // element's temperature. Use the bulk value unless the link is
                // current
                if ( tCalculator->element() != nullptr
                     && tCalculator->element()->id() == mCalculator->element()->id() )
                {
                    const Matrix< real > & tN = tCalculator->N( aK ) ;
                    const Vector < real > & T = tCalculator->q() ;

                    Temp = norm(tN * T) ;
                }
            }

            // beta = pi/2 is the REQUIRED convention here, not a placeholder:
            // a bulk HTS element has no tape normal to measure a field angle
            // against, so the solve side passes the same fixed value through
            // calculator::MaxwellData::beta_dummy() ( cl_FEM_Calculator.hpp ),
            // which asserts it is exactly 0.5*pi. Postproc and assembly must
            // agree on the jc argument or the written J/Jc stops matching the
            // resistivity that produced the solution -- that mismatch is
            // precisely the defect the thin-shell branch below was fixed for
            // ( bn_angle, 41b3b280 ). Do not "add theta dependency" here
            // without first giving bulk HTS a defined normal on BOTH sides.
            if ( mMaterial->type() == MaterialType::UserDefined )
            {
                // a user jc is a JcFunction of ( normB, angleNxB ) or
                // ( normB, angleNxB, T ), or a plain constant that both
                // jc() overloads fall back to; the T-only callback path is
                // rejected at registration, so only the arity is open
                if ( mMaterial->depends( MaterialProperty::jc, MaterialDependency::T ) )
                {
                    mJJc /= mMaterial->jc( norm(mB), constant::pi/2, Temp  ) ;
                }
                else
                {
                    mJJc /= mMaterial->jc( norm(mB), constant::pi/2  ) ;
                }
            }
            else
            {
                if (mMaterial->have_defect())
                {
                    const Matrix < real > Coords = Matrix< real >(mCalculator->N(aK)*mCalculator->X()) ;
                    mJJc /= mMaterial->jc( norm(mB), constant::pi/2, Temp,
                                            Coords(0,0), Coords(0,1), Coords(0,2),
                                            mCalculator->group()->parent()->parent()->controller()->time()) ;
                }
                else
                {
                    mJJc /= mMaterial->jc( norm(mB), constant::pi/2, Temp  ) ;
                }
            }

            combine( mH, mB, mJ, mJJc, mY );

            return mY ;
        }

        const Vector< real > &
        MaxwellPostprocessor::compute_superconductor_ts( const uint aK )
        {
            // compute the magnetic field
            mH = mCalculator->E( aK ) * mDOFs ;

            // compute the  permeability and the magnetic flux density
            const real tMu = mMaterial->mu(  norm( mH ) );
            mB = tMu * mH ;

            // compute the current
            mJ = mCalculator->C( aK ) * mDOFs ;

            // divide the current through the critical current
            mJJc = mJ ;

            Kernel * tThermalKernel = mCalculator->group()->parent()->parent()->controller()->thermal_kernel() ;

            // aura-only blocks may be absent from the thermal dof manager
            bool tIsThermal = tThermalKernel != nullptr
                && tThermalKernel->dofmgr()->block_exists( mCalculator->group()->id() ) ;
            real Temp = gTbulk ;
            if ( tIsThermal )
            {
                Calculator * tCalculator =
                    tThermalKernel->dofmgr()->block( mCalculator->group()->id() )->calculator() ;

                // the thermal peer keeps its previous element when this one is
                // missing from the thermal group ( see
                // link_element_maxwell_thermal ), so its q() would be another
                // element's temperature. Use the bulk value unless the link is
                // current
                if ( tCalculator->element() != nullptr
                     && tCalculator->element()->id() == mCalculator->element()->id() )
                {
                    const Matrix< real > & tN = tCalculator->N( aK ) ;
                    const Vector < real > & T = tCalculator->q() ;

                    Temp = norm(tN * T) ;
                }
            }

            // dofs from master element
            Vector< real > & phi_m = mCalculator->vector("phi_m");

            // dofs from slave element
            Vector< real > & phi_s = mCalculator->vector("phi_s");

            // normal field component
            Vector< real > & bn = mCalculator->vector("bn");

            // total field
            Vector< real > & b = mCalculator->vector("b");

            // other calculator for facet, and the kind of each volume side
            bool tMasterIsConductor ;
            bool tSlaveIsConductor ;

            Calculator * tCalc = mCalculator->get_normal_calculator(
                    phi_m, phi_s, tMasterIsConductor, tSlaveIsConductor );

            Vector< real > & hk = mCalculator->vector("hk");
            Vector< real > & hm = mCalculator->vector("hm");
            Vector< real > & hs = mCalculator->vector("hs");

            // since we are always having linear elements, we only need to compute
            // the normal field component once; per-side dispatch phi-region /
            // h-conductor through the same helper as the solve's compute_hn
            compute_h_trace( tCalc, true,  tMasterIsConductor, phi_m, hk, hm );
            compute_h_trace( tCalc, false, tSlaveIsConductor,  phi_s, hk, hs );

            // we use the average of both fields
            bn = hm + hs ;
            bn *= 0.5*constant::mu0 ;

            //Normal vector ;
            const Vector< real > & n = tCalc->normal();

            // keep only the normal component of the air-side average,
            // exactly like compute_hn on the solve side — the unprojected
            // average would add the air-side tangential field to |b|
            bn = dot( bn, n ) * n ;

            // total field
            b = mB + bn ;

            // field magnitude and angle through the assembly's own convention
            // ( bn_angle, unfolded [ 0, pi ] since 2026-08-16 ), so the displayed
            // J/Jc uses the same jc evaluation as the resistivity in the solve
            real norm_b ;
            real beta = mCalculator->bn_angle( b, n, norm_b ) ;

            if ( mMaterial->type() == MaterialType::UserDefined )
            {
                // a user jc is a JcFunction of ( normB, angleNxB ) or
                // ( normB, angleNxB, T ), or a plain constant that both
                // jc() overloads fall back to; the T-only callback path is
                // rejected at registration, so only the arity is open
                if ( mMaterial->depends( MaterialProperty::jc, MaterialDependency::T ) )
                {
                    mJJc /= mMaterial->jc( norm_b, beta, Temp  ) ;
                }
                else
                {
                    mJJc /= mMaterial->jc( norm_b, beta  ) ;
                }
            }
            else
            {
                if (mMaterial->have_defect())
                {
                    const Matrix < real > Coords = Matrix< real >(mCalculator->N(aK)*mCalculator->X()) ;
                    mJJc /= mMaterial->jc( norm_b, beta, Temp,
                                            Coords(0,0), Coords(0,1), Coords(0,2),
                                            mCalculator->group()->parent()->parent()->controller()->time()) ;
                }
                else
                {
                    mJJc /= mMaterial->jc( norm_b, beta, Temp  ) ;
                }
            }

            // ship the recovered normal component with the node fields:
            // the flux keeps the air-side normal component ( [n.b] = 0 ),
            // the field divides it by the permeability
            mB += bn ;
            mH += ( 1.0 / tMu ) * bn ;

            combine( mH, mB, mJ, mJJc, mY );

            return mY ;
        }

        void
        MaxwellPostprocessor::collect_element_data()
        {
            comm_barrier() ;

            if ( mCommRank == 0 )
            {
                Cell< Matrix< real > > tElementData( mCommSize, {} );
                collect( tElementData );

                // loop over all fields
                uint f=0;
                for ( string & tLabel : mElementTargetFields )
                {
                    Vector< real > & tField = mMesh->field_data( tLabel );

                    // proc 0's own elements
                    for ( index_t e=0; e<mMyOwnedElementIndices.size(); ++e )
                    {
                        tField( mMyOwnedElementIndices( e ) ) = mElementData( f, e );
                    }

                    // remote procs
                    for ( proc_t tProc = 1; tProc < mCommSize; ++tProc )
                    {
                        const Matrix< real > & tData = tElementData( tProc );
                        const Vector< index_t > & tIndices = mAllOwnedElementIndices( tProc );
                        index_t tNumElements = tIndices.length();

                        for ( index_t e=0; e<tNumElements; ++e )
                        {
                            tField( tIndices( e ) ) = tData( f, e );
                        }
                    }

                    ++f ;
                }
            }
            else
            {
               send( mElementData );
            }

            comm_barrier() ;
        }

        DomainType
        MaxwellPostprocessor::select_domain_type( const MaxwellPostprocessorType aType )
        {
            switch ( aType )
            {
                case MaxwellPostprocessorType::Air :
                {
                    return DomainType::Air ;
                }
                case MaxwellPostprocessorType::Ferro :
                {
                    return DomainType::Ferro ;
                }
                case MaxwellPostprocessorType::Conductor :
                case MaxwellPostprocessorType::SuperConductor :
                {
                    return DomainType::Conductor ;
                }
                case  MaxwellPostprocessorType::ThinShellConductor :
                case  MaxwellPostprocessorType::ThinShellSuperConductor :
                {
                    return DomainType::ThinShell ;
                }
                case MaxwellPostprocessorType::SideConnector :
                {
                    // the wall instance selects its blocks explicitly ( left
                    // AND right ); this value only feeds the unused
                    // domain-driven fallback
                    return DomainType::LeftCoating ;
                }
                default:
                {
                    BELFEM_ERROR( false, "undefinded domain type" );
                    return DomainType::UNDEFINED ;
                }
            }
        }

    }
}
