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

#include "stringtools.hpp"
#include "fn_entity_type.hpp"
#include "assert.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_FEM_DofMgr_SideSetData.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
 //------------------------------------------------------------------------------

            FieldList::FieldList(
                       Map < string, uint > & aDofMap,
                       Cell< string > & aDofs,
                       Cell< string > & aNonDof,
                       Cell< string > & aHidden,
                       Cell< string > & aAll ) :
                    mDofMap( aDofMap ),
                    Dofs( aDofs ),
                    NonDof( aNonDof ),
                    Hidden( aHidden ),
                    All( aAll )
            {

            }

//------------------------------------------------------------------------------

            void
            FieldList::initialize( IWG_Maxwell * aIWG )
            {
                // we only do something if the dofs have not been set by the IWG
                if( Dofs.size() == 0 )
                {
                    // interfaces
                    for( const string & tDof : InterfaceCondAir )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : InterfaceCondFm )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : InterfaceFmFm )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : InterfaceFmAir )
                    {
                        Dofs.push( tDof );
                    }

                    // antisymmetry
                    for( const string & tDof : SymmetryAir )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : SymmetryFerro )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : SymmetryConductor )
                    {
                        Dofs.push( tDof );
                    }

                    // antisymmetry
                    for( const string & tDof : AntiSymmetryAir )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : AntiSymmetryFerro )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : AntiSymmetryConductor )
                    {
                        Dofs.push( tDof );
                    }

                    // boundary
                    for( const string & tDof : BoundaryAir )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : BoundaryFerro )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : BoundaryConductor )
                    {
                        Dofs.push( tDof );
                    }

                    // shell dofs
                    for( const string & tDof : ThinShell )
                    {
                        Dofs.push( tDof );
                    }

                    // superconductor dofs
                    for( const string & tDof : Conductor )
                    {
                        Dofs.push( tDof );
                        InterfaceCondFm.push( tDof );
                        InterfaceCondAir.push( tDof );
                        SymmetryConductor.push( tDof );
                        AntiSymmetryConductor.push( tDof );
                        BoundaryConductor.push( tDof );
                        ThinShell.push( tDof );
                        Ghost.push( tDof );
                    }

                    // ferro dofs
                    for( const string & tDof : Ferro )
                    {
                        Dofs.push( tDof );
                        InterfaceCondFm.push( tDof );
                        InterfaceFmAir.push( tDof );
                        InterfaceFmFm.push( tDof );
                        SymmetryFerro.push( tDof );
                        AntiSymmetryFerro.push( tDof );
                        BoundaryFerro.push( tDof );
                    }

                    // coils
                    for( const string & tDof : Coil )
                    {
                        Dofs.push( tDof );
                    }

                    // air
                    for( const string & tDof : Air )
                    {
                        Dofs.push( tDof );
                        InterfaceCondAir.push( tDof );
                        InterfaceFmAir.push( tDof );
                        ThinShell.push( tDof );

                        SymmetryAir.push( tDof );
                        AntiSymmetryAir.push( tDof );
                        BoundaryAir.push( tDof );
                    }

                    aIWG->unique_and_rearrange( Dofs, true );

                    aIWG->unique_and_rearrange( ThinShell );

                    aIWG->unique_and_rearrange( InterfaceCondFm );
                    aIWG->unique_and_rearrange( InterfaceFmFm );
                    aIWG->unique_and_rearrange( InterfaceCondAir );
                    aIWG->unique_and_rearrange( InterfaceFmAir );

                    aIWG->unique_and_rearrange( SymmetryConductor );
                    aIWG->unique_and_rearrange( SymmetryFerro );
                    aIWG->unique_and_rearrange( SymmetryAir );
                    aIWG->unique_and_rearrange( AntiSymmetryConductor );
                    aIWG->unique_and_rearrange( AntiSymmetryFerro );
                    aIWG->unique_and_rearrange( AntiSymmetryAir );

                    aIWG->unique_and_rearrange( BoundaryConductor );
                    aIWG->unique_and_rearrange( BoundaryFerro );
                    aIWG->unique_and_rearrange( BoundaryAir );

                    aIWG->unique_and_rearrange( Ghost );

                    // create other dof lists
                    for( string tDof : Dofs )
                    {
                        if( tDof.size() >= 6 )
                        {
                            if( string_to_lower( tDof.substr( 0, 6 ) ) == "lambda" )
                            {
                                Lambda.push( tDof );
                            }
                        }

                        if ( tDof.c_str()[ 0 ] == '_' )
                        {
                            Hidden.push( tDof );
                        }
                    }

                    for(  string tDof : MagneticFieldDensity )
                    {
                        NonDof.push( tDof );
                    }
                    for( string tDof : CurrentDensity )
                    {
                        NonDof.push( tDof );
                    }
                    for( string tDof : CurrentBC )
                    {
                        NonDof.push( tDof );
                        Hidden.push( tDof );
                    }

                    aIWG->unique_and_rearrange( NonDof );
                    aIWG->unique_and_rearrange( Hidden );
                }
                else
                {
                    aIWG->unique_and_rearrange( Dofs, true );
                }

                // assemble all
                for( string tDof : Dofs )
                {
                    All.push( tDof );
                }
                for( string tDof : NonDof )
                {
                    All.push( tDof );
                }

                aIWG->unique_and_rearrange( All );

            }

//-----------------------------------------------------------------------------

            void
            FieldList::collect_block_dofs(
                    const Vector< id_t >             & aBlockIDs,
                    const Map< id_t, DomainType >    & aBlockTypeMap,
                          Cell< Vector< index_t > >  & aBlockDofs )
            {
                // determine the number of blocks that are used
                uint tNumBlocks = aBlockIDs.length();

                // allocate the memory
                aBlockDofs.set_size( tNumBlocks, Vector< id_t> () );

                // loop over all blocks
                for( uint k=0; k<tNumBlocks; ++k )
                {

                    // check the type of the block
                    switch( aBlockTypeMap( aBlockIDs( k ) ) )
                    {
                        // set the block specific dof types
                        case DomainType::Conductor :
                        case DomainType::ThinShell :
                        case DomainType::LeftCoating :
                        case DomainType::RightCoating :
                        {
                            this->create_doftable(
                                    Conductor,
                                    aBlockDofs( k ) );

                            break ;
                        }
                        case DomainType::Coil :
                        {
                            this->create_doftable(
                                    Coil,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        case DomainType::Ferro :
                        {
                            this->create_doftable(
                                    Ferro,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        case DomainType::Air :
                        case DomainType::Buffer :
                        {
                            this->create_doftable(
                                    Air,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid block type for block %lu : %s",
                                ( long unsigned int ) aBlockIDs( k ), to_string(  aBlockTypeMap( aBlockIDs( k ) ) ).c_str() );
                        }
                    }
                }
            }

//-----------------------------------------------------------------------------

            void
            FieldList::collect_sideset_dofs(
                    const Vector< id_t >              & aSideSetIDs,
                    const Map< id_t, DomainType >     & aSideSetTypeMap,
                    Cell< Vector< index_t > >         & aSideSetDofs )
            {
                // determine the number of blocks that are used
                uint tNumSideSets = aSideSetIDs.length();

                // allocate the memory
                aSideSetDofs.set_size( tNumSideSets, Vector< id_t> () );

                for( uint k=0; k<tNumSideSets; ++k )
                {
                    // check the type of the sideset
                    switch( aSideSetTypeMap( aSideSetIDs( k ) ) )
                    {
                        case DomainType::InterfaceCondFerro :
                        {
                            this->create_doftable(
                                    InterfaceCondFm,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::InterfaceCondAir :
                        {

                            this->create_doftable(
                                    InterfaceCondAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::InterfaceFerroAir :
                        {
                            this->create_doftable(
                                    InterfaceFmAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::ConductorSymmetry :
                        {
                            this->create_doftable(
                                    SymmetryConductor,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::AirSymmetry :
                        case DomainType::BufferSymmetry :
                        {
                            this->create_doftable(
                                    SymmetryAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::FerroSymmetry :
                        {
                            this->create_doftable(
                                    SymmetryFerro,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::ConductorAntiSymmetry :
                        {
                            this->create_doftable(
                                    AntiSymmetryConductor,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::FerroAntiSymmetry :
                        {
                            this->create_doftable(
                                    AntiSymmetryFerro,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::AirAntiSymmetry :
                        case DomainType::BufferAntiSymmetry :
                        {
                            this->create_doftable(
                                    AntiSymmetryAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::Cut :
                        {
                            this->create_doftable(
                                    Cut,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::ThinShell :
                        {
                            this->create_doftable(
                                    ThinShell,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::BackgroundField :
                        {
                            // todo: add wave BCs for background field here

                            this->create_doftable(
                                    BoundaryAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::EnrichedInterface :
                        {
                            this->create_doftable(
                                    Air,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::Ghost :
                        {
                            this->create_doftable(
                                    Ghost,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case DomainType::Inactive :
                        case DomainType::GeometryOnly :
                        {
                            // pass
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid sideset type for sideset %lu",
                                          ( long unsigned int ) aSideSetIDs( k ) );
                        }
                    }

                }
            }

//-----------------------------------------------------------------------------

            uint
            FieldList::doftype( const string & aLabel ) const
            {
                return mDofMap( aLabel );
            }

//-----------------------------------------------------------------------------

            void
            FieldList::create_doftable( const Cell< string >    & aDofList,
                                              Vector< index_t > & aDofTable )
            {
                // get size of list
                uint tSize = aDofList.size();

                // allocate memory
                aDofTable.set_size( tSize );

                // create table
                for( uint k=0; k<tSize; ++k )
                {
                    aDofTable( k ) = mDofMap( aDofList( k ) );
                }
            }

//-----------------------------------------------------------------------------
        } /* end namespace maxwell */
    } /* end namespace fem */
} /* end namespace belfem */
