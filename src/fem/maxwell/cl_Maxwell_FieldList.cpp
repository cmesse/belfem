//
// Created by christian on 12/20/21.
//

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
                    for( const string & tDof : InterfaceScAir )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : InterfaceScFm )
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
                    for( const string & tDof : SymmetrySc )
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
                    for( const string & tDof : AntiSymmetrySc )
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
                    for( const string & tDof : BoundarySc )
                    {
                        Dofs.push( tDof );
                    }
                    for( const string & tDof : Farfield )
                    {
                        Dofs.push( tDof );

                        string tLastDof = tDof + "0";
                        NonDof.push(tLastDof );
                        Hidden.push( tLastDof );
                    }

                    // shell dofs
                    for( const string & tDof : ThinShell )
                    {
                        Dofs.push( tDof );
                        string tLastDof = tDof + "0";

                        ThinShellLast.push(tLastDof );
                        NonDof.push(tLastDof );
                        Hidden.push( tLastDof );
                    }

                    // superconductor dofs
                    for( const string & tDof : Superconductor )
                    {
                        Dofs.push( tDof );
                        InterfaceScFm.push( tDof );
                        InterfaceScAir.push( tDof );
                        SymmetrySc.push( tDof );
                        AntiSymmetrySc.push( tDof );
                        BoundarySc.push( tDof );
                        ThinShell.push( tDof );
                    }

                    // ferro dofs
                    for( const string & tDof : Ferro )
                    {
                        Dofs.push( tDof );
                        InterfaceScFm.push( tDof );
                        InterfaceFmAir.push( tDof );
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
                        InterfaceScAir.push( tDof );
                        InterfaceFmAir.push( tDof );
                        ThinShell.push( tDof );

                        SymmetryAir.push( tDof );
                        AntiSymmetryAir.push( tDof );
                        BoundaryAir.push( tDof );
                        Farfield.push( tDof );
                    }

                    // ghost
                    for( string tDof : Ghost )
                    {
                        NonDof.push( tDof );
                    }

                    // check if cut dofs exist
                    if( Cut.size() > 0 )
                    {
                        // cuts
                        for( const string & tDof : Cut )
                        {
                            Dofs.push( tDof );
                        }

                        // add air dofs to cut
                        for( const string & tDof : Air )
                        {
                            Cut.push( tDof );
                        }
                    }


                    aIWG->unique_and_rearrange( Dofs, true );

                    aIWG->unique_and_rearrange( ThinShell );
                    aIWG->unique_and_rearrange( InterfaceScFm );
                    aIWG->unique_and_rearrange( InterfaceScAir );
                    aIWG->unique_and_rearrange( InterfaceFmAir );

                    aIWG->unique_and_rearrange( SymmetrySc );
                    aIWG->unique_and_rearrange( SymmetryFerro );
                    aIWG->unique_and_rearrange( SymmetryAir );

                    aIWG->unique_and_rearrange( AntiSymmetrySc );
                    aIWG->unique_and_rearrange( AntiSymmetryFerro );
                    aIWG->unique_and_rearrange( AntiSymmetryAir );

                    aIWG->unique_and_rearrange( BoundarySc );
                    aIWG->unique_and_rearrange( BoundaryFerro );
                    aIWG->unique_and_rearrange( BoundaryAir );

                    aIWG->unique_and_rearrange( Farfield );

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

                        if ( tDof == "ax" )
                        {
                            NonDof.push("ax0");
                            Hidden.push( "ax0");
                            FerroLast.push( "ax0" );
                        }
                        else if ( tDof == "ay" )
                        {
                            NonDof.push("ay0");
                            Hidden.push( "ay0");
                            FerroLast.push( "ay0" );
                        }
                        else if ( tDof == "az" )
                        {
                            NonDof.push("az0");
                            Hidden.push( "az0");
                            FerroLast.push( "az0" );
                        }
                        else if ( tDof == "phi" )
                        {
                            NonDof.push("phi0");
                            Hidden.push( "phi0");
                        }
                        else if ( tDof == "edge_h" )
                        {
                            NonDof.push("edge_h0");
                        }
                        else if ( tDof == "face_h" )
                        {
                            NonDof.push("face_h0");
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
                        case( DomainType::SuperConductor ) :
                        {
                            // set dofs to superconductor
                            this->create_doftable(
                                    Superconductor,
                                    aBlockDofs( k ) );

                            break ;
                        }
                        case( DomainType::Coil ) :
                        {
                            this->create_doftable(
                                    Coil,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        case( DomainType::FerroMagnetic ) :
                        {
                            this->create_doftable(
                                    Ferro,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        case( DomainType::Air ) :
                        {
                            this->create_doftable(
                                    Air,
                                    aBlockDofs( k ) );
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid block type" );
                        }
                    }
                }
            }

//-----------------------------------------------------------------------------

            void
            FieldList::collect_sideset_dofs(
                    const Vector< id_t >              & aSideSetIDs,
                    const Map< id_t, DomainType >     & aSideSetTypeMap,
                    const Map< id_t, MagfieldBcType > & aSideSetSubTypeMap,
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
                        case( DomainType::InterfaceScFm ) :
                        {
                            this->create_doftable(
                                    InterfaceScFm,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::InterfaceScAir ) :
                        {

                            this->create_doftable(
                                    InterfaceScAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::InterfaceFmAir ) :
                        {
                            this->create_doftable(
                                    InterfaceFmAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::SymmetrySc ) :
                        {
                            this->create_doftable(
                                    SymmetryAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::SymmetryFm ) :
                        {
                            this->create_doftable(
                                    SymmetryFerro,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::AntiSymmetrySc ) :
                        {
                            this->create_doftable(
                                    AntiSymmetrySc,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::AntiSymmetryFm ) :
                        {
                            this->create_doftable(
                                    AntiSymmetryFerro,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::AntiSymmetryAir ) :
                        {
                            this->create_doftable(
                                    AntiSymmetryAir,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::Cut ) :
                        {
                            this->create_doftable(
                                    Cut,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::ThinShell ) :
                        {
                            this->create_doftable(
                                    ThinShell,
                                    aSideSetDofs( k ) );
                            break ;
                        }
                        case( DomainType::Boundary ) :
                        {
                            // todo: add also symmetry BCs here

                            // check domain type
                            switch( aSideSetSubTypeMap( aSideSetIDs( k ) ) )
                            {
                                case( MagfieldBcType::Wave ) :
                                {
                                    this->create_doftable(
                                            BoundaryAir,
                                            aSideSetDofs( k ) );
                                    break ;
                                }
                                case( MagfieldBcType::Farfied ) :
                                {
                                    this->create_doftable(
                                            Farfield,
                                            aSideSetDofs( k ) );
                                    break ;
                                }
                                default :
                                {
                                    BELFEM_ERROR( false, "This Magfield BC type is not implemented yet" );
                                }
                            }

                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid sideset type" );
                        }
                    }

                }

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