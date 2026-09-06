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
#include "cl_FEM_Domain.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "en_SolverEnums.hpp"
#include "en_DomainType.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace fem
    {
        Domain::Domain( const input::Section * aSection )
        {
            mType = domain_type( aSection->type() );

            mLabel = aSection->label() == "" ?
                to_string( mType  ): aSection->label();

            switch( mType )
            {
                case DomainType::Air :
                case DomainType::Coil :
                {
                    this->read_groups( aSection, true );
                    break ;
                }
                case DomainType::Ferro :
                case DomainType::Conductor :
                case DomainType::Buffer :
                case DomainType::Default :
                {
                    this->read_groups( aSection, true );
                    BELFEM_ERROR( aSection->key_exists( "material" ),
                        "no material assigned to block %s", mLabel.c_str() );

                    mMaterialLabel = aSection->get_string( "material" );
                    break ;
                }
                case DomainType::AirSymmetry :
                case DomainType::BufferSymmetry :
                case DomainType::FerroSymmetry :
                case DomainType::ConductorSymmetry :
                case DomainType::AirAntiSymmetry :
                case DomainType::BufferAntiSymmetry :
                case DomainType::FerroAntiSymmetry :
                case DomainType::ConductorAntiSymmetry :
                case DomainType::BackgroundField :
                {
                    this->read_groups( aSection, false );
                    mIsBlock = false ;
                    mIsSideSet = true ;
                    break ;
                }
                case DomainType::ThinShell :
                {
                    // thin-shell sidesets may be signed ( orientation flip,
                    // see MaxwellFactory::read_signed_sidesets ); the domain
                    // only needs the absolute ids
                    this->read_groups( aSection, false, true );
                    mIsBlock = false ;
                    mIsSideSet = true ;
                    break ;
                }
                case DomainType::Curve :
                {
                    mIsBlock = false ;
                    break ;
                }
                default:
                {
                    break ;
                }
            }

        }

        void
        Domain::read_groups( const input::Section * aSection,
                             const bool aUseBlocks,
                             const bool aAllowSigns )
        {
            string tGroup  = aUseBlocks ? "block" : "sideset";
            string tGroups = aUseBlocks ? "blocks" : "sidesets";
            // let's be a little tolerant here with the names
            // input errors, not logic bugs: these must fire in release too
            BELFEM_ERROR( aSection->key_exists( tGroup ) || aSection->key_exists( tGroups ),
                 "Error in input file: group %s : no %s defined",
                 mLabel.c_str(), tGroups.c_str() );

            BELFEM_ERROR( aSection->key_exists( tGroup ) xor aSection->key_exists( tGroups ),
                             "Error in input file: group %s : need either %s or %s, not both",
                             mLabel.c_str(), tGroup.c_str(), tGroups.c_str() );

            const string & tKey =
                aSection->key_exists( tGroup ) ? tGroup : tGroups ;

            if ( aAllowSigns )
            {
                // strip gmsh-style orientation signs before the unsigned
                // parse; the signs themselves are consumed by
                // MaxwellFactory::read_signed_sidesets, this reader only
                // needs the absolute ids
                aSection->ids_from_string(
                    search_and_replace( aSection->get_string( tKey ), "-", "" ),
                    mGroupIDs );
            }
            else
            {
                aSection->get_ids( tKey, mGroupIDs );
            }
        }

    }

}
