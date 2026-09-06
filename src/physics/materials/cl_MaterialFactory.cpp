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

#include "assert.hpp"
#include "fn_check_unit.hpp"
#include "fn_material_data_path.hpp"
#include "cl_MaterialFactory.hpp"

#include "cl_JcFunction.hpp"
#include "cl_Material_Alloy.hpp"

#include "cl_JcFunction_ModifiedKim.hpp"
#include "cl_JcFunction_Database.hpp"

#include "cl_Material_Aluminum.hpp" // 13
#include "cl_Material_Chromium.hpp" // 24
#include "cl_Material_Iron.hpp"     // 26
#include "cl_Material_Nickel.hpp"   // 28
#include "cl_Material_Copper.hpp"   // 29
#include "cl_Material_Silver.hpp"   // 47
#include "cl_Material_Indium.hpp"   // 49
#include "cl_Material_WhiteTin.hpp" // 50
#include "cl_Material_Lead.hpp"     // 82

#include "cl_Material_Magnesia.hpp"

#include "cl_Material_HastelloyC276.hpp"
#include "cl_Material_YBCO.hpp"

#include "cl_Material_BhSplineCurve.hpp"
#include "cl_Material_UserDefined.hpp"

namespace belfem
{
    MaterialFactory::MaterialFactory( const input::Section* aSection )
    {
        uint n = aSection->num_sections();

        for ( uint d=0; d<n; ++d )
        {
            const input::Section * tMatSection = aSection->section( d ) ;
            string tMatLabel = tMatSection->type();

            Material * tMat = nullptr ;

            // Ferromagnetic material: pull a named B-H curve from an HDF5
            // database, e.g.
            //     iron { builtin ; bhfile: bhdata.hdf5 ; curve: RoxieIron ; }
            // The 'curve' key selects the group inside the database file; the
            // optional 'bhfile' key overrides the default database path.
            if ( tMatSection->key_exists( "curve" ) )
            {
                const string tBhFile = tMatSection->key_exists( "bhfile" )
                    ? tMatSection->get_string( "bhfile" ) : "bhdata.hdf5" ;

                tMat = new material::Iron( BELFEM_QUIET_NAN, false ) ;
                tMat->load_bh_curve(
                    this->create_bh_curve( tBhFile, tMatSection->get_string( "curve" ) ) ) ;

                // The curve branch leaves the loop below, so it needs
                // its own check. 'builtin' is tolerated and ONLY here: a
                // valueless 'builtin ;' beside 'curve' is a losing selector
                // that two shipped decks carry and that input_schema.yaml
                // already documents as inert. The exemption is exactly that
                // one key -- 'builtin iron' is a different key name and stays
                // an error, because it does nothing and reads as if it did
                check_unused_input( tMatSection, tMatLabel,
                                    "b-h curve ferromagnet",
                                    { "curve", "bhfile", "builtin" },
                                    {} ) ;

                BELFEM_ERROR( ! mMaterialsMap.key_exists( tMatLabel ),
                    "Material '%s' is defined more than once. The later "
                    "definition used to win silently and leak the earlier one.",
                    tMatLabel.c_str() ) ;

                mMaterialsMap[ tMatLabel ] = tMat ;
                continue ;
            }

            //Access built-in materials
            string tMatType ;
            bool tIsBuiltin = false ;

            if ( tMatSection->key_exists( "builtin" ) )
            {
                tMatType = tMatSection->get_string( "builtin" ) ;
                tIsBuiltin = true ;
            }
            else if ( string_to_lower( tMatSection->label() ) == "builtin" )
            {
                tMatType = tMatLabel ;
                tIsBuiltin = true ;
            }
            else if ( aSection->key_exists( tMatLabel ) )
            {
                if ( string_to_lower( aSection->get_string( tMatLabel ) ) == "builtin" )
                {
                    tMatType = tMatLabel ;
                    tIsBuiltin = true ;
                }
            }
            if ( tIsBuiltin )
            {
                // RRR must be known at construction: an Alloy derives rho_0
                // from it in its constructor and does not override set_RRR()
                // ( the base version aborts ). The pure
                // metals apply a non-NaN RRR in their own constructors, so
                // forwarding it here also builds the rho database only once.
                const real tRRR = tMatSection->key_exists( "RRR" )
                        ? tMatSection->get_real( "RRR" ) : BELFEM_QUIET_NAN ;

                tMat = create_material( tMatType, tRRR ) ;
                switch ( tMat->type() )
                {
                    case MaterialType::HTS :
                    {
                        // HTS materials accept either a database file
                        // (Jc(B,angle,T) / n(B,angle,T) tables) or constant
                        // values for jc, n, and optionally ec.
                        if ( tMatSection->key_exists( "file" ) )
                        {
                            tMat->set_jc_function( this->create_jc_function( tMatSection->get_string( "file" ), "jc" ) ) ;
                            tMat->set_n_function( this->create_jc_function( tMatSection->get_string( "file" ), "n" )  ) ;
                        }
                        else
                        {
                            BELFEM_ERROR(
                                tMatSection->key_exists( "jc" ) && tMatSection->key_exists( "n" ),
                                "HTS material '%s' requires either a 'file' key or constants 'jc' and 'n'",
                                tMatLabel.c_str() ) ;

                            tMat->set_constant( MaterialProperty::jc,
                                                tMatSection->get_value( "jc", "A/m^2" ).first ) ;

                            tMat->set_constant( MaterialProperty::n,
                                                tMatSection->get_real( "n" ) ) ;

                            // ec defaults to 1e-4 V/m to match the file-path
                            // behaviour in Material::set_jc_function
                            tMat->set_constant( MaterialProperty::ec,
                                                tMatSection->key_exists( "ec" )
                                                    ? tMatSection->get_value( "ec", "V/m" ).first
                                                    : 1e-4 ) ;
                        }

                        if (tMatSection->section_exists( "defect" ))
                        {
                            const input::Section * tDefectSection = tMatSection->section( "defect" ) ;
                            BELFEM_ERROR( tDefectSection->key_exists( "file" ), "File to defect must be defined" ) ;
                            BELFEM_ERROR( tDefectSection->key_exists( "label" ), "Label of defect must be defined" ) ;
                            tMat->read_defect( tDefectSection->get_string( "file" ), tDefectSection->get_string( "label" )  ) ;

                        }

                        //Check the resistivity function type (power-law or piecewise)
                        if (tMatSection->key_exists( "resistivity type" ))
                        {
                            const string tResType = tMatSection->get_string( "resistivity type" ) ;
                            if (tResType == "powerlaw" || tResType == "power-law")
                            {
                                tMat->set_piecewise( false ) ;
                            }
                            else if (tResType == "piecewise")
                            {
                                tMat->set_piecewise( true ) ;
                            }
                            else if (tResType == "riva")
                            {
                                tMat->set_resistivity_law( ResistivityLaw::Riva ) ;
                            }
                            else
                            {
                                BELFEM_ERROR( false, "Unknown resistivity type" ) ;
                            }
                        }
                    }
                    default:
                    {
                        break ;
                    }
                }
            }
            else if (tMatSection->section_exists( "usermat" ))
            {
                // The label "buffer" is reserved: when it appears as a layer
                // material name in a thin-shell stack, ThinShellFactory tags
                // the corresponding block as DomainType::Buffer and routes
                // the magnetic kernel through the scalar-phi formulation
                // (no diffusion). Allowing a user-defined or builtin material
                // to also register under this name would silently override
                // that special handling, so we reject it at the source.
                BELFEM_ERROR( string_to_lower( tMatLabel ) != "buffer",
                    "Material name 'buffer' is reserved for the thin-shell "
                    "phi-formulation buffer layer and cannot be used as a "
                    "user-defined or builtin material label." );

                const input::Section * tUserMatSection = tMatSection->section("usermat") ;
                if ( tUserMatSection->key_exists( "file" ) )
                {
                    //Access the user defined material library
                    string tMatLib = tUserMatSection->get_string( "file" ) ;
                    BELFEM_ERROR( tUserMatSection->key_exists( "label" ), "A 'usermat' section requires a 'label' key naming the plugin's init symbol." ) ;
                    tMat = create_material( tUserMatSection->get_string( "file" ), tUserMatSection->get_string( "label" ) ) ;

                    if ( tMat->have(MaterialProperty::jc) ) // check if is userdefined superconductor
                    {
                        if (tMatSection->section_exists( "defect" ))
                        {
                            const input::Section * tDefectSection = tMatSection->section( "defect" ) ;
                            BELFEM_ERROR( tDefectSection->key_exists( "file" ), "File to defect must be defined" ) ;
                            BELFEM_ERROR( tDefectSection->key_exists( "label" ), "Label of defect must be defined" ) ;
                            tMat->read_defect( tDefectSection->get_string( "file" ), tDefectSection->get_string( "label" )  ) ;

                        }

                        //Check the resistivity function type (power-law or piecewise)
                        if (tMatSection->key_exists( "resistivity type" ))
                        {
                            const string tResType = tMatSection->get_string( "resistivity type" ) ;
                            if (tResType == "powerlaw" || tResType == "power-law")
                            {
                                tMat->set_piecewise( false ) ;
                            }
                            else if (tResType == "piecewise")
                            {
                                tMat->set_piecewise( true ) ;
                            }
                            else if (tResType == "riva")
                            {
                                tMat->set_resistivity_law( ResistivityLaw::Riva ) ;
                            }
                            else
                            {
                                BELFEM_ERROR( false, "Unknown resistivity type" ) ;
                            }
                        }

                        // the plugin has already run its _init here. riva
                        // reads n at every superconducting evaluation; with
                        // no n source a release build would run fully
                        // normal via NaN, silently
                        if ( tMat->resistivity_law() == ResistivityLaw::Riva )
                        {
                            BELFEM_ERROR( tMat->have( MaterialProperty::n ),
                                "resistivity type 'riva' on material '%s': the plugin "
                                "registered jc but no n source ( function or constant ).",
                                tMatLabel.c_str() ) ;
                        }
                    }
                }
                else
                {
                    BELFEM_ERROR( false, "A 'usermat' section requires a 'file' key naming the plugin library. ( Defining material properties inline in the input file is not implemented. )" ) ;
                }
            }
            else if (tMatSection->section_exists( "custom" ))
            {
                // The subsection was called 'custom' until 2026-08-29. Without
                // this arm such a deck matches no shape at all and dies on the
                // message below, which never mentions the subsection sitting in
                // the same block. This is a DIAGNOSTIC, not an alias: the old
                // spelling still does not load a material.
                BELFEM_ERROR( false,
                    "Material '%s' has a 'custom' subsection. That subsection was renamed to 'usermat' "
                    "( BELFEM 0.9.0 ) -- rename it in the input file.",
                    tMatLabel.c_str() ) ;
            }
            else
            {
                BELFEM_ERROR( false, "Material '%s' must have either a 'usermat' section, a 'builtin' key, or a 'curve' key", tMatLabel.c_str() ) ;
            }

            // The temperature above which the piecewise and riva resistivity
            // laws skip the superconducting branch and return rho_n. The plain
            // power law does not read it. The correct value belongs to the
            // jc/n data rather than to the carrier material, so a plugin that
            // ships its own fits must be able to state it.
            //
            // Builtins are deliberately NOT overridable. Their T_crit is part
            // of a calibration already consumed at construction: YBCO samples
            // its Callaway lambda spline in the constructor
            // ( cl_Material_YBCO.cpp:122 ) and that sampling copies T_crit into
            // the model parameters ( :507 ). Moving the constant afterwards
            // would shift the resistivity gate while leaving the thermal
            // conductivity built around the old value. Refuse rather than
            // produce a silently inconsistent material.
            if ( tMatSection->key_exists( "critical temperature" ) )
            {
                BELFEM_ERROR( ! tIsBuiltin,
                    "Key 'critical temperature' cannot be used on the builtin material '%s'. "
                    "It applies only to a user-defined material, which owns its own jc and n data.",
                    tMatLabel.c_str() ) ;

                BELFEM_ERROR( tMat->have( MaterialProperty::jc ),
                    "Key 'critical temperature' applies only to a superconductor, "
                    "but material '%s' has no jc.",
                    tMatLabel.c_str() ) ;

                const real tTcrit = tMatSection->get_value( "critical temperature", "K" ).first ;

                BELFEM_ERROR( std::isfinite( tTcrit ) && tTcrit > 0.0,
                    "Critical temperature of material '%s' must be finite and positive ( got %g K ).",
                    tMatLabel.c_str(), ( double ) tTcrit ) ;

                tMat->set_constant( MaterialProperty::T_crit, tTcrit ) ;
            }

            if ( tMatSection->key_exists( "density correction" ) )
            {
                tMat->set_constant( MaterialProperty::density_correction, tMatSection->get_real( "density correction" ) ) ;
            }

            // artificial volumetric heat load, W/m³ as a function of ( x, t ):
            // any builtin or user material, superconductor or not
            if ( tMatSection->section_exists( "heating" ) )
            {
                const input::Section * tHeatSection = tMatSection->section( "heating" ) ;
                BELFEM_ERROR( tHeatSection->key_exists( "file" ), "File to heating plugin must be defined" ) ;
                BELFEM_ERROR( tHeatSection->key_exists( "label" ), "Label of heating plugin must be defined" ) ;
                tMat->read_heating( tHeatSection->get_string( "file" ), tHeatSection->get_string( "label" ) ) ;
            }

            // One allow-list per resolved shape. The sub-shapes are
            // not cosmetic: 'file' and the jc/n/ec constants are an exclusive
            // choice above, so a single list holding all four would accept a
            // deck that sets both and silently drop the constants -- which is
            // the very class this check exists to close
            if ( tIsBuiltin )
            {
                // whether RRR survives is a property of the constructor, not
                // of the name: Alloy forwards it and reports PureMetal, while
                // HastelloyC276, Magnesia and YBCO all drop it ( YBCO then
                // hardcodes 50 ). Gate on the type so a new material inherits
                // the right answer instead of a stale name list
                Cell< string > tKeys = { "builtin", "density correction" } ;
                Cell< string > tSections = { "heating" } ;

                if ( tMat->type() == MaterialType::PureMetal )
                {
                    tKeys.push( "RRR" ) ;
                }

                if ( tMat->type() == MaterialType::HTS )
                {
                    tKeys.push( "resistivity type" ) ;
                    tSections.push( "defect" ) ;

                    if ( tMatSection->key_exists( "file" ) )
                    {
                        tKeys.push( "file" ) ;
                    }
                    else
                    {
                        tKeys.push( "jc" ) ;
                        tKeys.push( "n" ) ;
                        tKeys.push( "ec" ) ;
                    }
                }

                check_unused_input( tMatSection, tMatLabel,
                                    "builtin material", tKeys, tSections ) ;
            }
            else
            {
                // plugin. Everything past the plugin itself is gated on the
                // plugin actually carrying jc, exactly as the branch above is
                Cell< string > tKeys = { "density correction" } ;
                Cell< string > tSections = { "usermat", "heating" } ;

                if ( tMat->have( MaterialProperty::jc ) )
                {
                    tKeys.push( "resistivity type" ) ;
                    tKeys.push( "critical temperature" ) ;
                    tSections.push( "defect" ) ;
                }

                check_unused_input( tMatSection, tMatLabel,
                                    "user-defined material", tKeys, tSections ) ;

                check_unused_input( tMatSection->section( "usermat" ), tMatLabel,
                                    "usermat subsection",
                                    { "file", "label" }, {} ) ;
            }

            if ( tMatSection->section_exists( "defect" ) )
            {
                check_unused_input( tMatSection->section( "defect" ), tMatLabel,
                                    "defect subsection",
                                    { "file", "label" }, {} ) ;
            }

            if ( tMatSection->section_exists( "heating" ) )
            {
                check_unused_input( tMatSection->section( "heating" ), tMatLabel,
                                    "heating subsection",
                                    { "file", "label" }, {} ) ;
            }

            BELFEM_ERROR( ! mMaterialsMap.key_exists( tMatLabel ),
                "Material '%s' is defined more than once. The later definition "
                "used to win silently and leak the earlier one.",
                tMatLabel.c_str() ) ;

            mMaterialsMap[tMatLabel] = tMat ;
        }
    }

    void
    MaterialFactory::check_unused_input(
            const input::Section * aSection,
            const string         & aLabel,
            const string         & aShape,
            const Cell< string > & aKeys,
            const Cell< string > & aSections )
    {
        // key() hands back the name as STORED, which the parser lowercased
        // ( cl_Input_Section.cpp:111 ). Callers of key_exists() never notice,
        // because that lowercases its argument for them -- but a direct
        // compare does, so normalize both sides. Spelling an allow-list entry
        // 'RRR' and comparing it raw would reject every copper deck in the
        // tree
        index_t tNumKeys = aSection->num_keys() ;

        for ( index_t k=0; k<tNumKeys; ++k )
        {
            const string tKey = string_to_lower( aSection->key( k ) ) ;

            bool tKnown = false ;

            for ( const string & tAllowed : aKeys )
            {
                if ( tKey == string_to_lower( tAllowed ) )
                {
                    tKnown = true ;
                    break ;
                }
            }

            BELFEM_ERROR( tKnown,
                "Material '%s' is a %s, which never reads the key '%s'. "
                "Remove the key, or define the material differently.",
                aLabel.c_str(), aShape.c_str(), tKey.c_str() ) ;
        }

        index_t tNumSections = aSection->num_sections() ;

        for ( index_t d=0; d<tNumSections; ++d )
        {
            const input::Section * tSub = aSection->section( d ) ;

            // a labelled subsection is never reachable: Section only maps the
            // UNLABELLED ones by type ( cl_Input_Section.cpp:82-84 ), and every
            // consumer here looks them up by type alone. Refusing it is the
            // whole point -- 'defect : mine { }' silently did nothing
            BELFEM_ERROR( tSub->label().empty(),
                "Material '%s' carries a labelled subsection '%s : %s'. "
                "Subsections here are looked up by type alone, so a labelled "
                "one is never read. Drop the label.",
                aLabel.c_str(), tSub->type().c_str(), tSub->label().c_str() ) ;

            bool tKnown = false ;

            for ( const string & tAllowed : aSections )
            {
                if ( string_to_lower( tSub->type() ) == string_to_lower( tAllowed ) )
                {
                    tKnown = true ;
                    break ;
                }
            }

            BELFEM_ERROR( tKnown,
                "Material '%s' is a %s, which never reads the subsection '%s'. "
                "Remove the subsection, or define the material differently.",
                aLabel.c_str(), aShape.c_str(), tSub->type().c_str() ) ;
        }
    }

//------------------------------------------------------------------------------

    Material *
    MaterialFactory::create_material( const string & aLabel,  const real aRRR, const bool aBuildTables )
    {
        string tLabel = string_to_lower( aLabel );

        if ( tLabel == "aluminum" ||tLabel == "aluminium"  || tLabel == "al" ) return new material::Aluminum( aRRR, aBuildTables ) ;
        if ( tLabel == "chromium"  || tLabel == "cr" ) return new material::Chromium( aRRR, aBuildTables ) ;

        if ( tLabel == "iron" || tLabel == "ferro" || tLabel == "fe" )  return new material::Iron( aRRR, aBuildTables ) ;

        if ( tLabel == "nickel" || tLabel == "ni")   return new material::Nickel( aRRR, aBuildTables );

        if ( tLabel == "copper" || tLabel == "cu")   return new material::Copper( aRRR, aBuildTables );

        if ( tLabel == "silver" || tLabel == "ag") return new material::Silver( aRRR, aBuildTables ) ;

        if ( tLabel == "indium" || tLabel == "in") return new material::Indium( aRRR, aBuildTables ) ;

        if ( tLabel == "tin" || tLabel == "sn" ) return new material::WhiteTin( aRRR, aBuildTables ) ;

        if ( tLabel == "lead" || tLabel == "pb") return new material::Lead( aRRR, aBuildTables ) ;


        if ( tLabel == "hastelloy" || tLabel == "hastelloyc276"  ) return new material::HastelloyC276 ;
        if ( tLabel == "ybco" ) return new material::YBCO ;

        if ( tLabel == "magnesia" || tLabel == "mgo" || tLabel == "buffer" ) return new material::Magnesia ;



        Cell< std::pair< string, real > > tComponents ;
        to_pair( aLabel, tComponents );

        if ( tComponents.size() > 0 )
        {
            return new material::Alloy( aLabel, tComponents, aRRR, aBuildTables ) ;
        }

        BELFEM_ERROR( false, "Material %s is not known", aLabel.c_str() );
        return nullptr ;
    }

    void
    MaterialFactory::print_material_list( std::ostream & aStream )
    {
        Cell< std::pair< uint, string > > tElements ;
        tElements.push( std::pair< uint, string >( 13, "aluminum" ) ) ;
        tElements.push( std::pair< uint, string >( 24, "chromium" ) ) ;
        tElements.push( std::pair< uint, string >( 26, "iron" ) ) ;
        tElements.push( std::pair< uint, string >( 28, "nickel" ) ) ;
        tElements.push( std::pair< uint, string >( 29, "copper" ) ) ;
        tElements.push( std::pair< uint, string >( 47, "silver" ) ) ;
        tElements.push( std::pair< uint, string >( 49, "indium" ) ) ;
        tElements.push( std::pair< uint, string >( 50, "tin" ) ) ;
        tElements.push( std::pair< uint, string >( 82, "lead" ) ) ;

        Cell< string > tAlloys ;
        tAlloys.push( "hastelloyc276" ) ;

        Cell< string > tOther ;
        tOther.push( "magnesia" ) ;
        tOther.push( "ybco" ) ;
        aStream << std::endl ;
        aStream << "    Elements : " << std::endl ;
        aStream << "    ----------" << std::endl ;
        for ( auto tPair : tElements )
        {
            aStream << "      * " << tPair.second << std::endl ;
        }

        aStream << std::endl ;
        aStream << "    Alloys : " << std::endl ;
        aStream << "    -------" << std::endl ;

        for ( string & tAlloy : tAlloys )
        {
            aStream << "      * " << tAlloy << std::endl ;
        }
        aStream << std::endl ;

        aStream << "    Other : " << std::endl ;
        aStream << "    -------" << std::endl ;
        for ( string & tMat : tOther )
        {
            aStream << "      * " << tMat << std::endl ;
        }
        aStream << std::endl ;

        aStream << "    Solder (estimates) : " << std::endl ;
        aStream << "    --------------------" << std::endl ;
        aStream << "      * <Element><percentage>, e.g. Sn40Pb60 or Sn54Pb26In20" << std::endl ;
        aStream << std::endl ;
    }

    Material *
    MaterialFactory::create_material( const string & aLibraryPath, const string & aLabel )
    {
        return new material::UserDefinedMaterial(
            material::data_file( aLibraryPath ), aLabel ) ;
    }

    material::BhCurve *
    MaterialFactory::create_bh_curve( const string & aPath, const string & aLabel )
    {
        // BhCurve is an abstract interface; the spline-backed BhSplineCurve is
        // the concrete HDF5-loading implementation.
        return new material::BhSplineCurve(
            material::data_file( aPath ), aLabel ) ;
    }

    material::JcFunction *
    MaterialFactory::create_jc_function( const real aJc0, const real aB, const real aBc, const real aK )
    {
        return new material::JcFunctionModifiedKim( aJc0, aB, aBc, aK );
    }

    material::JcFunction *
    MaterialFactory::create_jc_function( const string & aPath, const string & aLabel )
    {
        return new material::JcFunctionDatabase(
            material::data_file( aPath ), aLabel );
    }

}
