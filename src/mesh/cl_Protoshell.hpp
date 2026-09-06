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

#ifndef CL_PROTOSHELL_HPP
#define CL_PROTOSHELL_HPP
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Curve.hpp"

namespace belfem
{
    /**
     * @brief Thin shell configuration object for electromagnetic simulations
     * 
     * Protoshells act as templates/blueprints for creating thin shell structures
     * in Maxwell simulations. They store all metadata needed to create the actual
     * thin shell elements during mesh processing.
     * 
     * A Protoshell represents a thin conducting layer (like coil insulation or 
     * shielding) that is geometrically thin compared to the overall problem size
     * but electrically significant.
     * 
     * The class separates configuration (stored here) from mesh modification 
     * (handled by CutFactory), making the code more modular.
     * 
     * Usage workflow:
     * 1. Created by MaxwellFactory::read_thin_shell_data() from input file
     * 2. Passed to CutFactory for preprocessing terminals and creating cuts
     * 3. Results in "tapes" - the processed thin shell surfaces in the mesh
     */
    class Protoshell
    {
        const id_t     mID ;           ///< Unique identifier for this protoshell
        string         mLabel ;        ///< Human-readable name (e.g., "coil_insulation")
        
        // Geometry definition
        Vector< id_t > mSideSetIDs ;   ///< IDs of sidesets that form this thin shell surface
        
        // Electrical connections
        Vector< id_t > mTerminals ;    ///< IDs of terminal points where current can enter/exit
        
        // Material properties (supports multi-layered shells)
        Vector< real > mThicknesses ;  ///< Layer thicknesses in meters
        Cell< string > mMaterials ;    ///< Material names for each layer

        // Curve definitions (note: curves are deleted by mesh class)
        Cell< mesh::Curve * > mSideCurves ;     ///< Curves defining shell boundaries
        Cell< mesh::Curve * > mTerminalCurves ; ///< Curves defining terminal locations

        // Edge coating (surround plating on the tape slit edges, 3D only)
        bool mEdgeCoating = false ;      ///< opt-in: create side connector walls
        real mEdgeCoatingWidth = 0.0 ;   ///< wall width in m; 0 = derive from outer layer thickness

        // Sidesets listed with a NEGATIVE sign in the deck, gmsh style:
        //     sidesets : -5, -6, 7:20 ;
        // A sign requests that the orientation of every facet of that sideset
        // be flipped ( master/slave swap, winding rewrite ) after the master
        // normalization has run. This exists because the layer stack is laid
        // along the facet normal and the domain-type master rule is
        // mirror-symmetric: a stack bounded by air on both faces cannot come
        // out uniform — one outer tape always flips. Which side the layers
        // face is user intent that no geometry rule can derive ( a corc wrap
        // has no meaningful mean normal at all ), so the user declares it per
        // sideset. Stored as positive ids; mSideSetIDs carries the absolute
        // values, so every other consumer is sign-agnostic.
        Vector< id_t > mFlippedSideSets ; ///< sidesets to flip, empty when none
    public:

        Protoshell( const id_t aID ) :
            mID( aID )
        {

        }

        ~Protoshell() = default ;

        string &
        label() ;

        Vector< id_t > &
        sidesets() ;

        Vector< id_t > &
        terminals() ;

        Vector< real > &
        thicknesses() ;

        Cell< string > &
        materials() ;

        Cell< mesh::Curve * > &
        side_curves() ;

        Cell< mesh::Curve * > &
        terminal_curves() ;

        bool &
        edge_coating() ;

        real &
        edge_coating_width() ;

        Vector< id_t > &
        flipped_sidesets() ;

        id_t
        id() const ;

    };

    inline string &
    Protoshell::label()
    {
        return mLabel ;
    }

    inline Vector< id_t > &
    Protoshell::sidesets()
    {
        return mSideSetIDs ;
    }

    inline Vector< id_t > &
    Protoshell::terminals()
    {
        return mTerminals ;
    }

    inline Vector< real > &
    Protoshell::thicknesses()
    {
        return mThicknesses ;
    }

    inline Cell< string > &
    Protoshell::materials()
    {
        return mMaterials;
    }

    inline Cell< mesh::Curve * > &
    Protoshell::side_curves()
    {
        return mSideCurves;
    }

    inline Cell< mesh::Curve * > &
    Protoshell::terminal_curves()
    {
        return mTerminalCurves;
    }

    inline bool &
    Protoshell::edge_coating()
    {
        return mEdgeCoating;
    }

    inline real &
    Protoshell::edge_coating_width()
    {
        return mEdgeCoatingWidth;
    }

    inline Vector< id_t > &
    Protoshell::flipped_sidesets()
    {
        return mFlippedSideSets ;
    }

    inline id_t
    Protoshell::id() const
    {
        return mID;
    }

}
#endif //CL_PROTOSHELL_HPP
