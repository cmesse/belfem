//
// Created by christian on 10/4/21.
//

#ifndef BELFEM_EN_FEM_DOMAINTYPE_HPP
#define BELFEM_EN_FEM_DOMAINTYPE_HPP

namespace belfem
{
    namespace fem
    {
        //
        enum class DomainType
        {
            Default = 0,         // default setting for blocks and sidesets
            Air             =  1,  // Maxwell Specific, do not change order!
            Coil            =  2,  // Maxwell Specific, do not change order!
            Ferro           =  3,  // Maxwell Specific, do not change order!
            Conductor       =  4,  // Maxwell Specific, do not change order!
            InterfaceScAir  =  5,  // Maxwell Specific
            InterfaceFmAir  =  6,  // Maxwell Specific
            InterfaceScFm   =  7,  // Maxwell Specific
            InterfaceFmFm   =  8,  // Maxwell Specific
            Symmetry        =  9,  // Maxwell Specific
            AntiSymmetry    =  10,  // Maxwell Specific
            SymmetryConductor      = 11,  // Maxwell Specific
            SymmetryFerro      = 12,  // Maxwell Specific
            SymmetryAir     = 13,  // Maxwell Specific
            AntiSymmetryConductor  = 14,  // Maxwell Specific
            AntiSymmetryFerro  = 15,  // Maxwell Specific
            AntiSymmetryAir = 16,  // Maxwell Specific
            Cut             = 17,  // Maxwell Specific
            Boundary        = 18,  // Maxwell Specific
            ThinShell       = 19,  // Maxwell Specific
            Ghost           = 20,  // Maxwell Specific
            UNDEFINED = 21
        };

//------------------------------------------------------------------------------

        string
        to_string( const DomainType aDomaintype );

//------------------------------------------------------------------------------

        DomainType
        domain_type( const string & aString );

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_EN_FEM_DOMAINTYPE_HPP
