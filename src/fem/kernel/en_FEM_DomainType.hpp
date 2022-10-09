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
            Symmetry        =  8,  // Maxwell Specific
            AntiSymmetry    =  9,  // Maxwell Specific
            SymmetryConductor      = 10,  // Maxwell Specific
            SymmetryFerro      = 11,  // Maxwell Specific
            SymmetryAir     = 12,  // Maxwell Specific
            AntiSymmetryConductor  = 13,  // Maxwell Specific
            AntiSymmetryFerro  = 14,  // Maxwell Specific
            AntiSymmetryAir = 15,  // Maxwell Specific
            Cut             = 16,  // Maxwell Specific
            Boundary        = 17,  // Maxwell Specific
            ThinShell       = 18,  // Maxwell Specific
            Ghost           = 19,  // Maxwell Specific
            UNDEFINED = 20
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
