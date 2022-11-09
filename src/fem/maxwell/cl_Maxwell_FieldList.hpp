//
// Created by christian on 12/20/21.
//

#ifndef BELFEM_CL_MAXWELL_FIELDLIST_HPP
#define BELFEM_CL_MAXWELL_FIELDLIST_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "en_FEM_DomainType.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell ;

        namespace maxwell
        {
            /**
             * the dof list is a struct that contains the dof lists
             * for blocks, sidesets and interfaces
             */
            class FieldList
            {
            public:
//-----------------------------------------------------------------------------
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // these fields actually belong to the parent IWG.
                // we relink them here fore better code readability
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                //! this map associates unique numbers with the dofs
                Map< string, index_t > & mDofMap ;

                //! all dofs, points to mDofFields of parent IWG
                Cell< string > & Dofs ;

                //! non dof fields
                Cell< string > & NonDof ;

                //! fields that are hidden from exodus
                Cell< string > & Hidden ;

                //! all fields
                Cell< string > & All ;

                //! langrange multiplicators
                Cell< string > Lambda ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // main dofs, set by child IWG
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                //! dofs per superconducting block
                Cell< string > Superconductor ;

                //! dofs per coil
                Cell< string > Coil ;

                //! dofs per iron block
                Cell< string > Ferro ;

                //! dofs per air
                Cell< string > Air ;

                //! dofs along Cut-Interface
                Cell< string > Cut ;

                //! dofs along superconducting shell
                //! must contain edge dofs, face dofs and lambda dofs for air interface
                Cell< string > ThinShell ;

                //! dofs per iron block
                Cell< string > FerroLast ;

                Cell< string > ThinShellLast ;

                //! dofs on Ghost Elements
                Cell< string > Ghost ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // interface dofs, set by child IWG
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                //! dofs along air-superconductor-interface
                Cell< string > InterfaceScAir ;

                //! dofs along superconductor-ferro-interface
                Cell< string > InterfaceScFm ;

                Cell< string > InterfaceFmFm ;

                Cell< string > InterfaceFmAir ;

                Cell< string > InterfaceThinShell ;

                Cell< string > SymmetryAir ;
                Cell< string > SymmetryFerro ;
                Cell< string > SymmetryConductor ;

                Cell< string > AntiSymmetryAir ;
                Cell< string > AntiSymmetryFerro ;
                Cell< string > AntiSymmetryConductor ;

                Cell< string > BoundaryAir ;
                Cell< string > BoundaryFerro ;
                Cell< string > BoundarySc ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // non-dof Fields
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                //! magnetic field density
                Cell< string > MagneticFieldDensity ;

                //! current density
                Cell< string > CurrentDensity ;

                //! current bc for cut
                Cell< string > CurrentBC ;


                //! special list for printout
                Cell< string > Labels ;

//-----------------------------------------------------------------------------
            private:
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
            public:
//-----------------------------------------------------------------------------

                FieldList(
                        Map < string, uint > & aDofMap,
                        Cell< string > & aDofs,
                        Cell< string > & aNonDof,
                        Cell< string > & aHidden,
                        Cell< string > & aAll );

                ~FieldList() = default ;

//-----------------------------------------------------------------------------

                /**
                 * This subroutine is called by the initialize function of the IWG
                 */
                void
                initialize( IWG_Maxwell * aIWG );

//-----------------------------------------------------------------------------
                /**
                 * This subroutine is called by the initialize function of the IWG
                 * It assigns the block-specific dofs based on the given IDs
                 * @param aBlockIDs           ids of selected blocks
                 * @param aBlockTypeMap       the types that are associated with each block
                 * @param aBlockDofs          the dofs that are associated with each block
                 */
                void
                collect_block_dofs( const Vector< id_t >             & aBlockIDs,
                                    const Map< id_t, DomainType >    & aBlockTypeMap,
                                          Cell< Vector< index_t > >  & aBlockDofs );

//-----------------------------------------------------------------------------

                void
                collect_sideset_dofs(
                        const Vector< id_t >              & aSideSetIDs,
                        const Map< id_t, DomainType >     & aSideSetTypeMap,
                        const Map< id_t, MagfieldBcType > & aSideSetSubTypeMap,
                        Cell< Vector< index_t > >         & aSideSetDofs );

//-----------------------------------------------------------------------------
            private:
//-----------------------------------------------------------------------------

                /**
                 * converts the selected cell of strings into the dof indices
                 * @param aDofList
                 * @param aDofTable
                 */
                void
                create_doftable( const Cell< string > & aDofList, Vector< index_t > & aDofTable );

            };
        } /* end namespace maxwell */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_MAXWELL_FIELDLIST_HPP
