//
// Created by christian on 7/9/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_DOFDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_DOFDATA_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"
#include "cl_Bitset.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_DofMgr_Parameters.hpp"
#include "cl_IWG.hpp"

namespace belfem
{
    class Mesh ;

    namespace fem
    {
        class Dof ;
        class DofManager ;
        class Kernel ;

        namespace dofmgr
        {
            class DofData
            {
                //! the parent object
                DofManager * mParent;

                //! the kernel
                Kernel * mKernel ;

                //! the mesh this problem runs on
                Mesh * mMesh ;

                //! parameter object
                Parameters * mParams ;

                //! the rank of this proc
                const proc_t mMyRank ;

                //! DOFs used by this proc
                Cell< Dof * > mDOFs;

                //! this map converts doftypes to field indices on the mesh
                Map< index_t , index_t > mDofTypeToField ;

                //! this map links the dofs to a unique identifier
                Map< id_t, Dof * > mDofMap ;

                id_t mEdgeDofOffset   = gNoID ;
                id_t mFaceDofOffset   = gNoID ;
                id_t mCellDofOffset   = gNoID ;
                id_t mLambdaDofOffset = gNoID ;

                // DOFs used by each proc ( master only )
                Cell< Vector< id_t > > mDofIDs;

                // indices of dofs per proc, with respect to master
                Cell< Vector< index_t > > mDofIndices;

                // system wide dofs
                index_t mNumberOfFreeDofs;
                index_t mNumberOfFixedDofs;
                id_t mNumDofTypes = BELFEM_UINT_MAX ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                DofData( DofManager * aParent,
                         Parameters * aParams ) ;

//------------------------------------------------------------------------------

                ~DofData() ;

//------------------------------------------------------------------------------

                void
                create_dofs( IWG * aIWG );

//------------------------------------------------------------------------------

                void
                create_field_map( IWG  * aIwg );

//------------------------------------------------------------------------------

                /**
                 * expose the dof container
                 */
                 Cell< Dof * > &
                 dofs();

//----------------------------------------------------------------------------

                /**
                 * check if a dof exists, needed for bearing creation
                 */
                 bool
                 dof_exists( const id_t aID ) const;

//------------------------------------------------------------------------------

                index_t
                doftype_to_field_index( const index_t aDofType );

//------------------------------------------------------------------------------

                // access one dof by its id
                Dof *
                dof( const id_t aID );

//------------------------------------------------------------------------------

                id_t
                node_dof_id(
                        const id_t aNodeID,
                        const uint aDofType )  const ;

//------------------------------------------------------------------------------

                id_t
                edge_dof_id(
                        const id_t aEdgeID,
                        const uint aDofType )  const ;

//------------------------------------------------------------------------------

                id_t
                face_dof_id(
                        const id_t aFaceID,
                        const uint aDofType )  const ;

//------------------------------------------------------------------------------

                id_t
                cell_dof_id(
                        const id_t aCellID,
                        const uint aDofType )  const ;

//------------------------------------------------------------------------------

                id_t
                lambda_dof_id(
                        const id_t aFacetID,
                        const uint aDofType )  const ;

//------------------------------------------------------------------------------

                void
                init_dof_values();

//------------------------------------------------------------------------------

                void
                init_dirichlet_bcs();

//------------------------------------------------------------------------------

                void
                compute_dof_indices();

//------------------------------------------------------------------------------

                const index_t &
                number_of_free_dofs() const ;

//------------------------------------------------------------------------------

                const index_t &
                number_of_fixed_dofs() const ;

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------

                /**
                 * computes the number of dofs per element on block
                 */
                 uint
                 num_dofs_per_element( const id_t aBlockID ) const;

//------------------------------------------------------------------------------

                /**
                 * computes the number of dofs per element on sideset
                 */
                uint
                num_dofs_per_facet( const id_t aSideSetID ) const;

//------------------------------------------------------------------------------

                const Vector< index_t > &
                dof_indices( const uint aProcIndex ) const;

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

                void
                create_dof_map();

//------------------------------------------------------------------------------

                /**
                 * help function to be called by master
                 */
                void
                compute_max_ids( Vector< id_t > & aMaxEntityIDs );

//-------------------------------------------------------------------------------

                index_t
                count_node_dofs( IWG  * aIWG,
                                 Vector< id_t > & aEntityIDs,
                                 Vector< index_t > & aDofTypes,
                                 Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags );

//-----------------------------------------------------------------------------

                index_t
                count_edge_dofs( IWG  * aIWG,
                                 Vector< id_t > & aEntityIDs,
                                 Vector< index_t > & aDofTypes,
                                 Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags );

//-----------------------------------------------------------------------------

                index_t
                count_face_dofs( IWG  * aIWG,
                                 Vector< id_t > & aEntityIDs,
                                 Vector< index_t > & aDofTypes,
                                 Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags );

//-----------------------------------------------------------------------------

                index_t
                count_cell_dofs( IWG  * aIWG,
                                 Vector< id_t > & aEntityIDs,
                                 Vector< index_t > & aDofTypes,
                                 Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags );

//-----------------------------------------------------------------------------

                index_t
                count_lambda_dofs( IWG  * aIWG,
                                 Vector< id_t > & aEntityIDs,
                                 Vector< index_t > & aDofTypes,
                                 Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags );

//-----------------------------------------------------------------------------

                index_t
                count_dofs_for_proc(
                        const index_t                              aProc,
                        const Vector< id_t >                     & aDofIDs,
                        const Vector< id_t >                     & aEntityIDs,
                        const Vector< index_t >                  & aDofTypes,
                        const Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags,
                        Vector< id_t >                           & aProcDofIDs,
                        Vector< id_t >                           & aProcEntityIDs,
                        Vector< index_t >                        & aProcDofTypes );

//-------------------------------------------------------------------------------

                /**
                 * a help function that initializes the offsets for the
                 * dof computation function calculate_dof_id()
                 */
                void
                compute_dof_offsets(  IWG  * aIWG  );

//-------------------------------------------------------------------------------

                /**
                 * check which IDs are used on this proc and send the info
                 * to the master
                 */
                 void
                 create_dof_table();

//-------------------------------------------------------------------------------
                /**
                 * identifies which dofs are connected to the current dof
                 * and connects the graph
                 */
                 //void
                 //compute_adjacency(
                 //        Dof * aDof,
                 //        Cell< Dof * > & aWork,
                 //        const bool aFixedFlag );

//-------------------------------------------------------------------------------
            };

//------------------------------------------------------------------------------

            inline Cell< Dof * > &
            DofData::dofs()
            {
                return mDOFs ;
            }

//------------------------------------------------------------------------------

            inline Dof *
            DofData::dof( const id_t aID )
            {
                return mDofMap( aID );
            }

//------------------------------------------------------------------------------

            inline bool
            DofData::dof_exists( const id_t aID ) const
            {
                return mDofMap.key_exists( aID );
            }

//------------------------------------------------------------------------------

            inline index_t
            DofData::doftype_to_field_index( const index_t aDofType )
            {
                return mDofTypeToField( aDofType );
            }

//------------------------------------------------------------------------------

            inline id_t
            DofData::node_dof_id(
                    const id_t aNodeID,
                    const uint aDofType )  const
            {
                return aNodeID * mNumDofTypes + aDofType ;
            }

//------------------------------------------------------------------------------

            inline id_t
            DofData::edge_dof_id(
                    const id_t aEdgeID,
                    const uint aDofType )  const
            {
                return mEdgeDofOffset + aEdgeID * mNumDofTypes + aDofType ;
            }

//------------------------------------------------------------------------------

            inline id_t
            DofData::face_dof_id(
                    const id_t aFaceID,
                    const uint aDofType )  const
            {
                return mFaceDofOffset + aFaceID * mNumDofTypes + aDofType ;
            }

//------------------------------------------------------------------------------

            inline id_t
            DofData::cell_dof_id(
                    const id_t aCellID,
                    const uint aDofType )  const
            {
                return mCellDofOffset + aCellID * mNumDofTypes + aDofType ;
            }

//------------------------------------------------------------------------------

            inline id_t
            DofData::lambda_dof_id(
                    const id_t aFacetID,
                    const uint aDofType )  const
           {
                return mLambdaDofOffset + aFacetID * mNumDofTypes + aDofType ;
           }

//------------------------------------------------------------------------------

            inline const index_t &
            DofData::number_of_free_dofs() const
            {
                return mNumberOfFreeDofs ;
            }

//------------------------------------------------------------------------------

            inline const index_t &
            DofData::number_of_fixed_dofs() const
            {
                return mNumberOfFixedDofs ;
            }
//------------------------------------------------------------------------------

            inline const Vector< index_t > &
            DofData::dof_indices( const uint aProcIndex ) const
            {
                return mDofIndices( aProcIndex );
            }

//------------------------------------------------------------------------------

        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_DOFDATA_HPP
