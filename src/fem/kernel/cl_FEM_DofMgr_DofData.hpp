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
                const proc_t mCommRank ;
                const proc_t mCommSize ;

                //! DOFs used by this proc
                Cell< Dof * > mDOFs;

                //! Abstract dofs, needed for fem controller
                Cell< Dof * > mAbstractDOFs;

                //! Hanging DOFs
                Cell< Dof * > mHangingDOFs ;

                //! this map converts doftypes to field indices on the mesh
                Map< index_t , index_t > mDofTypeToField ;

                //! this map links the dofs to a unique identifier
                Map< id_t, Dof * > mDofMap ;

                id_t mEdgeDofOffset   = gNoID ;
                id_t mFaceDofOffset   = gNoID ;
                id_t mCellDofOffset   = gNoID ;
                id_t mLambdaDofOffset = gNoID ;

                Cell< Vector< index_t > > mDofIndexTables ;

                // system wide dofs. Initialized here: SolverData binds const
                // references to these at DofManager construction, BEFORE
                // create_dofs() runs — an early read must see 0, not garbage
                index_t mNumberOfFreeDofs = 0 ;
                index_t mNumberOfFixedDofs = 0 ;
                index_t mNumberOfHangingDofs = 0 ;
                index_t mMyNumberOfFreeDofs = 0 ;
                index_t mMyNumberOfFixedDofs = 0 ;
                index_t mMyNumberOfHangingDofs = 0 ;
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

//------------------------------------------------------------------------------

                /**
                 * expose the hanging dof container
                 */
                Cell< Dof * > &
                hanging_dofs();

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

                //! copy values between dofs and fields. Free dofs read
                //! field -> dof; fixed dofs write dof -> field, unless
                //! aFreeDofsOnly is set ( seeding mode: the field is the
                //! truth, e.g. a restored memdump, and a fixed dof may
                //! still hold a factory dummy )
                void
                init_dof_values( const bool aFreeDofsOnly = false );

//------------------------------------------------------------------------------

                void
                split_dof_container(
                    Cell< graph::Vertex * > & aFreeDofs,
                    Cell< graph::Vertex * > & aFixedDofs );

//------------------------------------------------------------------------------

                void
                restore_dof_container(
                    Cell< graph::Vertex * > & aFreeDofs,
                    Cell< graph::Vertex * > & aFixedDofs );

//------------------------------------------------------------------------------

                void
                synchronize_dirichlet_bcs();

//------------------------------------------------------------------------------

                void
                reorder_dofs(
                    const Vector< id_t > & aGraphData,
                    Graph & aFreeDofs,
                    Graph & aFixedDofs );

//------------------------------------------------------------------------------

                const index_t &
                number_of_free_dofs() const ;

//------------------------------------------------------------------------------

                const index_t &
                number_of_fixed_dofs() const ;

                const index_t &
                number_of_hanging_dofs() const ;

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
                dof_indices( const uint aProc ) const;

//------------------------------------------------------------------------------

                void
                connect_dofs_to_mesh();

//------------------------------------------------------------------------------

                void
                disconnect_dofs_from_mesh();

//------------------------------------------------------------------------------

                void
                create_dofwise_t_matrices_master();

//-------------------------------------------------------------------------------

                void
                collect_hanging_dofs();

//-------------------------------------------------------------------------------

                Cell< Dof * > &
                abstract_dofs();

//-------------------------------------------------------------------------------

                void
                extract_abstract_dofs_from_mesh();

//-------------------------------------------------------------------------------

                const index_t &
                my_number_of_free_dofs() const ;

//-------------------------------------------------------------------------------

                const index_t &
                my_number_of_fixed_dofs() const ;
//-------------------------------------------------------------------------------

                const index_t &
                my_number_of_hanging_dofs() const ;

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

                void
                remove_hanging_dofs_from_container();

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

            inline Cell< Dof * > &
            DofData::hanging_dofs()
            {
                return mHangingDOFs ;
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

            inline const index_t &
            DofData::number_of_hanging_dofs() const
            {
                return mNumberOfHangingDofs ;
            }
//------------------------------------------------------------------------------

            inline const Vector< index_t > &
            DofData::dof_indices( const uint aProc ) const
            {
                BELFEM_ASSERT( mDofIndexTables.size() > 0,
                    "DofTables have not been allocated yet");

                return mDofIndexTables( aProc );
            }

//------------------------------------------------------------------------------

            template< typename T >
            void
            allocate_dof_containers( Cell< T * > aBasis )
            {
                for ( T * tBasis : aBasis )
                {
                    tBasis->allocate_dof_container();
                }
            }
//------------------------------------------------------------------------------

            inline Cell< Dof * > &
            DofData::abstract_dofs()
            {
                BELFEM_ASSERT( mCommRank == 0, "Abstract dofs can only be accessed by master");
                return mAbstractDOFs;
            }

//-------------------------------------------------------------------------------

            inline const index_t &
            DofData::my_number_of_free_dofs() const
            {
                return mMyNumberOfFreeDofs;
            }

//-------------------------------------------------------------------------------

            inline const index_t &
            DofData::my_number_of_fixed_dofs() const
            {
                return mMyNumberOfFixedDofs;
            }

//-------------------------------------------------------------------------------

            inline const index_t &
            DofData::my_number_of_hanging_dofs() const
            {
                return mMyNumberOfHangingDofs;
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */


//------------------------------------------------------------------------------


    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_FEM_DOFMGR_DOFDATA_HPP
