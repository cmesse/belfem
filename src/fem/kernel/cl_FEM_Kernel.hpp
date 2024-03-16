//
// Created by Christian Messe on 24.10.19.
//

#ifndef BELFEM_CL_FEM_KERNEL_HPP
#define BELFEM_CL_FEM_KERNEL_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Material.hpp"
#include "en_Materials.hpp"

#include "cl_FEM_Dof.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_KernelParameters.hpp"

#include "cl_FEM_DofManager.hpp"

#include "en_IWGs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class IwgFactory ;

        class Kernel
        {
            // pointer to parameter list
            KernelParameters * mParams;

            // pointer to mesh
            Mesh             * mMesh;

            // rank of this proc
            const proc_t     mMyRank;

            // rank of master proc
            const proc_t     mMasterRank;

            // how many procs contribute to this kernel
            const proc_t     mNumberOfProcs;


            // offset for fields, in case fields already exist on mesh
            // when Kernel is generated
            const uint mFieldOffset ;

            // flag telling if the kernel destroys the parameters on exit
            bool mOwnParameters = false ;

            Vector< proc_t >    mCommTable;
            Map< proc_t, uint > mCommMap;

            proc_t           mMyCommIndex;

            Vector< index_t > mNumElementsPerBlock;

            // submesh
            Mesh * mSubMesh = nullptr;
            Map< id_t, mesh::Node*    > mNodeMap;
            Map< id_t, mesh::Element* > mElementMap;
            Map< id_t, mesh::Edge * >   mEdgeMap ;
            Map< id_t, mesh::Face * >   mFaceMap ;
            Map< id_t, mesh::Element* > mVertexMap;
            Map< id_t, mesh::Facet * >  mFacetMap ;
            // Map< id_t, mesh::Facet * >  mConnectorMap ;

            Cell< DofManager * > mDofManagers ;
            Cell< IWG * >        mIWGs ;

            Cell< Material * > mMaterials;
            Map< MaterialType, Material * > mMaterialMap;

            // list with nodes per proc, master only
            Cell< Vector< index_t > > mNodeTable;
            Cell< Vector< index_t > > mElementTable;
            Cell< Vector< index_t > > mEdgeTable;
            Cell< Vector< index_t > > mFaceTable;
            Cell< Vector< index_t > > mFacetTable;
            Cell< Vector< index_t > > mConnectorTable;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Kernel( KernelParameters * aKernelParameters );

//------------------------------------------------------------------------------

            ~Kernel();

//------------------------------------------------------------------------------

            /**
             * get the rank of the master
             */
             proc_t
             master();

//------------------------------------------------------------------------------

            /**
             * expose parameter object
             */
             const  KernelParameters *
             params();

//------------------------------------------------------------------------------

            /**
             * expose communication table
             */
             const Vector< proc_t > &
             comm_table() const;

            /**
             * get id of individual proc
             */
             proc_t
             comm_table( const index_t aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * expose a field
             */
            DofManager *
            dofmgr( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * returns the material object. Creates one if it didn't exist
             * already. This function is called by a block during Block::set_material
             */
             Material *
             get_material( const MaterialType aMaterial );

//------------------------------------------------------------------------------

            /**
             * manually add a material to the kernel. The material is now owned
             * and destroyed by the kernel
             */
            void
            add_material( Material * aMaterial );

//------------------------------------------------------------------------------

            /**
             * return the list of nodes from a specific proc other than master
             */
             const Vector< index_t > &
             node_table( const uint aProcIndex );

//------------------------------------------------------------------------------

            /**
             * return the list of elements from a specific proc other than master
             */
            const Vector< index_t > &
            element_table( const uint aProcIndex );

//------------------------------------------------------------------------------

            /**
             * return the list of edges from a specific proc other than master
             */
             const Vector< index_t > &
             edge_table( const uint aProcIndex );

//------------------------------------------------------------------------------

            /**
             * return the list of faces from a specific proc other than master
             */
            const Vector< index_t > &
            face_table( const uint aProcIndex );

//------------------------------------------------------------------------------

             /**
              * return the list of facets from a specific proc other than master
              */
             const Vector< index_t > &
             facet_table( const uint aProcIndex );

//------------------------------------------------------------------------------

             const proc_t &
             number_of_procs() const ;

//------------------------------------------------------------------------------

            /**
             * returns the number of fields that existed before Kernel was created
             */
            const uint &
            field_offset() const ;

//------------------------------------------------------------------------------

            /**
             * returns the mesh on the master and the submesh on other
             * @param aTarget
             */
            Mesh *
            mesh();

//------------------------------------------------------------------------------

            IWG *
            create_equation( const IwgType aEquationType,
                             const Vector< id_t > aBlocks={} );

//------------------------------------------------------------------------------

            /**
             * special function to add an equation that has already
             * been created. The equation is now owned and destroyed
             * by the kernel
             */
            void
            add_equation( IWG * aEquation );

//------------------------------------------------------------------------------

            DofManager *
            create_field( IWG * aEquation );

//------------------------------------------------------------------------------

            bool
            is_master() const ;

//------------------------------------------------------------------------------

            /**
             * if this flag is set, the kernel will destroy the parameters
             */
             void
             claim_parameter_ownership( const bool aFlag = true );


//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            //void
            //create_communication_list();

//------------------------------------------------------------------------------

            /**
             * tells a proc if it takes part in this kernel
             */
            void
            activate_procs();

//------------------------------------------------------------------------------

            void
            symrcm();

//------------------------------------------------------------------------------

            /**
             * count the number of elements on each block plus aura
             */
            void
            distribute_mesh();

//------------------------------------------------------------------------------

            void
            send_submesh( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            send_nodes( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            send_elements( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            send_facets( const proc_t & aTarget, const bool aConnectorSwitch = false );

//------------------------------------------------------------------------------

            void
            send_vertices( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            send_edges( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            send_faces( const proc_t & aTarget );

//------------------------------------------------------------------------------

            void
            receive_edges();

//------------------------------------------------------------------------------

            void
            receive_faces();

//------------------------------------------------------------------------------

            void
            receive_submesh();

//------------------------------------------------------------------------------

            void
            receive_nodes();

//------------------------------------------------------------------------------

            void
            receive_elements();

//------------------------------------------------------------------------------

            void
            receive_facets( const bool aConnectorSwitch = false  );
            
//------------------------------------------------------------------------------

            void
            receive_vertices();

//------------------------------------------------------------------------------

            /**
             * free some memory
             */
            void
            delete_maps();

//------------------------------------------------------------------------------

            /**
             * this step is needed to get the field synchronization working
             */
            void
            merge_facet_and_connector_tables();

//------------------------------------------------------------------------------

            /**
             * creates the dof connectivity graph that is used for the
             * dof numeration
             * @param aGraph
             */
            void
            create_dof_graph( Cell< Dof * > & aGraph );

//------------------------------------------------------------------------------

            /**
             * returns the index of the proc on the comm table
             * @param aProcID
             * @return
             */
            uint
            proc_index( const proc_t & aProcID );

//------------------------------------------------------------------------------

            /**
             * help function for parallel edge and face creation.
             * this function grabs elements based on gived blocks and
             * sidesets
             *
             * @param aBlockIDs
             * @param aSideSetIDs
             * @param aElements
             */
            void
            collect_elements(
                    const Vector< id_t >    & aBlockIDs,
                    const Vector< id_t >    & aSideSetIDs,
                    Cell< mesh::Element * > & aElements );

//------------------------------------------------------------------------------

            // make sure that all elements have positive volume
            void
            compute_element_volumes();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const proc_t &
        Kernel::number_of_procs() const
        {
            return mNumberOfProcs ;
        }
//------------------------------------------------------------------------------

        inline const uint &
        Kernel::field_offset() const
        {
            return mFieldOffset ;
        }

//------------------------------------------------------------------------------

        inline Mesh *
        Kernel::mesh()
        {
            return mMyRank == mMasterRank ? mMesh : mSubMesh ;
        }

//------------------------------------------------------------------------------

        inline uint
        Kernel::proc_index( const proc_t & aProcID )
        {
            return mCommMap( aProcID );
        }

//------------------------------------------------------------------------------

        inline bool
        Kernel::is_master() const
        {
            return mMyRank == mMasterRank ;
        }

//------------------------------------------------------------------------------

    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_FEM_KERNEL_HPP
