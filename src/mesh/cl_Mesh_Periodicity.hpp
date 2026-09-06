/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MESH_PERIODICITY_HPP
#define BELFEM_CL_MESH_PERIODICITY_HPP

#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Face.hpp"
#include "cl_SideSet.hpp"
#include "cl_Bitset.hpp"

namespace belfem
{
    class Mesh;


    namespace mesh
    {
        class PeriodicityFactory ;

        /**
         * Stores matched periodic boundary entity pairs (nodes, edges, faces).
         * For each entity type, master(i) corresponds to slave(i).
         *
         * Node pairs are populated by PeriodicityFactory. Edge and face pairs
         * are derived from the node mapping by update().
         */
        class Periodicity
        {
            Mesh * mMesh;

            // matched node pairs
            Cell< Node * > mMasterNodes ;
            Cell< Node * > mSlaveNodes ;

            // matched edge pairs
            Cell< Edge * > mMasterEdges ;
            Cell< Edge * > mSlaveEdges ;

            // matched face pairs
            Cell< Face * > mMasterFaces ;
            Cell< Face * > mSlaveFaces ;

            // matched facet pairs
            Cell< Facet * > mMasterFacets ;
            Cell< Facet * > mSlaveFacets ;

            // the three nodes that define the source plane
            Cell< Node * > mMasterPlane ;

            // the three nodes that define the target plane
            Cell< Node * > mSlavePlane ;

            // the hesse forms
            Vector< real > mMasterHesse ;
            Vector< real > mSlaveHesse ;

            Bitset< 8 > mPopulated ;

            Cell< std::pair< Node * , Node * > > mNodePairBackup ;

            //! set once restore_node_pairs() has consumed the backup; guards
            //! against a second geometric re-match of a mesh that contains
            //! coincident duplicate nodes
            bool mNodePairsRestored = false ;
            real mMaxTolerance = BELFEM_MESH_EPSILON ;

            friend class PeriodicityFactory;

            public:

                Periodicity( Mesh * aMesh );

                ~Periodicity();

                Cell< Node * > &
                master_nodes();

                Cell< Node * > &
                slave_nodes();

                Cell< Edge * > &
                master_edges();

                Cell< Edge * > &
                slave_edges();

                Cell< Face * > &
                master_faces();

                Cell< Face * > &
                slave_faces();

                Cell< Facet * > &
                master_facets();

                Cell< Facet * > &
                slave_facets();

                Cell< Node * > &
                master_plane();

                Cell< Node * > &
                slave_plane();

                /**
                 * Re-runs the factory on this object: resets all pairs,
                 * re-collects facets and nodes, restores the node pairs from
                 * the backup ( or re-matches them ), then derives edge, facet
                 * and face pairs. Edge and face pairs are only built if the
                 * mesh has edges / faces.
                 */
                void
                update();

                void
                set_entity_dependencies();

                void
                reset_nodes();

                void
                reset_edges();

                void
                reset_faces();

                void
                reset_facets();

                bool
                is_flagged( const EntityType aType ) const ;

                void
                backup_node_pairs();

                void
                restore_node_pairs();

                void
                clear_node_pair_backup();

                void
                add_node_pair_to_backup(
                    Node * aNodeA,
                    Node * aNodeB,
                    const real aTolerance = BELFEM_MESH_EPSILON );

                const Vector< real > &
                master_hesse() const ;

                const Vector< real > &
                slave_hesse() const ;

                bool
                has_node_pair_backup() const;

                bool
                node_pairs_restored() const;

            protected:

                void
                flag( const EntityType aType );

                void
                unflag( const EntityType aType );

                Vector< real > &
                master_hesse() ;

                Vector< real > &
                slave_hesse() ;

                int
                pair_orientation( Node * aNodeA, Node * aNodeB, const real aTolerance = BELFEM_MESH_EPSILON  ) const;

        };

        inline Cell< Node * > &
        Periodicity::master_nodes()
        {
            return mMasterNodes ;
        }

        inline Cell< Node * > &
        Periodicity::slave_nodes()
        {
            return mSlaveNodes ;
        }

        inline Cell< Edge * > &
        Periodicity::master_edges()
        {
            return mMasterEdges ;
        }

        inline Cell< Edge * > &
        Periodicity::slave_edges()
        {
            return mSlaveEdges ;
        }

        inline Cell< Face * > &
        Periodicity::master_faces()
        {
            return mMasterFaces ;
        }

        inline Cell< Face * > &
        Periodicity::slave_faces()
        {
            return mSlaveFaces ;
        }

        inline Cell< Facet * > &
        Periodicity::master_facets()
        {
            return mMasterFacets ;
        }

        inline Cell< Facet * > &
        Periodicity::slave_facets()
        {
            return mSlaveFacets ;
        }

        inline Cell< Node * > &
        Periodicity::master_plane()
        {
            return mMasterPlane ;
        }

        inline Cell< Node * > &
        Periodicity::slave_plane()
        {
            return mSlavePlane ;
        }

        inline bool
        Periodicity::is_flagged( const EntityType aType ) const
        {
            return mPopulated.test( static_cast< uint >( aType ) );
        }

        inline void
        Periodicity::flag( const EntityType aType )
        {
            mPopulated.set( static_cast< uint >( aType ) );
        }

        inline void
        Periodicity::unflag( const EntityType aType )
        {
            mPopulated.reset( static_cast< uint >( aType ) );
        }

        inline Vector< real > &
        Periodicity::master_hesse()
        {
            return mMasterHesse ;
        }

        inline Vector< real > &
        Periodicity::slave_hesse()
        {
            return mSlaveHesse ;
        }

        inline const Vector< real > &
        Periodicity::master_hesse() const
        {
            return mMasterHesse ;
        }

        inline const Vector< real > &
        Periodicity::slave_hesse() const
        {
            return mSlaveHesse ;
        }

        inline
        int
        Periodicity::pair_orientation( Node * aNodeA, Node * aNodeB, const real aTolerance ) const
        {
            real tAM =     std::abs( aNodeA->x() * mMasterHesse( 0 )
                                     + aNodeA->y() * mMasterHesse( 1 )
                                     + aNodeA->z() * mMasterHesse( 2 )
                                     - mMasterHesse( 3 ) );

            real tBM =     std::abs( aNodeB->x() * mMasterHesse( 0 )
                                     + aNodeB->y() * mMasterHesse( 1 )
                                     + aNodeB->z() * mMasterHesse( 2 )
                                     - mMasterHesse( 3 ) );

            real tAS =     std::abs( aNodeA->x() * mSlaveHesse( 0 )
                                     + aNodeA->y() * mSlaveHesse( 1 )
                                     + aNodeA->z() * mSlaveHesse( 2 )
                                     - mSlaveHesse( 3 ) );

            real tBS =     std::abs( aNodeB->x() * mSlaveHesse( 0 )
                                     + aNodeB->y() * mSlaveHesse( 1 )
                                     + aNodeB->z() * mSlaveHesse( 2 )
                                     - mSlaveHesse( 3 ) );


            if ( tAM < aTolerance && tBS < aTolerance )
            {
                return 1 ;
            }
            if ( tAS < aTolerance && tBM < aTolerance )
            {
                return -1 ;
            }
            return 0 ;
        }

        inline bool
        Periodicity::has_node_pair_backup() const
        {
            return mNodePairBackup.size() > 0 ;
        }

        inline bool
        Periodicity::node_pairs_restored() const
        {
            return mNodePairsRestored ;
        }

    }
}
#endif //BELFEM_CL_MESH_PERIODICITY_HPP