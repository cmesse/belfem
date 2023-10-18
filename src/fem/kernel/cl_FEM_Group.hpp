//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_CL_FEM_GROUP_HPP
#define BELFEM_CL_FEM_GROUP_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Mesh.hpp"
#include "meshtools.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_Material.hpp"
#include "en_FEM_DomainType.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    class Material;
    enum class MaterialType;

    namespace fem
    {
        class DofManagerBase;
        class Element;
        class Block ;
        class BoundaryCondition ;

//------------------------------------------------------------------------------

        class Group
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // pointer to parent
            DofManagerBase * mParent ;

            // pointer to swap object
            Calculator * mCalc = nullptr ;

            // either block or sideset
            const GroupType mType ;

            const ElementType mElementType;

            const id_t mID;

            const index_t mNumberOfElements;

            // flag that tells if elements are destoyed by destructor
            const bool mOwnElements ;

            const proc_t  mMyRank;

            // detailed domain type, mainly used for Maxwell
            DomainType mDomainType = DomainType::Default ;

            // container for node coordinates
            Matrix< real > mNodeCoords ;

            // must be set by child
            uint mNumberOfNodesPerElement = BELFEM_UINT_MAX ;

            bool mIsIsogeometric = true ;

            IntegrationData * mIntegrationData = nullptr ;
            IntegrationData * mGeometryIntegrationData = nullptr ;

            // shape function for element boundaries
            Cell< Cell< Matrix < real > > > mBoundaryN ;

            // edge shape function
            Matrix< real >         mWorkE ;
            Matrix< real >         mWorkDEDXi;
            Matrix< real >         mWorkDEDX;

            // Element container
            Cell< Element * > mElements;

            Matrix< real > mWorkDNDX;    // work matrix for dNdX

            // Work matrices for geometry transformation
            // remember index for which X was evaluate
            uint           mWorkIndex = BELFEM_UINT_MAX;


            Matrix< real > mWorkG;   // work matrix for Jacobian ( dof x dof ) ( eg TATCAD element )
            Matrix< real > mWorkH;   // work matrix for Jacobian ( dof x dof ) ( eg TATCAD element )


            Matrix< real > mWorkinvJ; // work matrix for inverse Jacobian  ( Nédélec )
            real mWorkDetInvJ;        // contains determinant of inverse Jacobian
            real mWorkDetJ;           // contains abs of det of  Jacobian

            Vector< real > mWorkRhs;  // work vector for rhs

            Matrix< real > mWorkB;    // work matrix for B-Matrix
            Matrix< real > mWorkC;    // work matrix for C-Matrix ( curl h )
            Matrix< real > mWorkD;    // work matrix for C-Matrix ( curl a )

            Matrix< real > mWorkJ;    // work matrix for first derivatives
            Matrix< real > mWorkK;    // work matrix for second derivatives
            Matrix< real > mWorkL;    // work matrix for second derivatives
            Matrix< real > mWorkM;    // work matrix for masses

            Matrix< real > mWorkN;      // work matrix for N-Matrix

            Matrix< real > mWorkX;    // work matrix for special purpose node coordinates
            Matrix< real > mWorkXm;   // work matrix for special purpose node coordinates
            Matrix< real > mWorkXs;   // work matrix for special purpose node coordinates

            Vector< real > mWorkphi;  // work vector for nodal field data
            Matrix< real > mWorkPhi;  // work matrix for nodal field data

            Vector< real > mWorkpsi; // special purpose vector, assign by iwg child
            Matrix< real > mWorkPsi; // special purpose matrix, assign by iwg child

            Vector< real > mWorkchi; // special purpose vector, assign by iwg child
            Matrix< real > mWorkChi; // special purpose matrix, assign by iwg child

            Vector< real > mWorksigma; // special purpose vector, assign by iwg child
            Matrix< real > mWorkSigma; // special purpose matrix, assign by iwg child

            Vector< real > mWorktau; // special purpose vector, assign by iwg child
            Matrix< real > mWorkTau; // special purpose matrix, assign by iwg child

            Vector< real > mWorkgamma; // special purpose vector, assign by iwg child
            Matrix< real > mWorkGamma; // special purpose matrix, assign by iwg child

            Vector< real > mWorktheta; // special purpose vector, assign by iwg child

            Vector< real > mWorkNedelec; // special purpose vector for edge data

            Vector< real > mWorkgeo ; // vector for calculation of geometry Jacobian

            Vector< real > mWorkNormal ;

            // pointer to material ( owned by kernel )
            Material * mMaterial = nullptr;

            // pointer to boundary condition, if set
            const BoundaryCondition * mBoundaryCondition = nullptr ;

            Map< id_t, Element * > mElementMap ; // map to access element by ID

            // pointer to material functions
            void
            ( Group:: * mLinearElasticity )( Matrix< real > & aC, const real aT  ) const;

            void
            ( Group:: * mThermalConductivity )( Matrix< real > & aLambda, const real aT  ) const;

            // empty sidesets and blocks have the id zero.
            // the fake ID helps to acces the underlying objects on the mesh
            id_t mMeshID ;

            //! flag telling if block has a right hand side
            bool mHasRHS = true ;

            //! flag telling if this group is used for the computation
            bool mIsActive = true ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Group(  DofManagerBase      * aParent,
                    const GroupType     aGroupType,
                    const ElementType   aElementType,
                    const        id_t   aID,
                    const     index_t   aNumberOfElements,
                    const        bool   aOwnElements=true );

//------------------------------------------------------------------------------

            virtual ~Group() = default ;

//------------------------------------------------------------------------------

            /**
             * return the type of this group
             */
             GroupType
             type() const ;

//------------------------------------------------------------------------------

            /**
             * return the id of the group
             */
             id_t
             id() const;

//------------------------------------------------------------------------------

            /**
             * return the element type of this block
             */
            virtual ElementType
            element_type() const;

//------------------------------------------------------------------------------

            /**
             * return the field of this block
             */
            DofManagerBase *
            parent();

//------------------------------------------------------------------------------

            /*
             * return the number of elements of this block
             */
            index_t
            number_of_elements() const;

//------------------------------------------------------------------------------

            /*
             * return the number of nodes per element
             */
            uint
            number_of_nodes_per_element() const;

//------------------------------------------------------------------------------

            /*
             * return the number of edges per element
             */
            uint
            number_of_edges_per_element() const;

//------------------------------------------------------------------------------

            /*
             * return the number of edges per element
             */
            uint
            number_of_faces_per_element() const;

//------------------------------------------------------------------------------

            /**
             * expose the integration data
             */
             const IntegrationData *
             integration() ;

//------------------------------------------------------------------------------

            /**
             * expose integration weights
             */
            const Vector< real > &
            integration_weights() const;

//------------------------------------------------------------------------------

            /**
             * expose integration points
             */
            const Matrix< real > &
            integration_points() const;

//------------------------------------------------------------------------------

            virtual const IntegrationData *
            thinshell_integration() const ;

//------------------------------------------------------------------------------

            /*
             * expose element container
             */
            Cell< fem::Element * > &
            elements();

//------------------------------------------------------------------------------

            /**
             * expose calculator object
             */
             Calculator *
             calculator() ;

//------------------------------------------------------------------------------

            /**
             * expose container for node coords
             */
             Matrix< real > &
             node_coords() ;

//------------------------------------------------------------------------------

            /**
             * return the shape function at given point
             */
            const Matrix< real > &
            N( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the shape function at given point for 3D vector field
             */
            const Matrix< real > &
            Nvector( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the shape function at given point, bit as vector
             */
            const Vector< real > &
            n( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the first derivative of shape function at given point
             */
            const Matrix< real > &
            dNdXi( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the second derivative of shape function at given point
             */
            const Matrix< real > &
            d2NdXi2( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the geometry function at given point
             */
            const Matrix< real > &
            G( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the first derivative of geometry function at given point
             */
            const Matrix< real > &
            dGdXi( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * return the second derivative of geometry function at given point
             */
            const Matrix< real > &
            d2GdXi2( const uint aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * the index of the integration point at which the geometry
             * Jacobian was calculated. Remembering this is a precaution to get
             * second derivatives right
             */
            uint &
            work_index();

//------------------------------------------------------------------------------

            /**
             * matrix for edge shape functions
             */
            Matrix< real > &
            work_E();

//------------------------------------------------------------------------------

            /**
             * matrix for edge shape function derivatives
             */
            Matrix< real > &
            work_dEdXi();

//------------------------------------------------------------------------------

            /**
             * matrix for edge shape function derivatives
             */
            Matrix< real > &
            work_dEdX();

//------------------------------------------------------------------------------

            /**
             * special purpose help matrix
             */
            Matrix< real > &
            work_G();

//------------------------------------------------------------------------------

            /**
             * special purpose help matrix
             */
             Matrix< real > &
             work_H();

//------------------------------------------------------------------------------


            /**
             * there are three work matrices X, Y, Z, which are used to
             * transform first and second derivatives onto the geometry space
             * the rules are:
             *
             * J * dNdX   = dNdXi
             * K * d2Ndx2 = d2NdXi2 - L * dNdX
             */
            Matrix< real > &
            work_J();

            Matrix< real > &
            work_K();

            Matrix< real > &
            work_L();

//------------------------------------------------------------------------------

            /**
             * help matrix used for mass matrix
             * @return
             */
            Matrix< real > &
            work_M();

//------------------------------------------------------------------------------

            /**
             * help matrix that stores invese of the Jacobian
             */
            Matrix< real > &
            work_invJ();

//------------------------------------------------------------------------------

            /**
             * container for determinant
             * @return
             */
            real &
            work_det_J();
//------------------------------------------------------------------------------

            /**
             * returns the absolute value of det(J)
             * @return
             */
            real
            dx() const ;

//------------------------------------------------------------------------------

            /**
             * container for inverse of determinant
             * @return
             */
             real &
             work_det_invJ();

//------------------------------------------------------------------------------

            Matrix< real > &
            work_dNdX();

//------------------------------------------------------------------------------
// containers for nodal data
//------------------------------------------------------------------------------

            Vector< real > &
            work_phi();

            Matrix< real > &
            work_Phi();

            Vector< real > &
            work_psi();

            Matrix< real > &
            work_Psi();

            Vector< real > &
            work_gamma();

            Vector< real > &
            work_theta();

            Matrix< real > &
            work_Gamma();

            Vector< real > &
            work_chi();

            Matrix< real > &
            work_Chi();

            Vector< real > &
            work_sigma();

            Matrix< real > &
            work_Sigma();

            Vector< real > &
            work_tau();

            Matrix< real > &
            work_Tau();

            Matrix< real > &
            work_X();

            Matrix< real > &
            work_Xm();

            Matrix< real > &
            work_Xs();

            Vector< real > &
            work_nedelec();

            Vector< real > &
            work_geo();

            Vector< real > &
            work_normal();

//------------------------------------------------------------------------------
// matrices for linear elasticity
//------------------------------------------------------------------------------

            Matrix< real > &
            work_C();


            Matrix< real > &
            work_D();

            Matrix< real > &
            work_N();

            Matrix< real > &
            work_B();

            Vector< real > &
            work_rhs();

//------------------------------------------------------------------------------

            /**
             * access scal
             * ars from the parent
             */
             Vector< real > &
             field_data( const string & aLabel );


//------------------------------------------------------------------------------

            void
            set_material( const MaterialType aMaterial );

//------------------------------------------------------------------------------

            void
            set_material( Material * aMaterial );

//------------------------------------------------------------------------------

            const  Material *
            material() const;

//------------------------------------------------------------------------------

            // needed if you want to access element by id
            void
            create_element_map() ;

//------------------------------------------------------------------------------

            // get an element using its ID instead of index
            Element *
            element( const id_t aID );

//------------------------------------------------------------------------------
// Material Functions
//------------------------------------------------------------------------------

            void
            linear_elasticity(
                    Matrix< real > & aC,
                    const         real aT=BELFEM_TREF ) const;

            void
            thermal_conductivity(
                    Matrix< real > & aLambda,
                    const         real aT=BELFEM_TREF ) const;

//------------------------------------------------------------------------------
// RHS flags
//------------------------------------------------------------------------------

            // flag telling field if an RHS side exists
            bool
            has_rhs() const ;

            // set or unset the rhs flag
            void
            set_rhs_flag( const bool aFlag );

//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

            /**
             * if the sideset or block is empty, ID is zero.
             * use mesh_id if you really need the ID of the corresponding
             * mesh object
             * @param aID
             */
            const id_t &
            mesh_id() const ;

            void
            set_mesh_id( const id_t & aID );

//------------------------------------------------------------------------------

            /**
             * for special purpose integration, sideset only
             */
            virtual const IntegrationData *
            master_integration( const uint aSideSetIndex );

            virtual const IntegrationData *
            slave_integration( const uint aSideSetIndex );

            void
            delete_pointers();

//------------------------------------------------------------------------------

            const InterpolationFunction *
            interpolation_function() const ;

//------------------------------------------------------------------------------

            void
            set_domain_type( const DomainType aType );

//------------------------------------------------------------------------------

            DomainType
            domain_type() const ;

//------------------------------------------------------------------------------

            /**
             * returns type of sideset master
             */
            virtual ElementType
            master_type() const ;
//------------------------------------------------------------------------------

            /**
             * returns type of sideset slave
             */
            virtual ElementType
            slave_type() const ;

//------------------------------------------------------------------------------

            /**
             * set the pointer for the boundary condition
             */
            void
            set_boundary_condition( const BoundaryCondition * aBoundaryCondition );

//------------------------------------------------------------------------------

            /**
             * return the boundary condition
             */
            const BoundaryCondition *
            boundary_condition() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual uint
            number_of_thin_shell_layers() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual uint
            number_of_ghost_sidesets() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual const Material *
            thin_shell_material( const uint aLayerIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual real
            thin_shell_thickness( const uint aLayerIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual real
            thin_shell_thickness() const ;

//------------------------------------------------------------------------------

            /**
             * function telling if values form this group are computed
             * ( default : true )
             */
             bool
             is_active() const ;

//------------------------------------------------------------------------------

            /**
             * set the active flag of the group
             */
            void
            activate( bool aFlag );

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------


            void
            assume_isogeometry();

//------------------------------------------------------------------------------
// Material Helpers
//------------------------------------------------------------------------------

            // called by constructor
            void
            link_material_functions();

//------------------------------------------------------------------------------

            void
            linear_elasticity_PlaneStress(
                    Matrix< real > & aC,
                    const         real   aT ) const;

//------------------------------------------------------------------------------

            void
            linear_elasticity_3d(
                    Matrix< real > & aC,
                    const         real   aT ) const;

//------------------------------------------------------------------------------

            void
            thermal_conductivity_2d(
                    Matrix< real > & aLambda,
                    const         real   aT ) const;

//------------------------------------------------------------------------------

            void
            thermal_conductivity_3d(
                    Matrix< real > & aLambda,
                    const         real   aT ) const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline DofManagerBase *
        Group::parent()
        {
            return mParent;
        }

//------------------------------------------------------------------------------

        inline GroupType
        Group::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        inline id_t
        Group::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline index_t
        Group::number_of_elements() const
        {
            return mNumberOfElements;
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_nodes_per_element() const
        {
            return mNumberOfNodesPerElement ;
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_edges_per_element() const
        {
            return mesh::number_of_edges( this->element_type() );
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_faces_per_element() const
        {
            return mesh::number_of_faces( this->element_type() );
        }

//------------------------------------------------------------------------------

        inline Cell< fem::Element * > &
        Group::elements()
        {
            return mElements;
        }
//------------------------------------------------------------------------------

        inline const IntegrationData *
        Group::integration()
        {
            return mIntegrationData ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Group::integration_weights() const
        {
            return mIntegrationData->weights() ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::integration_points() const
        {
            return mIntegrationData->points() ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Group::n( const uint aIndex ) const
        {
            return mIntegrationData->phi( aIndex );
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::node_coords()
        {
            return mNodeCoords ;
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::N( const uint aIndex ) const
        {
            return mIntegrationData->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::Nvector( const uint aIndex ) const
        {
            return mIntegrationData->Nvector( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::dNdXi( const uint aIndex ) const
        {
            return mIntegrationData->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::d2NdXi2( const uint aIndex ) const
        {
            return mIntegrationData->d2NdXi2( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::G( const uint aIndex ) const
        {
            return mGeometryIntegrationData->N( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::dGdXi( const uint aIndex ) const
        {
            return mGeometryIntegrationData->dNdXi( aIndex );
        }

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        Group::d2GdXi2( const uint aIndex ) const
        {
            return mGeometryIntegrationData->d2NdXi2( aIndex );
        }

//------------------------------------------------------------------------------

        inline uint &
        Group::work_index()
        {
            return mWorkIndex;
        }


//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_H()
        {
            return mWorkH;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_E()
        {
            return mWorkE;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_dEdXi()
        {
            return mWorkDEDXi ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_dEdX()
        {
            return mWorkDEDX ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_M()
        {
            return mWorkM;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_G()
        {
            return mWorkG;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_J()
        {
            return mWorkJ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_invJ()
        {
            return mWorkinvJ;
        }

//------------------------------------------------------------------------------

        inline real &
        Group::work_det_J()
        {
            return mWorkDetJ ;
        }

//------------------------------------------------------------------------------

        inline real
        Group::dx() const
        {
            return mWorkDetJ ;
        }

//------------------------------------------------------------------------------

        inline real &
        Group::work_det_invJ()
        {
            return mWorkDetInvJ ;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_K()
        {
            return mWorkK;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_L()
        {
            return mWorkL;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_dNdX()
        {
            return mWorkDNDX;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_C()
        {
            return mWorkC;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_D()
        {
            return mWorkD;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_N()
        {
            return mWorkN;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_B()
        {
            return mWorkB;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_rhs()
        {
            return mWorkRhs;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_phi()
        {
            return mWorkphi;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Phi()
        {
            return mWorkPhi;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_psi()
        {
            return mWorkpsi;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Psi()
        {
            return mWorkPsi;
        }
//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_chi()
        {
            return mWorkchi;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Chi()
        {
            return mWorkChi;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_sigma()
        {
            return mWorksigma;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Sigma()
        {
            return mWorkSigma;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_tau()
        {
            return mWorktau;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Tau()
        {
            return mWorkTau;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_gamma()
        {
            return mWorkgamma;
        }
//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_theta()
        {
            return mWorktheta;
        }
//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Gamma()
        {
            return mWorkGamma;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_X()
        {
            return mWorkX;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Xm()
        {
            return mWorkXm;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::work_Xs()
        {
            return mWorkXs;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_nedelec()
        {
            return mWorkNedelec ;
        }

//------------------------------------------------------------------------------

        // get an element using its ID instead of index
        inline Element *
        Group::element( const id_t aID )
        {
            return mElementMap( aID );
        }

//------------------------------------------------------------------------------

        inline const Material *
        Group::material() const
        {
            return mMaterial;
        }

//------------------------------------------------------------------------------
// Material functions
//------------------------------------------------------------------------------

        inline void
        Group::linear_elasticity(
                Matrix< real > & aC,
                const         real aT ) const
        {
            ( this->*mLinearElasticity )( aC, aT );
        }

//------------------------------------------------------------------------------

        // flag telling field if an RHS side exists
        inline bool
        Group::has_rhs() const
        {
            return mHasRHS ;
        }

//------------------------------------------------------------------------------

        // set or unset the rhs flag
        inline void
        Group::set_rhs_flag( const bool aFlag )
        {
            mHasRHS = aFlag ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::thermal_conductivity(
                Matrix< real > & aLambda,
                const         real aT ) const
        {
            ( this->*mThermalConductivity )( aLambda, aT );
        }

//------------------------------------------------------------------------------

        inline void
        Group::linear_elasticity_PlaneStress(
                Matrix< real > & aC, const real aT ) const
        {
            mMaterial->C_ps( aC, aT );
        }

//------------------------------------------------------------------------------

        inline void
        Group::linear_elasticity_3d(
                Matrix< real > & aC, const real aT ) const
        {
            mMaterial->C( aC, aT );
        }

//------------------------------------------------------------------------------

        inline void
        Group::thermal_conductivity_2d(
                Matrix< real > & aLambda, const real aT ) const
        {
            mMaterial->lambda_p( aLambda, aT );
        }

//------------------------------------------------------------------------------

        inline void
        Group::thermal_conductivity_3d(
                Matrix< real > & aLambda, const real aT ) const
        {
            mMaterial->lambda_3d( aLambda, aT );
        }

//------------------------------------------------------------------------------

        inline const id_t &
        Group::mesh_id() const
        {
            return mMeshID ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::set_mesh_id( const id_t & aID )
        {
            mMeshID = aID ;
        }

//------------------------------------------------------------------------------

        inline const InterpolationFunction *
        Group::interpolation_function() const
        {
            return mIntegrationData->function() ;
        }

//------------------------------------------------------------------------------


        inline void
        Group::set_domain_type( const DomainType aType )
        {
            mDomainType = aType ;
        }

//------------------------------------------------------------------------------

        inline DomainType
        Group::domain_type() const
        {
            return mDomainType ;
        }

//------------------------------------------------------------------------------

        inline bool
        Group::is_active() const
        {
            return mIsActive ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::activate( bool aFlag )
        {
            mIsActive = aFlag ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group:: work_geo()
        {
            return mWorkgeo ;
        }

//------------------------------------------------------------------------------

        inline Vector< real > &
        Group::work_normal()
        {
            return mWorkNormal ;
        }

//------------------------------------------------------------------------------

        inline Calculator *
        Group::calculator()
        {
            return mCalc ;
        }



//------------------------------------------------------------------------------

        inline ElementType
        Group::element_type() const
        {
            return mElementType;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_GROUP_HPP
