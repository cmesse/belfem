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

#include "globals.hpp"
#include "cl_ThermalFactory.hpp"
#include "cl_FEM_Controller.hpp"
#include "cl_IWG_MaxwellThermal.hpp"
#include "fn_check_unit.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {

        ThermalFactory::ThermalFactory( const string & aInputFile, Mesh * aMesh ) :
            mCommRank( comm_rank() ),
            mInputFile( new InputFile( aInputFile ) ),
            mBoundaryConditionFactory( this->create_bc_factory() ),
            mMesh( aMesh )
        {
            // create a parameter object
            mKernelParameters = new KernelParameters( aMesh );
        }

        ThermalFactory::ThermalFactory( const string & aInputFile, Kernel * aMagneticKernel ) :
            mCommRank( comm_rank() ),
            mInputFile( new InputFile( aInputFile ) ),
            mBoundaryConditionFactory( this->create_bc_factory() ),
            mMesh( aMagneticKernel->mesh() )

        {
            // create a parameter object
            mKernelParameters = new KernelParameters( aMagneticKernel );
        }

        ThermalFactory::~ThermalFactory()
        {
            delete mInputFile;

            // the factory object is ours; its boundary conditions belong to the
            // kernel once create_thermal_kernel() handed them over
            delete mBoundaryConditionFactory;

            for ( Domain *tDomain : mDomains )
            {
                delete tDomain;
            }

            if ( mOwnKernelParameters && mKernelParameters != nullptr )
            {
                delete mKernelParameters;
            }
            if ( mOwnThermalEquation && mThermalEquation != nullptr )
            {
                delete mThermalEquation;
            }
        }

//------------------------------------------------------------------------------

        bool
        ThermalFactory::block_is_supported( const DomainType aType )
        {
            switch ( aType )
            {
                case DomainType::Conductor :
                case DomainType::ThinShell :
                case DomainType::Ferro :
                case DomainType::Buffer :
                {
                    return true ;
                }
                default :
                {
                    return false ;
                }
            }
        }

//------------------------------------------------------------------------------

        std::shared_ptr< Kernel >
        ThermalFactory::create_thermal_kernel()
        {
            BELFEM_ERROR( mMesh != nullptr, "No mesh defined" );

            Vector< id_t > tBlocks ;

            if ( mCommRank == 0 )
            {
                index_t tCount = 0 ;
                for ( mesh::Block * tBlock : mMesh->blocks() )
                {
                    if ( ThermalFactory::block_is_supported( tBlock->domain_type() ) )
                    {
                        ++tCount ;
                    }
                }
                tBlocks.set_size( tCount );
                tCount = 0 ;
                for ( mesh::Block * tBlock : mMesh->blocks() )
                {
                    if ( ThermalFactory::block_is_supported( tBlock->domain_type() ) )
                    {
                        tBlocks( tCount++ ) = tBlock->id() ;
                    }
                }
                comm_barrier() ;
                share( tBlocks );
            }
            else
            {
                comm_barrier() ;
                Cell< id_t > tAllBlocks ;
                receive( tAllBlocks );
                uint tCount = 0 ;

                for ( id_t tID : tAllBlocks )
                {
                    if ( mMesh->block_exists( tID ) )
                    {
                        ++tCount ;
                    }
                }
                tBlocks.set_size( tCount );
                tCount = 0 ;
                for ( id_t tID : tAllBlocks )
                {
                    if ( mMesh->block_exists( tID ) )
                    {
                        tBlocks( tCount++ ) = tID ;
                    }
                }
            }

            mKernelParameters->select_blocks( tBlocks );

            // construct in place: the temporary-plus-copy form invoked the
            // implicit copy constructor ( Kernel declares a destructor, so
            // C++17 suppresses the move ) and copied uninitialized members
            mThermalKernel = std::make_shared<Kernel>( mKernelParameters );

            mThermalKernel->claim_parameter_ownership( true );
            mOwnKernelParameters = false ;


            // TODO(cm): only TwoD and ThreeD are selected here; AxSymmX and
            // AxSymmY need a deck key
            mDimensionality = mMesh->number_of_dimensions() == 2
                ? ModelDimensionality::TwoD
                : ModelDimensionality::ThreeD;

            //Create the equation
            mThermalEquation =  new IWG_MaxwellThermal( mDimensionality, IwgType::MaxwellThermal, IwgMode::Iterative );

            mThermalEquation->select_blocks( tBlocks );

            comm_barrier();

            // the kernel owns the equation from here on; create_field() only
            // links it
            mThermalKernel->add_equation( mThermalEquation );
            mOwnThermalEquation = false ;
            mThermalField = mThermalKernel->create_field( mThermalEquation );

            // get the solver section from the input file
            const input::Section *tSolverSection = mInputFile->section(
                "solver" );

            // when we have thermal, we will be able to chose different settings for the thermal and the magnetic part
            const input::Section *tLinearSection = tSolverSection->
                section_exists( "linear thermal" )
                ? tSolverSection->section( "linear thermal" )
                : tSolverSection->section( "linear" );

            this->configure_solver( tLinearSection, mThermalField );

            // the dof manager holds the selected blocks only, so the
            // refusal cannot fire unless the selection and this dispatch
            // stop reading the same predicate
            for( Block * tGroup : mThermalKernel->dofmgr()->blocks() )
            {
                const DomainType tType = mMesh->block( tGroup->id() )->domain_type() ;

                BELFEM_ERROR( ThermalFactory::block_is_supported( tType ),
                    "block %lu reached the thermal kernel with domain type %s, which is not a thermal block type: the selection in create_thermal_kernel and this dispatch disagree",
                    ( long unsigned int ) tGroup->id(),
                    to_string( tType ).c_str() );

                tGroup->activate( true );
                tGroup->set_domain_type( tType );
            }

            // the thermal equation selects no sidesets, so the dof manager
            // holds none; a thermal boundary condition on a sideset is not
            // supported

            //Send the dofmngr to the boundary conditions
            mBoundaryConditionFactory->set_fields( mThermalField );

            //Send boundary conditions to the kernel
            for ( PhysicalBoundaryCondition *tBC : mBoundaryConditionFactory->
                  boundary_conditions() )
            {
                mThermalKernel->add_boundary_condition( tBC );
            }

            //Initialize boundary conditions
            for ( PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions() )
            {
                if ( tBC->type() == BoundaryConditionType::Bearing )
                {
                    for (id_t tID : tBC->domains())
                    {
                        mThermalKernel->dofmgr()->bearing( tID )->impose_dirichlet( 0.0 );
                    }
                }
                else if ( tBC->type() == BoundaryConditionType::Dirichlet )
                {
                    for (id_t tID: tBC->domains())
                    {
                        mThermalKernel->dofmgr()->sideset( tID )->impose_dirichlet( 0.0 ) ;
                    }
                }
            }

            mThermalField->initialize() ;

            BELFEM_ERROR( !std::isnan( gTbulk ), "Warning : Initial temperature was not set" ) ;

            // initialize temperature field
            mThermalKernel->mesh()->field_data("T").fill( gTbulk );
            mThermalKernel->dofmgr()->field_data( "T" ).fill( gTbulk ) ;

            // the initialize() above copied the field into the dofs while
            // the field was still zero; re-seed the free dofs now that T
            // is set, otherwise the first Picard residual sees x = 0 and
            // reports exactly 1.0
            mThermalField->seed_dof_values() ;

            // one mesh global per labeled thermal boundary condition, on the
            // SAME mesh the magnetic kernel writes ( they alias in the
            // coupled run ). RANK 0 ONLY — cf. the Maxwell twin
            if ( mCommRank == 0 )
            {
                for ( PhysicalBoundaryCondition * tBC :
                      mBoundaryConditionFactory->boundary_conditions() )
                {
                    if ( tBC->label().size() > 0 )
                    {
                        Mesh * tMesh = mThermalKernel->dofmgr()->mesh() ;

                        BELFEM_ERROR( ! tMesh->global_variable_exists( tBC->label() ),
                            "duplicate boundary condition global '%s' - label the sections to disambiguate",
                            tBC->label().c_str() );

                        tMesh->create_global_variable( tBC->label(), tBC->value() );
                    }
                }
            }

            return mThermalKernel ;
        }

//------------------------------------------------------------------------------

        void
        ThermalFactory::configure_solver( const input::Section *aSection,
                                          DofManager *aField )
        {
            // read solver settings
            SolverParameters tParams( aSection );

            // get the solver from the kernel
            aField->set_solver( tParams );

        }

//------------------------------------------------------------------------------

        ThermalBoundaryConditionFactory *
        ThermalFactory::create_bc_factory()
        {
            //Create the Thermal BC Factory
            BELFEM_ERROR(mInputFile->section_exists("boundary conditions" ),"No boundary conditions defined") ;
            if ( mInputFile->section( "boundary conditions" )->section_exists( "thermal" ) )
            {
                return new ThermalBoundaryConditionFactory( mInputFile->section( "boundary conditions" )->section( "thermal" ) ) ;
            }
            else
            {
                // if no thermal section exists, just initialize an empty BC Factory
                return new ThermalBoundaryConditionFactory();
            }

        }

//------------------------------------------------------------------------------

    }
}
