//
// Created by gregorygiard on 10/24/25.
//
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

        std::shared_ptr< Kernel >
        ThermalFactory::create_thermal_kernel()
        {
            BELFEM_ERROR( mMesh != nullptr, "No mesh defined" );

            Vector< id_t > tBlocks ;

            if ( mCommRank == 0 )
            {
                // count how many blocks are not air
                index_t tCount = 0 ;
                for ( mesh::Block * tBlock : mMesh->blocks() )
                {
                    switch ( tBlock->domain_type() )
                    {
                    case DomainType::Conductor :
                    case DomainType::ThinShell :
                    case DomainType::Ferro :
                    case DomainType::Buffer :
                    {
                        ++tCount ;
                        break ;
                    }
                    default:
                    {
                        break ;
                    }
                    }
                }
                tBlocks.set_size( tCount );
                tCount = 0 ;
                for ( mesh::Block * tBlock : mMesh->blocks() )
                {
                    switch ( tBlock->domain_type() )
                    {
                        case DomainType::Conductor :
                        case DomainType::ThinShell :
                        case DomainType::Ferro :
                        case DomainType::Buffer :

                        {
                            tBlocks( tCount++ ) = tBlock->id() ;
                            break ;
                        }
                        default:
                        {
                            break ;
                        }
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


            // todo: be aware that axisymmetry could also be a case
            mDimensionality = mMesh->number_of_dimensions() == 2
                ? ModelDimensionality::TwoD
                : ModelDimensionality::ThreeD;

            //Create the equation
            mThermalEquation =  new IWG_MaxwellThermal( mDimensionality, IwgType::MaxwellThermal, IwgMode::Iterative );

            mThermalEquation->select_blocks( tBlocks );

            comm_barrier();

            //Create the field
            mThermalField = mThermalKernel->create_field( mThermalEquation );
            mOwnThermalEquation = false ;

            // get the solver section from the input file
            const input::Section *tSolverSection = mInputFile->section(
                "solver" );

            // when we have thermal, we will be able to chose different settings for the thermal and the magnetic part
            const input::Section *tLinearSection = tSolverSection->
                section_exists( "linear thermal" )
                ? tSolverSection->section( "linear thermal" )
                : tSolverSection->section( "linear" );

            this->configure_solver( tLinearSection, mThermalField );

            for( Block * tGroup : mThermalKernel->dofmgr()->blocks() )
            {
                switch( mMesh->block( tGroup->id() )->domain_type() )
                {
                    case DomainType::Conductor :
                    case DomainType::ThinShell :
                    case DomainType::Air :
                    case DomainType::Ferro :
                    case DomainType::Buffer :
                    {
                        tGroup->activate( true );
                        tGroup->set_domain_type( mMesh->block( tGroup->id() )->domain_type() );
                        break ;
                    }
                    default:
                    {
                        tGroup->activate( false );
                    }
                }
            }

            // get enriched order
            uint tOrder = mUseEnrichment ? ( mMesh->max_element_order() == 1 ? 5 : 13 ) : 0 ;
            for( SideSet * tGroup : mThermalKernel->dofmgr()->sidesets())
            {
                switch ( mMesh->sideset( tGroup->id() )->domain_type() )
                {
                    case( DomainType::Dirichlet ) :
                    case( DomainType::Neumann ) :
                    {
                        tGroup->calculator()->set_integration_order( tOrder ) ;
                        tGroup->activate( true );
                        break;
                    }
                    default:
                    {
                        tGroup->activate( false );
                    }
                }
            }

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
