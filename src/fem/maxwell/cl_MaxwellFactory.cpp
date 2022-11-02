//
// Created by christian on 9/20/21.
//


#include "commtools.hpp"
#include "cl_MaxwellFactory.hpp"
#include "units.hpp"
#include "en_FEM_DomainType.hpp"
#include "fn_unique.hpp"
#include "fn_check_unit.hpp"
#include "cl_Block.hpp"
#include "cl_Mesh_Scissors.hpp"
#include "cl_Mesh_TapeRoller.hpp"
#include "cl_MaterialFactory.hpp"

#include "en_SolverEnums.hpp"
#include "en_FEM_MagfieldBcType.hpp"
#include "cl_ODE_Integrator.hpp"
#include "cl_BhExtrapolationOde.hpp"
#include "en_ODE_Type.hpp"

//#include "cl_IWG_Maxwell_HA_Tri3.hpp"
//#include "cl_IWG_Maxwell_HA_Tri6.hpp"

#include "cl_IWG_Maxwell_HPhi_Tri3.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri6.hpp"
//#include "cl_IWG_Maxwell_HPhi_Tet4.hpp"
//#include "cl_IWG_Maxwell_HPhi_Tet10.hpp"

#include "cl_IWG_Maxwell_L2_Current.hpp"
#include "cl_IWG_Maxwell_L2_Magfield.hpp"

#include "cl_OneDMapper.hpp"
#include "fn_linspace.hpp"
#include "cl_Spline.hpp"
#include "fn_create_beam_poly.hpp"
#include "fn_polyfit.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"

#include "fn_Mesh_compute_surface_normals.hpp"
#include "cl_Element_Factory.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        MaxwellFactory::MaxwellFactory( const InputFile & aInputFile )  :
            IwgFactory( this->load_mesh( aInputFile ) ),
            mInputFile( aInputFile ),
            mParameters( new KernelParameters( this->mesh() ) )
        {
            // create communication table
            this->create_commtable();

            this->create_domain_groups();
            this->select_blocks();

            this->read_formulation();

#ifdef BELFEM_FERRO_HPHIA
#ifdef BELFEM_FERRO_LINEAR
            if( mFormulation == IwgType::MAXWELL_HPHI_TRI6 )
            {
                if( comm_rank() == 0 )
                {
                    this->reduce_order_of_ferro_blocks_tri6();
                }
                comm_barrier() ;
            }
#endif
#endif
            // create the edges for this mesh
            mMesh->create_edges( false, mNedelecBlocks, mNedelecSideSets ) ;


            // check if this is a higher order mesh
            if( mMesh->max_element_order() > 1 )
            {
                // create faces
                mMesh->create_faces( false, mNedelecBlocks, mNedelecSideSets ) ;
            }


            // read the timestepping information
            this->read_timestep() ;

            this->create_bcs() ;

            // define which blocks and sidesets are used and what they
            // represent
            this->create_topology();

            // Dirichlet bcs do not need elements. We remove those
            this->blacklist_dirichlet_sidesets();

            this->create_kernel();

            // build the iwg and the dof manager for the magnetic fields
            this->create_magnetic_field();

            // add postprocessors to magnetic field
            if( mHaveJ )
            {
                this->create_current_projector();
            }
            if( mHaveB )
            {
                this->create_magfield_projector() ;
            }
        }

//------------------------------------------------------------------------------

        MaxwellFactory::~MaxwellFactory()
        {
            if( mOwnKernel )
            {
                delete mKernel ;
            }

            if( mOwnBoundaryConditions )
            {
                for ( BoundaryCondition * tBC : mBoundaryConditions )
                {
                    delete tBC ;
                }
            }

            for( Material * tMaterial : mMaterials )
            {
                // check if material is owned by the kernel
                if( ! tMaterial->is_flagged() )
                {
                    delete tMaterial ;
                }
            }


            for( DomainGroup * tGroup : mBlocks )
            {
                delete tGroup ;
            }
            mBlockMap.clear() ;


            for( DomainGroup * tGroup : mBoundaries )
            {
                delete tGroup ;
            }

            mBoundaryMap.clear() ;

            for( DomainGroup * tGroup : mRawSymmetries )
            {
                delete tGroup ;
            }
            for( DomainGroup * tGroup : mRawAntiSymmetries )
            {
                delete tGroup ;
            }
            for( DomainGroup * tGroup : mSymmetries )
            {
                delete tGroup ;
            }
            for( DomainGroup * tGroup : mAntiSymmetries )
            {
                delete tGroup ;
            }

            for( DomainGroup * tGroup : mInterfaces )
            {
                delete tGroup ;
            }

            for( DomainGroup * tGroup : mCuts )
            {
                delete tGroup ;
            }
            mCutMap.clear() ;

            for( DomainGroup * tGroup : mTapes )
            {
                delete tGroup ;
            }

            if ( mOwnMesh )
            {
                delete mMesh ;
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_commtable()
        {
            // create the comm table ( we don't have a kernel yet )
            proc_t tNumProcs = comm_size() ;

            if( comm_rank() == 0 && tNumProcs > 1 )
            {
                mCommTable.set_size( tNumProcs-1 );
                uint tCount = 0 ;
                for( proc_t p=1; p<tNumProcs; ++p )
                {
                    mCommTable( tCount++ ) = p ;
                }
            }
        }

//------------------------------------------------------------------------------

        Mesh *
        MaxwellFactory::load_mesh( const InputFile & aInputFile )
        {
            Mesh * aMesh = new Mesh( aInputFile.section( "mesh" )->get_string( "meshFile" ) );

            if( gComm.rank() == 0 )
            {

                // get mesh unit
                string tUnit = aInputFile.section( "mesh" )->get_string( "meshUnit" );

                // we only need to scale the mesh if it is not in meters
                if ( tUnit != "m" )
                {
                    // grab value
                    value tValue = unit_to_si( tUnit );

                    BELFEM_ERROR( check_unit( tValue, "m" ), "Invalid length unit: %s", tUnit.c_str());

                    aMesh->scale_mesh( tValue.first );
                }
            }

            return aMesh ;
        }
//------------------------------------------------------------------------------

        string
        MaxwellFactory::outfile() const
        {
            return mInputFile.section( "output")->get_string("outputFile");
        }

//------------------------------------------------------------------------------

        string
        MaxwellFactory::backupfile() const
        {
            return mInputFile.section("output")->key_exists("backupfile") ?
                   mInputFile.section("output")->get_string("backupfile") :
                   "memdump.hdf5" ;
        }

//------------------------------------------------------------------------------

        bool
        MaxwellFactory::restart() const
        {
            return mInputFile.section("timestepping")->key_exists("restart") ?
                   mInputFile.section("timestepping")->get_bool( "restart") :
                   false ;

        }

//------------------------------------------------------------------------------

        real
        MaxwellFactory::maxtime() const
        {
            return mInputFile.section("timestepping")->get_value("maxtime", "s").first ;
        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::read_formulation()
        {

            // only the master proc knows the element order,
            // we need to synchronize that information
            uint tElementOrder = mMesh->max_element_order() ;
            comm_barrier() ;
            broadcast( 0, tElementOrder );

            // read the tau parameter
            if( mInputFile.section( "maxwell")->section("solver")->section("nonlinear")->key_exists("stabilize") )
            {
                mTau = mInputFile.section("maxwell")->section("solver")->section("nonlinear")->get_real( "stabilize");
            }

            // get string from input file
            string tKey = string_to_lower(
                    mInputFile.section("maxwell")->get_string("formulation") );

            // get the number of dimensions from the mesh
            uint tNumDim = mMesh->number_of_dimensions() ;

            if( tKey == "h-a" )
            {
                if( tNumDim == 2 )
                {
                    switch ( tElementOrder )
                    {
                        case( 1 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HA_TRI3 ;
                            mProjectJ    = IwgType::MAXWELL_L2J_TRI3 ;
                            mProjectB    = IwgType::MAXWELL_HA_L2B_TRI3 ;
                            break ;
                        }
                        case( 2 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HA_TRI6 ;
                            mProjectJ    = IwgType::MAXWELL_L2J_TRI6 ;
                            mProjectB    = IwgType::MAXWELL_HA_L2B_TRI6 ;
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Mesh must be first or second order" );
                        }
                    }
                }
                else if ( tNumDim == 3 )
                {
                    switch ( tElementOrder )
                    {
                        case( 1 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HA_TET4 ;
                            mProjectJ    = IwgType::MAXWELL_L2J_TET4 ;
                            mProjectB    = IwgType::MAXWELL_HA_L2B_TET4 ;
                            break ;
                        }
                        case( 2 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HA_TET10 ;
                            mProjectJ    = IwgType::MAXWELL_L2J_TET10 ;
                            mProjectB    = IwgType::MAXWELL_HA_L2B_TET10 ;
                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Mesh must be first or second order" );
                        }
                    }
                }
                else
                {
                    BELFEM_ERROR( false, "Invalid mesh dimension %u", ( unsigned int ) tNumDim );
                }
            }
            else if ( tKey == "h-phi" || tKey == "h-φ" || tKey == "h-Φ" )
            {
                if( tNumDim == 2 )
                {
                    switch ( tElementOrder )
                    {
                        case ( 1 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HPHI_TRI3;
                            mProjectJ    = IwgType::MAXWELL_L2J_TRI3 ;
                            mProjectB    = IwgType::MAXWELL_HPHI_L2B_TRI3 ;
                            break;
                        }
                        case ( 2 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HPHI_TRI6;
                            mProjectJ    = IwgType::MAXWELL_L2J_TRI6 ;
                            mProjectB    = IwgType::MAXWELL_HPHI_L2B_TRI6 ;
                            break;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Mesh must be first or second order" );
                        }
                    }
                }
                else if ( tNumDim == 3 )
                {
                    switch ( tElementOrder )
                    {
                        case ( 1 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HPHI_TET4;
                            mProjectJ    = IwgType::MAXWELL_L2J_TET4 ;
                            mProjectB    = IwgType::MAXWELL_HPHI_L2B_TET4 ;
                            break;
                        }
                        case ( 2 ) :
                        {
                            mFormulation = IwgType::MAXWELL_HPHI_TET10;
                            mProjectJ    = IwgType::MAXWELL_L2J_TET10 ;
                            mProjectB    = IwgType::MAXWELL_HPHI_L2B_TET10 ;
                            break;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Mesh must be first or second order" );
                        }
                    }
                }
                else
                {
                    BELFEM_ERROR( false, "Invalid mesh dimension %u", ( unsigned int ) tNumDim );
                }
            }
            else
            {
                BELFEM_ERROR( false, "Unknown formulation : %s", tKey.c_str() );
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::read_timestep()
        {
            // get the section in the input file
            const input::Section * tSection = mInputFile.section("timestepping");

            // read the timestep
            mTimeStep = tSection->get_value( "timestep" , "s" ).first ;

            if( tSection->key_exists("theta") )
            {
                // reat the theta parameter
                string tString = string_to_lower( tSection->get_string( "theta" ));

                if ( tString == "galerkin" )
                {
                    mTheta = 2.0 / 3.0;
                }
                else if ( tString == "crank–nicolson" || tString == "cranknicolson" )
                {
                    mTheta = 0.5;
                }
                else if ( tString == "forward-euler" ||
                          tString == "forwardeuler" ||
                          tString == "euler-forward" ||
                          tString == "eulerforward" ||
                          tString == "explicit-euler" ||
                          tString == "expliciteuler" ||
                          tString == "euler-explicit" ||
                          tString == "eulerexplicit" )
                {
                    mTheta = 0.0;
                }
                else if ( tString == "backward-euler" ||
                          tString == "backwardeuler" ||
                          tString == "euler-backward" ||
                          tString == "eulerbackward" ||
                          tString == "implicit-euler" ||
                          tString == "impliciteuler" ||
                          tString == "euler-implicit" ||
                          tString == "eulerimplicit" )
                {
                    mTheta = 1.0;
                }
                else
                {

                    mTheta = to_real( tString );
                    BELFEM_ERROR( !std::isnan( mTheta ), "unknown theta: %s", tString.c_str());

                    // fix for galerkin
                    if ( mTheta == 0.66 )
                    {
                        mTheta = 2.0 / 3.0;
                    }
                }
            }
            else
            {
                // use backward euler as default
                mTheta = 1.0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::read_field_switches()
        {
            const input::Section * tSection = mInputFile.section("maxwell");

            if( tSection->key_exists( "output") )
            {
                // get entry
                string tString = tSection->get_string("output") ;
                tString = search_and_replace( tString,  ",", " " );
                tString = search_and_replace( tString,  ":", " " );
                tString = search_and_replace( tString,  "[", " " );
                tString = search_and_replace( tString,  "]", " " );
                tString = search_and_replace( tString,  "=", " " );

                Cell< string > tWords = string_to_words( clean_string( tString ) );

                // help flags for detecting smooth factor
                int tFlag = 0 ;

                for( string & tWord : tWords )
                {
                    if( tWord == "j" )
                    {
                        mHaveJ = true ;
                        tFlag = 1 ;
                    }
                    else if( tWord == "b" )
                    {
                        mHaveB = true ;
                        tFlag = 3 ;
                    }
                    else if( tWord == "smooth" ||  tWord == "s" )
                    {
                        ++tFlag ;
                    }
                    else if( tFlag == 2 )
                    {
                        mAlphaJ = to_real( tWord );
                        tFlag = 0 ;
                    }
                    else if ( tFlag == 4 )
                    {
                        mAlphaB = to_real( tWord );
                        tFlag = 0 ;

                    }
                    else if( tWord == "normb" || tWord == "norm_b" )
                    {
                        mComputeNormB = true ;
                        tFlag = 0 ;
                    }
                }
            }

            mHaveThermal = mInputFile.section_exists("thermal");
            mHaveMech    = mInputFile.section_exists("mech");

            // check if linear flag is set
            /*if( tSection->key_exists( "linear") )
            {
                mEnforceLinear = tSection->get_bool("linear") ;
            }

            // make parallel
            if( comm_rank() == 0 )
            {
                mEnforceLinear = mEnforceLinear && mMesh->max_element_order() > 1 ;
            }
            comm_barrier() ;

            uint tFlag = mEnforceLinear ? 1 : 0 ;
            broadcast( 0, tFlag );
            mEnforceLinear = tFlag == 1 ? 1 : 0 ; */
        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_materials()
        {
            // dandy protection
            for( Material * tMaterial : mMaterials )
            {
                delete tMaterial ;
            }
            mMaterialMap.clear() ;

            if( mInputFile.section_exists( "maxwell") )
            {
                // create the material factory
                MaterialFactory tFactory ;

                if( mInputFile.section("maxwell")->section_exists("materials") )
                {

                    // get materials section
                    const input::Section * tMaterials
                            = mInputFile.section( "maxwell" )->section( "materials" );

                    index_t tNumSections = tMaterials->num_sections();

                    for ( index_t k = 0; k < tNumSections; ++k )
                    {
                        // get material entry
                        const input::Section * tInput = tMaterials->section( k );


                        // create the material
                        Material * tMaterial = this->create_maxwell_material( tInput ) ;

                        // add material to list
                        mMaterials.push( tMaterial );

                        // add material to label map
                        string tLower = string_to_lower( tMaterial->label() ) ;
                        mMaterialLabelMap[ tMaterial->label()  ] = tMaterial ;

                        if( tMaterial->type() == MaterialType::Maxwell )
                        {
                            // get labels
                            const Cell< string > & tTargets
                                = reinterpret_cast< MaxwellMaterial * >( tMaterial )->domain_labels();

                            // add material to map
                            for ( index_t t = 0; t < tTargets.size(); ++t )
                            {
                                mMaterialMap[ tTargets( t ) ] = tMaterial;
                            }
                        }
                    }
                }

                // check if tape materials are defined
                if( mTapeMaterialLabels.size() > 0 )
                {
                    uint tCount = 0 ;

                    // allocate container for tape materials
                    mTapeMaterials.set_size( mTapeMaterialLabels.size(), nullptr );

                    for( const string & tLabel : mTapeMaterialLabels )
                    {
                        if ( ! mMaterialLabelMap.key_exists( tLabel ) )
                        {
                            // we check if the material is intrinsic, if not, we assume that it is a maxwell material
                            MaterialType tMaterialType = string_to_material_type( tLabel, false );

                            // create the material
                            Material * tMaterial = tFactory.create_material( tMaterialType );

                            // add material to container
                            mMaterials.push( tMaterial );

                            // add material to tape container
                            mTapeMaterials( tCount++ ) = tMaterial ;
                        }
                        else
                        {
                            // grab material from list and add to table
                            mTapeMaterials( tCount++ ) = mMaterialLabelMap( tLabel );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        Material *
        MaxwellFactory::create_maxwell_material( const input::Section * aInput )
        {
            if( aInput->key_exists( "intrinsic" ) )
            {
                // create a material factory
                MaterialFactory tFactory ;

                // return the new material from the database
                Material * aMaterial = tFactory.create_material( aInput->label() );

                if( aInput->key_exists( "rrr" ) )
                {
                    aMaterial->set_rrr( aInput->get_real( "rrr" ) );
                }

                return aMaterial ;
            }
            else
            {
                // create new material
                MaxwellMaterial * aMaterial = new MaxwellMaterial( aInput->label() );

                // legacy patch, resistivity sounds better than resistance
                string tType = "" ;
                if ( aInput->key_exists( "resistance" ) )
                {
                    tType = string_to_lower( aInput->get_string( "resistance" ));
                }
                else if ( aInput->key_exists( "resistivity" ) )
                {
                    tType = string_to_lower( aInput->get_string( "resistivity" ));
                }

                // check if resistivity exists
                if ( tType.size() > 0 )
                {

                    // sanity checks
                    if ( tType == "powerlaw-ej" || tType == "powerlaw-ejb" || tType == "powerlaw-ejt" ||
                         tType == "powerlaw-ejbt" )
                    {

                        BELFEM_ERROR(
                                ( aInput->key_exists( "ec" ) && !aInput->key_exists( "rho_c" ))
                                || ( !aInput->key_exists( "ec" ) && aInput->key_exists( "rho_c" )),
                                "must prescribe either ec or rho_c, not both" );

                        BELFEM_ERROR(
                                ( aInput->key_exists( "n" ) && !aInput->key_exists( "n0" ))
                                || ( !aInput->key_exists( "n" ) && aInput->key_exists( "n0" )),
                                "variable must be called either n or n0, not both" );

                        // don't throw this error if we provide a database
                        if( !  ( aInput->key_exists("database") && aInput->key_exists("dataset") ) )
                        {
                            BELFEM_ERROR(
                                    ( aInput->key_exists( "jc" ) && !aInput->key_exists( "jc0" ))
                                    || ( !aInput->key_exists( "jc" ) && aInput->key_exists( "jc0" )),
                                    "must provide either jc, jc0 or database and dataset" );
                        }
                    }



                    if ( tType == "constant" )
                    {

                        if ( aInput->key_exists( "rho_c" ))
                        {
                            aMaterial->set_rho_el_const( aInput->get_value( "rho_c", "Ohm*m" ).first );
                        }
                        if ( aInput->key_exists( "rho" ))
                        {
                            aMaterial->set_rho_el_const( aInput->get_value( "rho", "Ohm*m" ).first );
                        }
                        else
                        {
                            real tEc = aInput->get_value( "ec", "V/m" ).first;

                            real tJc = aInput->key_exists( "jc" ) ?
                                  aInput->get_value( "jc", "A/m^2" ).first :
                                  aInput->get_value( "jc0", "A/m^2" ).first;

                            aMaterial->set_rho_el_const( tEc / tJc );
                        }
                    }
                    else if ( tType == "powerlaw-ej" || tType == "powerlaw" )
                    {
                        real tJc = BELFEM_QUIET_NAN;

                        string tDatabase;
                        string tDataset;

                        if ( aInput->key_exists( "database" ) && aInput->key_exists( "dataset" ))
                        {
                            tDatabase = aInput->get_string( "database" );
                            tDataset = aInput->get_string( "dataset" );
                            tJc = 0;
                        }
                        else
                        {
                            tJc = aInput->key_exists( "jc0" ) ?
                                  aInput->get_value( "jc0", "A/m^2" ).first :
                                  aInput->get_value( "jc", "A/m^2" ).first;
                        }

                        real tN = aInput->key_exists( "n0" ) ?
                                  aInput->get_value( "n0", "-" ).first :
                                  aInput->get_value( "n", "-" ).first;

                        real tEc = aInput->key_exists( "rho_c" ) ?
                                   aInput->get_value( "rho_c", "Ohm*m" ).first * tJc :
                                   aInput->get_value( "ec", "V/m" ).first;

                        aMaterial->set_rho_el_ej( tEc, tJc, tN );

                        if( tJc == 0 )
                        {
                            aMaterial->set_j_crit( tDatabase, tDataset );
                        }

                    }
                    else if ( tType == "powerlaw-ejt" )
                    {

                        real tJc = aInput->key_exists( "jc0" ) ?
                                   aInput->get_value( "jc0", "A/m^2" ).first :
                                   aInput->get_value( "jc", "A/m^2" ).first;

                        real tN0 = aInput->key_exists( "n0" ) ?
                                   aInput->get_value( "n0", "-" ).first :
                                   aInput->get_value( "n", "-" ).first;

                        real tEc = aInput->key_exists( "rho_c" ) ?
                                   aInput->get_value( "rho_c", "Ohm*m" ).first * tJc :
                                   aInput->get_value( "ec", "V/m" ).first;


                        real tT0 = aInput->get_value( "T0", "K" ).first;
                        real tTc = aInput->get_value( "Tc", "K" ).first;

                        aMaterial->set_rho_el_ejt( tEc, tJc, tN0, tT0, tTc );
                    }
                    else if ( tType == "powerlaw-ejb" )
                    {
                        real tJc = aInput->key_exists( "jc0" ) ?
                                   aInput->get_value( "jc0", "A/m^2" ).first :
                                   aInput->get_value( "jc", "A/m^2" ).first;

                        real tEc = aInput->key_exists( "rho_c" ) ?
                                   aInput->get_value( "rho_c", "Ohm*m" ).first * tJc :
                                   aInput->get_value( "ec", "V/m" ).first;

                        real tN0 = aInput->key_exists( "n0" ) ?
                                   aInput->get_value( "n0", "-" ).first :
                                   aInput->get_value( "n", "-" ).first;

                        real tN1 = aInput->get_value( "n1", "-" ).first;

                        real tB0 = aInput->get_value( "b0", "T" ).first;

                        aMaterial->set_rho_el_ejb( tEc, tJc, tN0, tN1, tB0 );
                    }
                    else if ( tType == "powerlaw-ejbt" )
                    {
                        real tJc = aInput->key_exists( "jc0" ) ?
                                   aInput->get_value( "jc0", "A/m^2" ).first :
                                   aInput->get_value( "jc", "A/m^2" ).first;

                        real tEc = aInput->key_exists( "rho_c" ) ?
                                   aInput->get_value( "rho_c", "Ohm*m" ).first * tJc :
                                   aInput->get_value( "ec", "V/m" ).first;

                        real tN0 = aInput->key_exists( "n0" ) ?
                                   aInput->get_value( "n0", "-" ).first :
                                   aInput->get_value( "n", "-" ).first;

                        real tN1 = aInput->get_value( "n1", "-" ).first;

                        real tB0 = aInput->get_value( "b0", "T" ).first;

                        real tT0 = aInput->get_value( "T0", "K" ).first;
                        real tTc = aInput->get_value( "Tc", "K" ).first;

                        aMaterial->set_rho_el_ejbt( tEc, tJc, tN0, tN1, tB0, tT0, tTc );
                    }
                    else
                    {
                        BELFEM_ERROR( false, "unknown resistance law: %s", tType.c_str());
                    }
                }

                if ( aInput->key_exists( "permeability" ))
                {
                    // check type of permeability law
                    const string & tType = aInput->get_string( "permeability" );

                    if ( tType == "constant" )
                    {
                        if ( aInput->key_exists( "mu_r" ))
                        {
                            aMaterial->set_mu_r( aInput->get_value( "mu_r", "-" ).first );
                        }
                    }
                    else if ( tType == "textfile" )
                    {
                        const string & tPath = aInput->get_string( "path" );
                        const string tUnitB = aInput->key_exists( "unitb" ) ? aInput->get_string( "unitb" )
                                                                            : "T";
                        const string tUnitH = aInput->key_exists( "unith" ) ? aInput->get_string( "unith" )
                                                                            : "A/m";

                        const value tMaxB = aInput->get_value( "maxb", "T" );

                        real tM = BELFEM_QUIET_NAN;
                        aMaterial->set_nu_s( this->read_bhfile( tPath, tUnitB, tUnitH, tMaxB, tM ) );
                        aMaterial->set_m0( tM );
                    }
                    else
                    {
                        BELFEM_ERROR( false, "unknown permeability law: %s", tType.c_str());
                    }
                }
                if ( aInput->key_exists( "thermal" ) )
                {
                    aMaterial->set_thermal_material( aInput->get_string( "thermal" ) );
                }
                return aMaterial;
            }
        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_layers()
        {
            if( mInputFile.section_exists( "maxwell") )
            {
                if ( mInputFile.section( "maxwell")->section_exists( "tape" ))
                {
                    // get section
                    const input::Section * tSection = mInputFile.section( "maxwell")->section( "tape" );

                    uint tNumLayers = tSection->num_keys();
                    mTapeThicknesses.set_size( tNumLayers );

                    for ( uint l = 0; l < tNumLayers; ++l )
                    {
                        // get label
                        const string & tKey = tSection->key( l );
                        mTapeThicknesses( l ) = tSection->get_value( tKey, "m" ).first;

                        mTapeMaterialLabels.push( tKey  );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::link_materials( DofManager * aDofManager )
        {
            // loop over all blocks
            for( DomainGroup * tGroup : mBlocks )
            {
                if(    tGroup->type() == DomainType::Conductor
                    || tGroup->type() == DomainType::Ferro )
                {
                    // grab material
                    Material * tMaterial = mMaterialMap( tGroup->label() );

                    // check if material has been used already
                    // (each material must exist only once in the kernel)
                    if( ! tMaterial->is_flagged() )
                    {
                        // mark material as used
                        tMaterial->flag();

                        // add material to kernel
                        aDofManager->parent()->add_material( tMaterial );
                    }
                    // link material with block
                    for ( id_t tID: tGroup->groups() )
                    {
                        aDofManager->block( tID )->set_material( tMaterial );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        DomainGroup *
        MaxwellFactory::create_domain_group( const input::Section * aSection )
        {
            // get the label of the block
            const string & tString = aSection->label();

            // check if there is a space
            size_t tSplit = tString.find_first_of( " ", 0 );

            // get type of block
            string tTypeString = tString.substr( 0, tSplit );


            // get name of block
            const string tLabel = tSplit < tString.length() ?
                                  tString.substr( tSplit + 1, tString.length() - tSplit - 1 ) : tTypeString;


            // -  - - - - - - - - - - - - - - -

            // get raw data
            string tLine;

            for ( index_t k = aSection->start(); k < aSection->end(); ++k )
            {
                tLine += " " + mInputFile.line( k );
            }

            // convert to vector
            Vector< id_t > tData ;
            this->read_id_line( tLine, tData );

            // end get raw data
            // -  - - - - - - - - - - - - - - -

            if( string_to_lower( tTypeString ) == "farfield" )
            {
                tTypeString = "symmetry" ;
            }

            // create domain type
            DomainType tType = domain_type( tTypeString );

            switch ( tType )
            {
                case ( DomainType::Conductor ) :
                case ( DomainType::Ferro )  :
                case ( DomainType::Coil ) :
                {
                    return new DomainSolid( tType, tLabel, tData );
                }
                case ( DomainType::Cut ) :
                {
                    // make sure that data can be divided by three
                    BELFEM_ERROR( tData.length() % 3 == 0, "Number of cut data must be divisible by three" );

                    index_t tN = tData.length() / 3;

                    // populate data
                    Vector< id_t > tGroups( tN );
                    Vector< id_t > tPlus( tN );
                    Vector< id_t > tMinus( tN );

                    index_t tCount = 0;

                    for ( uint i = 0; i < tN; ++i )
                    {
                        tGroups( i ) = tData( tCount++ );
                        tPlus( i ) = tData( tCount++ );
                        tMinus( i ) = tData( tCount++ );
                    }

                    // create the group object
                    return new DomainCut( tLabel, tGroups,
                                          tPlus, tMinus, tType  );
                }
                case( DomainType::ThinShell ) :
                {

                    // make sure that data can be divided by three
                    BELFEM_ERROR( tData.length() % 3 == 0, "Number of tape data must be divisible by three" );

                    index_t tN = tData.length() / 3 ;

                    // populate data
                    Vector< id_t > tGroups( tN );
                    Vector< id_t > tPlus( tN );
                    Vector< id_t > tMinus( tN );

                    index_t tCount = 0;

                    for ( uint i = 0; i < tN; ++i )
                    {
                        tGroups( i ) = tData( tCount++ );
                        tPlus( i )   = tData( tCount++ );
                        tMinus( i )  = tData( tCount++ );
                    }

                    // create the group object
                    return new DomainThinShell( tLabel, tGroups,
                                                tPlus, tMinus, tPlus, tType  );
                }
                case ( DomainType::Air ) :
                case ( DomainType::Boundary ) :
                case( DomainType::Symmetry ) :
                case( DomainType::AntiSymmetry ) :
                {
                    return new DomainGroup( tType, tLabel, tData );
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Domain Type" );
                    return nullptr;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_bcs()
        {
            // grab section in input file
            const input::Section * tBCs
                = mInputFile.section("maxwell")->section("boundaryConditions") ;

            // get number of bcs
            index_t tNumBCs = tBCs->num_sections() ;

            Cell< BoundaryCondition * > tCurrentBCs ;
            Cell< BoundaryCondition * > tMagfieldBCs ;

            // check if formulation is h-phi
            bool tIsHA = mInputFile.section("maxwell")->get_string( "formulation") == "h-a" ;

            // loop over all bcs
            for( index_t b = 0 ; b<tNumBCs; ++b )
            {
                const input::Section * tBC = tBCs->section( b );

                // split words of secition
                Cell< string > tWords = string_to_words( tBC->label() );


                // get number of blocks
                index_t tNumWords = tWords.size() ;

                BELFEM_ERROR( tNumWords > 1, "current bc not associated with a block name" );

                // get type of bc
                string tBcType = string_to_lower( tWords( 0 ) );

                // get parameters
                real tAmplitude = BELFEM_QUIET_NAN ;
                real tPeriod    = BELFEM_QUIET_NAN ;
                real tPhase     = 0.0 ;
                real tPenalty   = tBC->key_exists("penalty") ? tBC->get_real( "penalty" ) : 1.0 ;

                const string tScalingType = tBC->key_exists("type") ?  string_to_lower( tBC->get_string("type") ) :
                        "none" ;

                BoundaryConditionScaling tType = BoundaryConditionScaling::UNDEFINED ;
                MagfieldBcType tMagfieldBcType = MagfieldBcType::UNDEFINED ;

                id_t tBearingID = gNoID ;

                // catch special case
                if( tScalingType == "none"  )
                {
                    continue;
                }
                else if( tScalingType == "farfield" ||  tScalingType == "symmetry"  )
                {
                    tMagfieldBcType = MagfieldBcType::Symmetry ;
                }
                else if( tScalingType == "antisymmetry"  )
                {
                    tMagfieldBcType = MagfieldBcType::AntiSymmetry ;
                }
                else if ( tScalingType == "bearing" || tScalingType == "anchor" || tScalingType == "fixedpoint" )
                {
                    tMagfieldBcType = MagfieldBcType::Bearing ;
                    tBearingID = tBC->get_int( "node" ) ;
                }
                else
                {
                    tMagfieldBcType = MagfieldBcType::Wave ;

                    // get current type
                    tType = bc_scaling_type( tScalingType );

                    BELFEM_ERROR( tType == BoundaryConditionScaling::Constant ||
                                 ( tBC->key_exists("period") && ! tBC->key_exists("frequency") ) ||
                                 ( ! tBC->key_exists("period") && tBC->key_exists("frequency") ),
                                 "Must impose either period or frequency, not both!" );


                    // read parameters from input file
                    switch ( tType )
                    {
                        case ( BoundaryConditionScaling::Constant ) :
                        {
                            break;
                        }
                        case ( BoundaryConditionScaling::Ramp ) :
                        {
                            tPeriod = tBC->key_exists( "period" ) ? tBC->get_value( "period", "s" ).first :
                                      1.0 / tBC->get_value( "frequency", "Hz" ).first;
                            break;
                        }
                        case ( BoundaryConditionScaling::Sine ) :
                        case ( BoundaryConditionScaling::Sawtooth ) :
                        case ( BoundaryConditionScaling::Square ) :
                        case ( BoundaryConditionScaling::Triangle ) :
                        case ( BoundaryConditionScaling::Sigmoid ) :
                        {
                            tPeriod = tBC->key_exists( "period" ) ? tBC->get_value( "period", "s" ).first :
                                      1.0 / tBC->get_value( "frequency", "Hz" ).first;

                            if ( tBC->key_exists( "phase" ) )
                            {
                                tPhase = tBC->get_value( "phase", "deg" ).first;
                            }
                            break;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "unsupported bc type" );
                        }
                    }

                }

                // check if this is a current bc
                if( tBcType == "current" )
                {
                    BELFEM_ERROR( tMagfieldBcType == MagfieldBcType::Wave , "illegal bctype set for current");

                    tAmplitude = tBC->get_value("amplitude", "A" ).first ;

                    BoundaryConditionImposing tImposing = tIsHA ? BoundaryConditionImposing::Neumann : BoundaryConditionImposing::Lambda ;

                    for( index_t w=1; w<tNumWords; ++w )
                    {
                        // get label
                        const string tLabel = string_to_lower( tWords( w ) );

                        // create new current bc
                        MaxwellBoundaryConditionCurrent * tCurrent = new MaxwellBoundaryConditionCurrent(
                                tImposing,
                                tWords( w ) ) ;

                        // check if this is a negative bc
                        real tSign = ( tLabel.find( "minus") < tLabel.length() || tLabel.find("-") < tLabel.length() ) ?
                                -1.0 : 1.0 ;

                        // set parameters for bc
                        tCurrent->set( tType, tSign * tAmplitude, tPeriod, tPhase );

                        // add bc to container
                        tCurrentBCs.push( tCurrent );

                        // add bc to map
                        mCurrentBcMap[ tWords( w ) ] = tCurrent ;
                    }

                }
                else if ( tBcType == "magfield" )
                {
                    if( tMagfieldBcType == MagfieldBcType::Bearing )
                    {
                        MaxwellBoundaryConditionMagfield * tMagfield = new MaxwellBoundaryConditionMagfield(
                                BoundaryConditionImposing::Dirichlet,
                                tWords( 1 ) );
                        tMagfield->set_bearing( tBearingID );
                        tMagfieldBCs.push( tMagfield );
                        mMagfieldBcMap[ tWords( 1 ) ] = tMagfield;
                    }
                    else if ( tMagfieldBcType == MagfieldBcType::Wave )
                    {
                        // get string
                        string tString = tBC->get_string( "amplitude" );

                        size_t i = tString.find( "(" ) + 1;
                        size_t j = tString.find( ")", i ) - 1;

                        // get value string
                        string tValStr = clean_string(
                                search_and_replace( tString.substr( i, j - i ), ",", " " ));

                        // convert string to numbers
                        Cell< string > tVals = string_to_words( tValStr );


                        // convert string to vector
                        Vector< real > tValues( tVals.size() );
                        for ( uint k = 0; k < tVals.size(); ++k )
                        {
                            tValues( k ) = std::stod( tVals( k ) );
                        }

                        // get imposing type
                        BoundaryConditionImposing tImposing =  ( norm( tValues ) < BELFEM_EPSILON ) || tIsHA ?
                                BoundaryConditionImposing::Dirichlet :
                                BoundaryConditionImposing::Weak ;

                        // get unit
                        string tUnit = clean_string( tString.substr( j + 2 ));

                        // make sure that unit is correct
                        value tScale = unit_to_si( tUnit );
                        BELFEM_ERROR( check_unit( tScale, "T" ), "invalid unit, expect T or equivalent" );

                        // scale the values
                        tValues *= tScale.first;

                        // create new bcs
                        for( index_t w=1; w<tNumWords; ++w )
                        {
                            // get label
                            const string tLabel = string_to_lower( tWords( w ) );

                            MaxwellBoundaryConditionMagfield * tMagfield = new MaxwellBoundaryConditionMagfield(
                                    tImposing,
                                    tWords( w ) );

                            // check sign
                            if ( tLabel.find( "minus" ) < tLabel.length() || tLabel.find( "-" ) < tLabel.length())
                            {
                                Vector< real > tMinus( tValues );
                                tMinus *= -1.0;
                                tMagfield->set( tType, tMinus, tPeriod, tPhase );
                            }
                            else
                            {
                                tMagfield->set( tType, tValues, tPeriod, tPhase );
                            }

                            // set penalty of magfield bc
                            tMagfield->set_penalty( tPenalty );

                            tMagfieldBCs.push( tMagfield );
                            mMagfieldBcMap[ tWords( w ) ] = tMagfield;
                        }
                    }
                    else
                    {
                        // get imposing type
                        BoundaryConditionImposing tImposing =  BoundaryConditionImposing::Weak ;

                        for( index_t w=1; w<tNumWords; ++w )
                        {
                            // get label
                            const string tLabel = string_to_lower( tWords( w ));

                            MaxwellBoundaryConditionMagfield * tMagfield
                                    = new MaxwellBoundaryConditionMagfield( tImposing, tWords( w ));

                            tMagfield->set( tMagfieldBcType );

                            // set penalty of magfield bc
                            tMagfield->set_penalty( tPenalty );

                            tMagfieldBCs.push( tMagfield );
                            mMagfieldBcMap[ tWords( w ) ] = tMagfield;
                        }
                    }


                } // end if magfield
                else if ( comm_rank() == 0 )
                {
                    std::cout << "Warning: bc type " << tWords( 0 )  << " not implemented"  << std::endl ;
                }


            } // end bc loop


            // BCs must be in this order
            for( BoundaryCondition * tBc : tCurrentBCs )
            {
                mBoundaryConditions.push( tBc );
            }
            for( BoundaryCondition * tBc : tMagfieldBCs )
            {
                mBoundaryConditions.push( tBc );
            }
        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_kernel()
        {
            // check how many fields we have
            this->read_field_switches() ;

            // count linear fields : field one for formulation
            uint tNumLinearFields = 1 ;

            // check if we compute J
            if( mHaveJ ) ++tNumLinearFields ;

            // check if we compute B

            if( mHaveB )
            {
                ++tNumLinearFields ;
                if( is_maxwell_hphi( mFormulation ) )
                {
                    ++tNumLinearFields ;
                }
            }

            uint tNumFields = tNumLinearFields ;

            // check if we have thermal
            if( mHaveThermal ) ++tNumFields ;

            // mech requires mech + sigma postprocess
            if( mHaveMech ) tNumFields += 2 ;

            // set flags for linear enforcement
            Vector< uint > tFlags( tNumFields, 0 );

            for( uint f=0; f<tNumLinearFields; ++f )
            {
                tFlags( f ) = mEnforceLinear ? 1 : 0 ;
            }

            // write flags to parameters
            mParameters->enforce_linear( tFlags );

            // mParameters->set_sideset_integration_orders( 7 );
            // mParameters->set_block_integration_orders( 7 );

            if( comm_rank() == 0 )
            {
                // count selected blocks ( blocks for visualization are not included )
                uint tCount = 0;
                for ( uint b = 0; b < mMesh->number_of_blocks(); ++b )
                {
                    if ( mMesh->blocks()( b )->id() <= mMaxBlockID )
                    {
                        ++tCount;
                    }
                }
                Vector< id_t > tSelectedBlocks( tCount );
                tCount = 0;
                for ( uint b = 0; b < mMesh->number_of_blocks(); ++b )
                {
                    id_t tID = mMesh->blocks()( b )->id();
                    if ( tID <= mMaxBlockID )
                    {
                        tSelectedBlocks( tCount++ ) = tID;
                    }
                }
                mParameters->select_blocks( tSelectedBlocks );
            }

            // create the kernel
            mKernel = new Kernel( mParameters );

            // claim parameter ownership
            mKernel->claim_parameter_ownership( true );
        }

//------------------------------------------------------------------------------

        IWG_Maxwell *
        MaxwellFactory::create_iwg( const IwgType aType )
        {
            switch( aType )
            {
                case( IwgType::MAXWELL_HPHI_TRI3 ) :
                {
                    return new IWG_Maxwell_HPhi_Tri3() ;
                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TRI6 ) :
                {
                    return new IWG_Maxwell_HPhi_Tri6() ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "unsoppurted type");
                    return nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::set_tape_materials( DofManager * aDofManager )
        {
            if( mTapeMaterialLabels.size() > 0 )
            {
                // loop over all sidesets on the mesh
                for( mesh::SideSet * tMeshSideSet : mKernel->mesh()->sidesets() )
                {
                    // check if sideset exists on dof manager
                    if( aDofManager->sideset_exists( tMeshSideSet->id() ) )
                    {
                        // grab sideset
                        SideSet * tSideSet = aDofManager->sideset( tMeshSideSet->id() );

                        // check if sideset is a thin shell
                        if( tSideSet->domain_type() == DomainType::ThinShell )
                        {
                            tSideSet->set_thin_shells( mTapeMaterials, mTapeThicknesses );
                        }
                    }
                }


                /*
                 * the lines below set the physical tag for the elements so that the
                 * material can be identified in ParaView. We only do this for the
                 * main proc
                 */

                if( comm_rank() > 0 )
                {
                    return;
                }
                // get the block IDs
                const Vector< id_t > & tBlockIDs = aDofManager->iwg()->ghost_blocks() ;

                uint tNumBlocks = tBlockIDs.length() ;

                // create a material map
                Map< string, uint > tMap ;

                uint tCount = 0 ;
                for( const string & tLabel : mTapeMaterialLabels )
                {
                    if( ! tMap.key_exists( tLabel ) )
                    {
                        tMap[ tLabel ] = ++tCount ;
                    }
                }

                for( uint b=0; b<tNumBlocks; ++b )
                {
                    // grab block on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockIDs( b ));

                    // get the material number
                    uint tNumber = tMap[ mTapeMaterialLabels( b ) ];

                    // set the physical tag of the elements
                    for ( mesh::Element * tElement: tBlock->elements())
                    {
                        tElement->set_physical_tag( tNumber );
                    }
                }

            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_magnetic_field()
        {
            // make a new maxwell equation
            IWG_Maxwell * tMaxwell = this->create_iwg( mFormulation );



            // the IDs for the cuts differ from the sidesets on the mesh,
            // we need to change this
            for( uint s=0; s<mSideSetIDs.length(); ++s )
            {
                if( mSideSetTypes( s ) == DomainType::Cut )
                {
                    mSideSetIDs( s ) = mSideSetToCutMap( mSideSetIDs( s ) );
                }
            }

            // wait for other procs
            comm_barrier();

            // flag curved elements on kernel mesh
            mKernel->mesh()->flag_curved_elements() ;

            // wait for other procs
            comm_barrier();

            tMaxwell->set_blocks( mBlockIDs, mBlockTypes );

            tMaxwell->set_sidesets( mSideSetIDs, mSideSetTypes );

            // also set the BC types
            for( BoundaryCondition * tBC : mBoundaryConditions )
            {
                // check if this is a magfield bc
                if( tBC->physics() == BoundaryConditionPhysics::Magfield )
                {
                    // loop over all sideset ids
                    for( id_t tID : tBC->sidesets() )
                    {
                        tMaxwell->set_magfield_bc_type( tID, reinterpret_cast< MaxwellBoundaryConditionMagfield * >( tBC )->subtype() );
                    }
                }
            }

            // the following lines are needed to tell the IWG from which sideset the dofs for the ghosts
            // are copied
            mGhostMaster = 0 ;

            for( uint k=0; k<mSideSetIDs.length(); ++k )
            {
                if( mSideSetTypes( k ) == DomainType::ThinShell )
                {
                    mGhostMaster = mSideSetIDs( k );
                    break ;
                }
            }

            if( mGhostMaster > 0 )
            {
                tMaxwell->set_ghost_sidesets( mGhostMaster, mGhostSideSets );
                tMaxwell->set_ghost_blocks( mGhostBlocks );

                // the maps associates the sideset facets with the index in the ghost block
                tMaxwell->create_ghost_map( mKernel->mesh(), mThinShellSideSets, mKernel->comm_table() );
            }

            for ( BoundaryCondition * tBC: mBoundaryConditions )
            {
                tMaxwell->add_boundary_condition( tBC );
            }

            // for thin shells
            tMaxwell->set_thin_shell_link_mode( SideSetDofLinkMode::MasterAndSlave );

            mKernel->add_equation( tMaxwell );

            // create dof manager for magnetic field
            mMagneticField = mKernel->create_field( tMaxwell );

            // create data arrays on mesh
           // mMagneticField->create_fields( tMaxwell );

            // link field with materials
            this->link_materials( mMagneticField );

            // set the tape materials, if tapes exist
            this->set_tape_materials( mMagneticField );

            // set timestepping info
            tMaxwell->set_timestep( mTimeStep );


            // set the tau parameter for stabilization
            tMaxwell->set_tau( mTau );

            // todo, must be always theta=1
            BELFEM_ERROR( mTheta = 1.0 , "must use Euler Implicit here!");
            tMaxwell->set_euler_method( 1.0 );


            // set the solver of this field based on the input
            this->set_solver( mMagneticField,
                              mInputFile.section("maxwell")
                              ->section("solver")->section("linear") );

            // add boundary conditions to field
            for( BoundaryCondition * tBC : mBoundaryConditions )
            {
                tBC->set_dof_manager( mMagneticField );
            }
            mOwnBoundaryConditions = false ;

            // set the current blocks, needed for elementwise current postprocess
            /*tMaxwellOld->set_current_blocks( mCurrentBlocks );

            // check the biot-savart flag
            if( mInputFile.section("maxwell")->key_exists("biotsavart") )
            {
                mComputeBiotSavart =  mInputFile.section("maxwell")->get_bool("biotsavart");

                if( mComputeBiotSavart )
                {
                    BELFEM_ERROR( tMaxwellOld->type() == IwgType::MAXWELL_HA2D ||
                                         tMaxwellOld->type() == IwgType::MAXWELL_HAAX ||
                                         tMaxwellOld->type() == IwgType::MAXWELL_HA3D,
                                 "Error in input file: Biot-Savart can only be computed if we use a h-a formulation" );

                }
            } */

            // write the domain types into mMagneticField
            this->set_domain_types();

            mMagneticField->fix_node_dofs_on_ghost_sidesets() ;

            comm_barrier() ;
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::set_domain_types()
        {
            for( uint b=0; b<mBlockIDs.length(); ++b )
            {
                if ( mMagneticField->block_exists( mBlockIDs( b ) ) )
                {
                    mMagneticField->block(
                            mBlockIDs( b ) )->set_domain_type( mBlockTypes( b ) );

                }
            }

            for( uint s=0; s<mSideSetIDs.length(); ++s )
            {
                if ( mMagneticField->sideset_exists( mSideSetIDs( s ) ) )
                {
                    mMagneticField->sideset( mSideSetIDs( s ) )
                            ->set_domain_type( mSideSetTypes( s ) );

                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::set_solver(  DofManager * aDofMaganer,
                                     const input::Section * aSection )
        {
            // get solver type
            SolverType tSolverType = solver_type( aSection->get_string("library") );

            aDofMaganer->set_solver( tSolverType );

            if( tSolverType == SolverType::PETSC )
            {
                KrylovMethod tMethod = aSection->key_exists("krylovmethod") ?
                                       krylov_method( aSection->get_string("krylovmethod" )) :
                                       KrylovMethod::GMRES ;

                Preconditioner tPreconditioner = aSection->key_exists("preconditioner") ?
                                                 preconditioner( aSection->get_string("preconditioner") ) :
                                                 Preconditioner::ASM ;

                real tEpsilon = aSection->key_exists("epsilon") ? aSection->get_real( "epsilon") : 1e-8 ;

                aDofMaganer->solver()->set_petsc( tPreconditioner, tMethod, tEpsilon );

            }
            if( tSolverType == SolverType::MUMPS )
            {
                SerialReodrdering   tSerial   = aSection->key_exists( "serialreordering") ?
                                                serial_reordering( aSection->get_string( "serialreordering") ) :
                                                SerialReodrdering::AUTOMATIC ;

                ParallelReodrdering tParallel = aSection->key_exists("parallelreordering") ?
                                                parallel_reordering( aSection->get_string( "parallelreordering") ) :
                                                ParallelReodrdering::AUTOMATIC ;

                aDofMaganer->solver()->set_mumps_reordering( tSerial, tParallel );

                BlockLowRanking    tBLR       = aSection->key_exists("blocklowranking") ?
                                                block_low_ranking( aSection->get_string("blocklowranking") ) :
                                                BlockLowRanking::Automatic ;

                real tEpsilon = aSection->key_exists("blrepsilon") ? aSection->get_real( "blrepsilon") : 0.0 ;

                aDofMaganer->solver()->set_mumps_blr( tBLR, tEpsilon );

            }
        }

//------------------------------------------------------------------------------

        NonlinearSettings::NonlinearSettings( const real aMinIter,
                                              const real aMaxIter,
                                              const real aPicardOmega,
                                              const real aPicardEpsilon,
                                              const real aNewtonOmega,
                                              const real aNewtonEpsilon ) :
                           minIter( aMinIter ),
                           maxIter( aMaxIter ),
                           picardOmega( aPicardOmega ),
                           picardEpsilon( aPicardEpsilon ),
                           newtonOmega( aNewtonOmega ),
                           newtonEpsilon( aNewtonEpsilon )
        {

        }

//------------------------------------------------------------------------------

        NonlinearSettings
        MaxwellFactory::nonlinear_settings()
        {
            // get key
            const input::Section * tSection = mInputFile.section("maxwell")->section("solver")->section("nonlinear");

            uint tMinIter = tSection->key_exists( "miniter") ? tSection->get_int("miniter") : 3 ;
            uint tMaxIter = tSection->key_exists( "maxiter") ? tSection->get_int("maxiter") : 100 ;

            real tNewtonEpsilon = tSection->section_exists("newton") ?
                                  tSection->section("newton")->key_exists("epsilon") ?
                                  tSection->section("newton")->get_real("epsilon") : 1e-7: 1e-7 ;

            real tNewtonOmega = tSection->section_exists("newton") ?
                                tSection->section("newton")->key_exists("omega") ?
                                tSection->section("newton")->get_real("omega") : 0.1 : 0.1 ;

            real tPicardOmega = tSection->section_exists("picard") ?
                                tSection->section("picard")->key_exists("omega") ?
                                tSection->section("picard")->get_real("omega") : 0.1 : 0.1 ;



            real tPicardEpsilon = tSection->section_exists("picard") ?
                                tSection->section("picard")->key_exists("epsilon") ?
                                tSection->section("picard")->get_real("epsilon") : tNewtonEpsilon : tNewtonEpsilon ;





            NonlinearSettings aData(
                    tMinIter,
                    tMaxIter,
                    tPicardOmega,
                    tPicardEpsilon,
                    tNewtonOmega,
                    tNewtonEpsilon );

            return aData ;
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_current_projector()
        {

            // collect all superconducting blocks
            uint tCount = 0.0 ;
            for(  uint b=0; b<mBlockTypes.size(); ++b )
            {
                if( mBlockTypes( b ) == DomainType::Conductor )
                {
                    ++tCount ;
                }
            }

            if( tCount > 0 )
            {
                Vector< id_t > tBlocks( tCount );
                tCount = 0;
                for ( uint b = 0; b < mBlockTypes.size(); ++b )
                {
                    if ( mBlockTypes( b ) == DomainType::Conductor )
                    {
                        tBlocks( tCount++ ) = mBlockIDs( b );
                    }
                }

                Cell< DomainType > tTypes( tCount, DomainType::Conductor );

                // get the element type
                ElementType tType = this->element_type( tBlocks( 0 ) );

                IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Current( tType, mAlphaJ );

                // set the blocks of the IWG
                tIWG->set_blocks( tBlocks, tTypes );

                // add iwg to kernel
                this->add_postprocessor_to_kernel( tIWG );

                // fix group types of IWG
                tIWG->update_group_types();
            }
        }

//------------------------------------------------------------------------------

        /* alternate version with one single projector
               void
               MaxwellFactory::create_magfield_projector()
               {
                   // count fields
                   uint tCount = 0 ;

                   for( uint s=0; s<mSideSetIDs.length(); ++s )
                   {
                       ++tCount ;
                   }

                   Vector< id_t > tCuts( tCount );
                   tCount = 0 ;
                   for( uint s=0; s<mSideSetIDs.length(); ++s )
                   {
                       if( mSideSetTypes( s ) == DomainType::Cut )
                       {
                           tCuts( tCount++ )  = mSideSetIDs( s );
                       }
                   }

                   // determine mode for projectors
                   Magfield_L2_Mode tMode = Magfield_L2_Mode::UNDEFINED ;

                   if( is_maxwell_hphi( mFormulation ) )
                   {
                       tMode = Magfield_L2_Mode::HPHI ;
                   }
                   else if ( is_maxwell_ha( mFormulation ) )
                   {
                       tMode = Magfield_L2_Mode::HA ;
                   }


                   // create the iwg
                   IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Magfield(
                           this->element_type( mBlockIDs( 0 ) ),
                           tMode,
                           mAlphaB );


                   // set the blocks of the IWG
                   tIWG->set_blocks( mBlockIDs, mBlockTypes );

                   if( is_maxwell_hphi( mFormulation ) && tCount > 0 )
                   {
                       Cell< DomainType > tCutTypes( tCuts.length(), DomainType::Cut );
                       tIWG->set_sidesets( tCuts, tCutTypes );
                   }

                   // add iwg to kernel
                   this->add_postprocessor_to_kernel( tIWG );

                   // fix group types of IWG
                   tIWG->update_group_types();

               } */

        void
        MaxwellFactory::create_magfield_projector()
        {
            // count fields
            uint tAirCount = 0 ;
            uint tFerroCount = 0 ;
            uint tConductorCount = 0 ;
            uint tCutCount = 0 ;

            for(  DomainType tType : mBlockTypes )
            {
                switch( tType )
                {
                    case( DomainType::Air ) :
                    case( DomainType::Coil ) :
                    {
                        ++tAirCount ;
                        break ;
                    }
                    case( DomainType::Ferro ) :
                    {
                        ++tFerroCount ;
                        break ;
                    }
                    case( DomainType::Conductor ) :
                    {
                        ++tConductorCount ;
                        break ;
                    }
                    default :
                    {
                        break ;
                    }
                }
            }
            for(  DomainType tType : mSideSetTypes )
            {
                if( tType == DomainType::Cut )
                {
                    ++tCutCount ;
                }
            }

            Vector< id_t > tAirBlocks( tAirCount );
            Vector< id_t > tFerroBlocks( tFerroCount );
            Vector< id_t > tConductorBlocks( tConductorCount );
            Vector< id_t > tCuts( tCutCount );

            tAirCount       = 0 ;
            tFerroCount     = 0 ;
            tConductorCount = 0 ;
            tCutCount       = 0 ;

            for(  uint b=0; b<mBlockTypes.size(); ++b )
            {
                switch( mBlockTypes( b ) )
                {
                    case( DomainType::Air ) :
                    case( DomainType::Coil ) :
                    {
                        tAirBlocks( tAirCount++ ) =  mBlockIDs( b );
                        break ;
                    }
                    case( DomainType::Ferro ) :
                    {
                        tFerroBlocks( tFerroCount++ ) =   mBlockIDs( b );
                        break ;
                    }
                    case( DomainType::Conductor ) :
                    {
                        tConductorBlocks( tConductorCount++ ) = mBlockIDs( b );
                        break ;
                    }
                    default :
                    {
                        break ;
                    }
                }
            }
            for( uint s=0; s<mSideSetIDs.length(); ++s )
            {
                if( mSideSetTypes( s ) == DomainType::Cut )
                {
                    tCuts( tCutCount++ )  = mSideSetIDs( s );
                }
            }

            // determine mode for projectors
            Magfield_L2_Mode tMode = Magfield_L2_Mode::UNDEFINED ;

            if( is_maxwell_hphi( mFormulation ) )
            {
                tMode = Magfield_L2_Mode::HPHI ;
            }
            else if ( is_maxwell_ha( mFormulation ) )
            {
                tMode = Magfield_L2_Mode::HA ;
            }

            // crate ferro projector
            if( tFerroCount > 0 )
            {
                // create the iwg
                IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Magfield(
                        this->element_type( tFerroBlocks( 0 )),
                        tMode,
                        mAlphaB );

                Cell< DomainType > tTypes( tFerroBlocks.length(), DomainType::Ferro );

                // set the blocks of the IWG
                tIWG->set_blocks( tFerroBlocks, tTypes );

                // add iwg to kernel
                this->add_postprocessor_to_kernel( tIWG );

                // fix group types of IWG
                tIWG->update_group_types();
            }

            // create conducting projector
            if( tConductorCount > 0 )
            {
                // create the iwg
                IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Magfield(
                        this->element_type( tConductorBlocks( 0 ) ),
                        tMode,
                        mAlphaB );

                Cell< DomainType > tTypes( tConductorBlocks.length(), DomainType::Conductor );

                // set the blocks of the IWG
                tIWG->set_blocks( tConductorBlocks, tTypes );

                // add iwg to kernel
                this->add_postprocessor_to_kernel( tIWG );

                // fix group types of IWG
                tIWG->update_group_types();
            }

            // create air projector
            if( tAirCount > 0 )
            {

                // create the iwg
                IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Magfield(
                        this->element_type( tAirBlocks( 0 )),
                        tMode,
                        mAlphaB );

                Cell< DomainType > tTypes( tAirBlocks.length(), DomainType::Air );

                // set the blocks of the IWG
                tIWG->set_blocks( tAirBlocks, tTypes );

                if( is_maxwell_hphi( mFormulation ) && tCutCount > 0 )
                {
                    Cell< DomainType > tCutTypes( tCuts.length(), DomainType::Cut );
                    tIWG->set_sidesets( tCuts, tCutTypes );
                }

                // add iwg to kernel
                DofManager * tProjector = this->add_postprocessor_to_kernel( tIWG );

                // fix group types of IWG
                tIWG->update_group_types();
            }
        }

// -----------------------------------------------------------------------------

        ElementType
        MaxwellFactory::element_type( const id_t aBlockID )
        {
            uint tVal = 0;
            // check element types
            if ( mKernel->is_master())
            {
                tVal = ( uint ) mMesh->block( aBlockID )->element_type();
                send( mKernel->comm_table(), tVal );
            }
            else
            {
                receive( mKernel->master(), tVal );
            }
            return ( ElementType ) tVal;
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::blacklist_dirichlet_sidesets()
        {
            // backup copy of sets
            Vector< id_t > tSideSetIDs( mSideSetIDs );

            // backup sideset types
            Cell< DomainType > tSideSetTypes( mSideSetTypes.size(), DomainType::UNDEFINED );
            uint tCount = 0 ;
            for( DomainType tType : mSideSetTypes )
            {
                tSideSetTypes( tCount++ ) = tType ;
            }

            // now we loop over all BCs and create a blacklist
            Map< id_t, uint > tBlacklist ;

            tCount = 0 ;

            for( BoundaryCondition * tBC : mBoundaryConditions )
            {
                if( tBC->physics() == BoundaryConditionPhysics::Magfield && tBC->imposing() == BoundaryConditionImposing::Dirichlet )
                {
                    // add sidesets to blacklist
                    for( id_t tID : tBC->sidesets() )
                    {
                        tBlacklist[ tID ] = tCount ++ ;
                    }
                }
            }

            // check if sidesets have been blacklisted
            if( tCount > 0 )
            {
                uint tNumAllSidesets = tSideSetIDs.length() ;

                uint tNumSideSets = tNumAllSidesets - tCount ;
                mSideSetIDs.set_size( tNumSideSets, 0 );
                mSideSetTypes.set_size( tNumSideSets, DomainType::UNDEFINED );

                tCount = 0 ;

                // only add sets that are not blacklisted
                for( uint s=0; s<tNumAllSidesets; ++s )
                {
                    if( ! tBlacklist.key_exists( tSideSetIDs( s ) ) )
                    {
                        mSideSetIDs( tCount ) = tSideSetIDs( s );
                        mSideSetTypes( tCount++ ) = tSideSetTypes( s );
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        DofManager *
        MaxwellFactory::add_postprocessor_to_kernel( IWG_Maxwell * aEquation )
        {
            // add equation to kernel
            mKernel->add_equation( aEquation );

            // create dof projector
            DofManager * aProjector = mKernel->create_field( aEquation );

            // link field with materials
            this->link_materials( aProjector );

            // make sure that data arrays on mesh exist
            aProjector->create_fields( aEquation );

            // set the solver of this field based on the input
            this->set_solver( aProjector,
                              mInputFile.section("maxwell")
                                      ->section("solver")->section("linear") );

            // catch problem with PARDISO
            if( aProjector->solver()->type() ==  SolverType::PARDISO )
            {
#ifdef BELFEM_STRUMPACK
                aProjector->set_solver( SolverType::STRUMPACK );
#elif BELFEM_MUMPS
                aProjector->set_solver( SolverType::MUMPS );
#else
                aProjector->set_solver( SolverType::UMFPACK );
#endif
            }

            // no BLR for L2, otherwise crash!
            aProjector->solver()->set_mumps_blr( BlockLowRanking::Off, 0.0 );

            aProjector->fix_node_dofs_on_ghost_sidesets() ;

            // add projector to list
            mMagneticField->postprocessors().push( aProjector );

            return aProjector ;
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::create_topology()
        {
            // moved to load mesh
            // create the domain groups based on the input file
            // this->create_domain_groups();
            // this->select_blocks();

            // populate physical tags in mesh
            if( comm_rank() == mMesh->master() )
            {
                this->write_block_phystags_to_mesh();
                this->fix_facet_masters() ;
                this->create_interfaces();
                this->create_symmetries( false ) ;
                this->create_symmetries( true ) ;
                this->write_sideset_phystags_to_mesh();
            }
            else
            {
                this->create_interfaces();
                this->create_symmetries( false ) ;
                this->create_symmetries( true ) ;
            }

            // crate the layers of the tape
            this->create_layers() ;

            // with the blocks in place, we can also create the materials
            this->create_materials() ;
            this->select_bcs_and_cuts() ;
            this->select_sidesets() ;
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_domain_groups()
        {
            // get blocks
            const input::Section * tBlocks = mInputFile.section( "topology" )->section( "blocks" );

            // get the number of blocks
            index_t tNumBlocks = tBlocks->num_sections() ;

            // loop over all blocks
            for( index_t b=0; b<tNumBlocks; ++b )
            {
                // create the new group
                DomainGroup * tGroup
                        = this->create_domain_group( tBlocks->section( b ) );

                // add group to list
                mBlocks.push( tGroup );
                mBlockMap[ tGroup->label() ] = tGroup ;
            }

            Map< id_t, DomainType > tSideSetTypeMap ;

            // get sidesets
            const input::Section * tSideSets = mInputFile.section( "topology" )->section( "sidesets" );

            // get the number of sidesets
            index_t tNumSideSets = tSideSets->num_sections() ;

            Cell< DomainGroup * > tNedelecGroups ;

            // check if this is a h-phi formulation. usually, we would use
            //  is_maxwell_hphi( mFormulation ) here, but mFormulation is not set jet;
            string tKey = string_to_lower(
                    mInputFile.section("maxwell")->get_string("formulation") );
            bool tIsHPhi = ( tKey == "h-phi" || tKey == "h-φ" || tKey == "h-Φ" ) ;

            // loop over all sidesets
            for( index_t s=0; s<tNumSideSets; ++s )
            {
                // create the new group
                DomainGroup * tGroup
                        = this->create_domain_group( tSideSets->section( s ) );

                switch( tGroup->type() )
                {
                    case( DomainType::ThinShell ) :
                    {
                        tNedelecGroups.push( tGroup );
                        mTapes.push( tGroup );
                        break ;
                    }
                    case( DomainType::Cut ) :
                    {
                        if( tIsHPhi )
                        {
                            mCuts.push( tGroup );
                            mCutMap[ tGroup->label() ] = tGroup ;
                        }
                        else
                        {
                            // ignore cuts if this is h-a
                            delete tGroup ;
                        }

                        break ;
                    }
                    case( DomainType::Boundary ) :
                    {
                        // add group to list
                        mBoundaries.push( tGroup );
                        mBoundaryMap[ tGroup->label() ] = tGroup;
                        break ;
                    }
                    case( DomainType::Symmetry ) :
                    {
                        // add group to list
                        mRawSymmetries.push( tGroup );
                        break ;
                    }
                    case( DomainType::AntiSymmetry ) :
                    {
                        // add group to list
                        mRawAntiSymmetries.push( tGroup );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "invalid sideset type" );
                    }
                }
            }

            // identify nedelec sidesets
            index_t tCount = 0 ;

            // count nedelec sidesets
            for(   DomainGroup * tGroup : tNedelecGroups )
            {
                tCount += tGroup->groups().length() ;
            }

            mNedelecSideSets.set_size( tCount );

            // reset counter
            tCount = 0 ;

            // populate container
            for( DomainGroup * tGroup : tNedelecGroups )
            {
                for( id_t g : tGroup->groups() )
                {
                    mNedelecSideSets( tCount++ ) = g ;
                }
            }

            // check input
            unique( mNedelecSideSets ) ;

            BELFEM_ERROR( mNedelecSideSets.length() == tCount, "List of thin shell sidesets is not unique" );

        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::select_blocks()
        {
            // initialize counter
            index_t tCount = 0 ;

            // count used blocks
            for( DomainGroup * tGroup : mBlocks )
            {
                tCount += tGroup->groups().length() ;
            }

            // allocate memory
            mBlockIDs.set_size( tCount );

            // reset counter
            tCount = 0 ;

            // counter for nedelec blocks
            uint tNedelecCount = 0 ;

            // temporary map for types
            Map< id_t, DomainType > tBlockTypeMap ;

            // loop over used blocks
            for( DomainGroup * tGroup : mBlocks )
            {
                const Vector< id_t > & tGroups = tGroup->groups() ;

                for( index_t k = 0 ; k<tGroups.length(); ++k  )
                {
                    mBlockIDs( tCount++ ) = tGroups( k );
                    if ( ! tBlockTypeMap.key_exists( tGroups( k ) ) )
                    {
                        if( tGroup->type() == DomainType::Conductor )
                        {
                            ++tNedelecCount ;
                        }
                        tBlockTypeMap[ tGroups( k ) ] = tGroup->type() ;
                    }
                }
            }

            // make sure that everything is unique
            unique( mBlockIDs );
            BELFEM_ERROR( mBlockIDs.length() == tCount, "List of Blocks is not unique" );

            // populate type container
            mBlockTypes.set_size( tCount, DomainType::UNDEFINED );
            tCount = 0 ;
            for( id_t tID : mBlockIDs )
            {
                mBlockTypes( tCount++ ) = tBlockTypeMap( tID );
            }

            // allocate memory
            mNedelecBlocks.set_size( tNedelecCount );
            tCount = 0 ;
            for( id_t tID : mBlockIDs )
            {
                if( tBlockTypeMap( tID ) == DomainType::Conductor )
                {
                    mNedelecBlocks( tCount++ ) = tID ;
                }
            }

            BELFEM_ASSERT( tCount == tNedelecCount, "Error while counting Nedelec Blocks" );

        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::write_block_phystags_to_mesh()
        {
            // reset the tags of the elements
            Cell< mesh::Element * > & tAllElements = mMesh->elements() ;
            for( mesh::Element * tElement : tAllElements )
            {
                tElement->set_physical_tag( 0 );
            }


            uint tNumZeros = ( uint ) std::log10( mMesh->number_of_blocks() );
            tNumZeros < 1 ? 1 : tNumZeros ;

            string tFormat = "%0" + std::to_string( tNumZeros ) + "u_%s" ;

            for( mesh::Block * tBlock : mMesh->blocks() )
            {
                tBlock->label() = sprint( tFormat.c_str(), ( uint ) tBlock->id(), "Block" );
            }

            for ( DomainGroup * tGroup : mBlocks )
            {
                // convert domain type to tag
                uint tTag = static_cast< uint > ( tGroup->type());

                for ( id_t tBlockID: tGroup->groups() )
                {
                    if ( mMesh->block_exists( tBlockID ))
                    {
                        // get block on mesh
                        mesh::Block * tBlock = mMesh->block( tBlockID );

                        // loop over all elements of this block
                        Cell< mesh::Element * > & tElements = tBlock->elements();

                        for ( mesh::Element * tElement: tElements )
                        {
                            tElement->set_physical_tag( tTag );
                        }

                        tBlock->label() = sprint( tFormat.c_str(),
                                                  ( uint ) tBlock->id(),
                                                  tGroup->label().c_str() );
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::write_sideset_phystags_to_mesh()
        {
            for( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                // generate new default label for sideset
                std::ostringstream tLabel ;
                if( mMesh->number_of_sidesets() > 9 && tSideSet->id() < 10 )
                {
                    tLabel << "0" ;
                }
                tLabel << tSideSet->id() << "_SideSet";
                tSideSet->label() = tLabel.str() ;
            }

            for( uint s=0; s<3; ++s )
            {
                Cell< DomainGroup * > & tSideSets = s==0 ? mTapes : s==1 ? mInterfaces : mBoundaries ;

                for ( DomainGroup * tGroup : tSideSets )
                {
                    // convert domain type to tag
                    uint tTag = static_cast< uint > ( tGroup->type() );

                    for( id_t tID : tGroup->groups() )
                    {
                        if ( mMesh->sideset_exists( tID ) )
                        {
                            // get block on mesh
                            mesh::SideSet * tSideSet = mMesh->sideset( tID );

                            // loop over all elements of this block
                            Cell< mesh::Facet * > & tFacets = tSideSet->facets();

                            for ( mesh::Facet * tFacet: tFacets )
                            {
                                tFacet->element()->set_physical_tag( tTag );
                            }

                            // generate new label for sideset
                            std::ostringstream tLabel;
                            if ( mMesh->number_of_sidesets() > 9 && tSideSet->id() < 10 )
                            {
                                tLabel << "0";
                            }
                            tLabel << tSideSet->id() << "_";
                            if ( string_to_lower( tGroup->label()) == string_to_lower( to_string( tGroup->type())))
                            {
                                tLabel << to_string( tGroup->type());
                            }
                            else
                            {
                                tLabel << to_string( tGroup->type()) << "_" << tGroup->label();
                            }
                            tSideSet->label() = tLabel.str();
                        }
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::fix_facet_masters()
        {
            // make sure that everything is unflagged
            mMesh->unflag_all_facets();

            // first, we check the interfaces
            for( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                // make sure that this is an internal sideset
                if( tSideSet->facet_by_index( 0 )->has_slave() )
                {
                    // get the sidesets on the mesh
                    Cell< mesh::Facet * > & tFacets = tSideSet->facets() ;

                    // loop over all facets
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        DomainType tMasterType
                                = static_cast< DomainType >( tFacet->master()->physical_tag() );

                        DomainType tSlaveType
                                = static_cast< DomainType >( tFacet->slave()->physical_tag());

                        // the rule is: conductor beats air beats ferro
#ifdef BELFEM_FERRO_HPHIA
                        // check if this sideset is an interface and must be swapped
                        if (   ( tMasterType == DomainType::Air && tSlaveType == DomainType::Conductor )
                            || ( tMasterType == DomainType::Ferro && tSlaveType == DomainType::Conductor )
                            || ( tMasterType == DomainType::Ferro && tSlaveType == DomainType::Air ))
                        {
                            tFacet->flag() ;
                        }
#elif BELFEM_FERROAIR_LAMBDA
                            if (   ( tMasterType == DomainType::Air && tSlaveType == DomainType::Conductor )
                            || ( tMasterType == DomainType::Ferro && tSlaveType == DomainType::Conductor )
                            || ( tMasterType == DomainType::Ferro && tSlaveType == DomainType::Air ))
                        {
                            tFacet->flag() ;
                        }
#else
                        // check if this sideset is an interface and must be swapped
                        if (   ( tMasterType == DomainType::Air && tSlaveType == DomainType::Conductor )
                            || ( tMasterType == DomainType::Ferro && tSlaveType == DomainType::Conductor ) )
                        {
                            tFacet->flag() ;
                        }
#endif
                        else if( tMasterType == tSlaveType )
                        {
                           // otherwise, we prefer the element with the lower ID
                           if( tFacet->master()->id() > tFacet->slave()->id() )
                           {
                               tFacet->flag() ;
                           }
                        }
                    }
                }
            }

            // get the sidesets on the mesh
            Cell< mesh::Facet * > & tFacets = mMesh->facets() ;

            // loop over all facets
            for( mesh::Facet * tFacet : tFacets )
            {
                // check if facet is flagged for swapping
                if( tFacet->is_flagged() )
                {
                    // grab master and slave elements
                    mesh::Element * tMaster = tFacet->master() ;
                    mesh::Element * tSlave  = tFacet->slave() ;

                    // grab indices
                    index_t tMasterIndex = tFacet->master_index() ;
                    index_t tSlaveIndex  = tFacet->slave_index() ;

                    // swap elements
                    tFacet->set_master( tSlave, tSlaveIndex );
                    tFacet->set_slave( tMaster, tMasterIndex );

                    // swap elements
                    tFacet->set_master( tSlave, tSlaveIndex );
                    tFacet->set_slave( tMaster, tMasterIndex );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::create_interfaces()
        {
            Vector< id_t > tSideSetIDs ;
            Vector< id_t > tMasterBlockIDs ;
            Vector< id_t > tSlaveBlockIDs ;
            Vector< uint > tSideSetTypes ;

            // cast domain types
            if( comm_rank() == mMesh->master() )
            {
                // count interfaces
                uint tCount = 0 ;

                // loop over all sidesets
                for( mesh::SideSet * tSideSet : mMesh->sidesets() )
                {
                    // make sure that sideset is not empty
                    if( tSideSet->number_of_facets() > 0 )
                    {
                        // get first facet
                        mesh::Facet * tFacet = tSideSet->facet_by_index( 0 );

                        // make sure that it has a slave
                        if( tFacet->has_slave() )
                        {
                            DomainType tMasterType
                                = static_cast< DomainType >( tFacet->master()->physical_tag() );

                            DomainType tSlaveType
                                = static_cast< DomainType >( tFacet->slave()->physical_tag() );
#ifdef BELFEM_FERRO_HPHIA
                            // check if this sideset is an interface
                            if ( ( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Air )
                               ||( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Ferro )
                               ||( tMasterType == DomainType::Air  && tSlaveType == DomainType::Ferro ) )
                            {
                                ++tCount ;
                            }
#elif BELFEM_FERROAIR_LAMBDA
                            // check if this sideset is an interface
                            if ( ( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Air )
                               ||( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Ferro )
                               ||( tMasterType == DomainType::Air  && tSlaveType == DomainType::Ferro ) )
                            {
                                ++tCount ;
                            }
#else
                            // check if this sideset is an interface
                            if ( ( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Air )
                               ||( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Ferro ) )
                            {
                                ++tCount ;
                            }
#endif
                        }
                    }
                }

                // allocate memory
                tSideSetIDs.set_size( tCount );
                tMasterBlockIDs.set_size( tCount );
                tSlaveBlockIDs.set_size( tCount );
                tSideSetTypes.set_size( tCount );

                // reset counter
                tCount = 0 ;

                // loop over all sidesets
                for( mesh::SideSet * tSideSet : mMesh->sidesets() )
                {
                    // make sure that sideset is not empty
                    if( tSideSet->number_of_facets() > 0 )
                    {
                        // get first facet
                        mesh::Facet * tFacet = tSideSet->facet_by_index( 0 );

                        // make sure that it has a slave
                        if( tFacet->has_slave() )
                        {
                            DomainType tMasterType
                                    = static_cast< DomainType >( tFacet->master()->physical_tag() );

                            DomainType tSlaveType
                                    = static_cast< DomainType >( tFacet->slave()->physical_tag() );

                            // check if this sideset is an interface
                            if ( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Air )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceScAir );
                            }
                            else if ( tMasterType == DomainType::Conductor && tSlaveType == DomainType::Ferro )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceScFm );
                            }
#ifdef BELFEM_FERRO_HPHIA
                            else if ( tMasterType == DomainType::Air && tSlaveType == DomainType::Ferro )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceFmAir );
                            }
#elif BELFEM_FERROAIR_LAMBDA
                            else if ( tMasterType == DomainType::Air && tSlaveType == DomainType::Ferro )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceFmAir );
                            }
#endif
                        }
                    }
                }
                comm_barrier() ;
                send_same( mCommTable, tSideSetIDs );
                send_same( mCommTable, tMasterBlockIDs );
                send_same( mCommTable, tSlaveBlockIDs );
                send_same( mCommTable, tSideSetTypes );
            }
            else
            {
                comm_barrier();
                receive( mMesh->master(), tSideSetIDs );
                receive( mMesh->master(), tMasterBlockIDs );
                receive( mMesh->master(), tSlaveBlockIDs );
                receive( mMesh->master(), tSideSetTypes );
            }

            // create the interfaces
            uint tNumInterfaces = tSideSetIDs.length() ;

            mInterfaces.set_size( tNumInterfaces, nullptr );

            for( uint k=0; k<tNumInterfaces; ++k )
            {
                mInterfaces( k ) = new DomainInterface(
                        static_cast< DomainType >( tSideSetTypes( k ) ),
                        tSideSetIDs( k ),
                        tMasterBlockIDs( k ),
                        tSlaveBlockIDs( k ) );

            }

            if( is_maxwell_ha( mFormulation ) )
            {
                // flag only sc-air interfaces
                for( DomainGroup * tInterface : mInterfaces )
                {
                    if( tInterface->type() == DomainType::InterfaceScAir )
                    {
                        tInterface->flag() ;
                    }
                }
            }
            else if ( is_maxwell_hphi( mFormulation ) )
            {
                // flag all interfaces
                for( DomainGroup * tInterface : mInterfaces )
                {
                    tInterface->flag();
                }
            }

            // add used interfaces to list
            // flag all interfaces
            for( DomainGroup * tInterface : mInterfaces )
            {
                if( tInterface->is_flagged() )
                {
                    mSideSets.push( tInterface );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::create_symmetries( const bool aAntiFlag )
        {
            Vector< id_t > tSideSetIDs ;
            Vector< uint > tMasterTypes ;

            Cell< DomainGroup * > & tRawSymmetries = aAntiFlag ? mRawAntiSymmetries : mRawSymmetries ;
            Cell< DomainGroup * > & tSymmetries = aAntiFlag ? mAntiSymmetries : mSymmetries ;

            const DomainType tAirType = aAntiFlag ?
                    DomainType::AntiSymmetryAir
                    : DomainType::SymmetryAir ;

            const DomainType tFerroType = aAntiFlag ?
                    DomainType::AntiSymmetryFerro :
                    DomainType::SymmetryFerro ;

            const DomainType tConductorType = aAntiFlag ?
                                          DomainType::AntiSymmetryConductor :
                                          DomainType::SymmetryConductor ;

            const string tAirLabel =  aAntiFlag ? "antiSymmetryAir" : "symmetryAir" ;
            const string tFerroLabel =  aAntiFlag ? "antiSymmetryFerro" : "symmetryFerro" ;
            const string tConductorLabel =  aAntiFlag ? "antiSymmetryConductor" : "symmetryConductor" ;

            // cast domain types
            if( comm_rank() == mMesh->master() )
            {
                // reset the counter
                uint tCount = 0 ;

                for( DomainGroup * tGroup : tRawSymmetries )
                {
                    tCount += tGroup->groups().length() ;
                }

                // allocate memory
                tSideSetIDs.set_size( tCount );
                tMasterTypes.set_size( tCount );

                // reset the counter
                tCount = 0 ;

                // loop over all groups
                for( DomainGroup * tGroup : tRawSymmetries )
                {
                    // loop over all sidesets
                    for( id_t s : tGroup->groups() )
                    {
                        mesh::SideSet * tSideSet = mMesh->sideset( s );

                        // get first facet
                        mesh::Facet * tFacet = tSideSet->facet_by_index( 0 );

                        // sanity check
                        BELFEM_ERROR( !tFacet->has_slave(),
                                      "sideset %lu is supposed to be a symmetry, but it is not at the edge of the mesh",
                                      ( long unsigned int ) s );

                        // remember sidesetr ID
                        tSideSetIDs( tCount ) = s;

                        // remember type
                        tMasterTypes( tCount++ ) = tFacet->master()->physical_tag() ;
                    }
                }


                comm_barrier() ;
                send_same( mCommTable, tSideSetIDs );
                send_same( mCommTable, tMasterTypes );
            }
            else
            {
                comm_barrier();
                receive( mMesh->master(), tSideSetIDs );
                receive( mMesh->master(), tMasterTypes );
            }


            uint tNumSideSets = tSideSetIDs.length() ;

            // count sideset types
            uint c = 0 ;
            uint f = 0 ;
            uint a = 0 ;

            for( uint s=0; s<tNumSideSets; ++s )
            {

                switch (  static_cast< DomainType >( tMasterTypes( s ) ) )
                {
                    case( DomainType::Conductor ) :
                    {
                        ++c ;
                        break ;
                    }
                    case( DomainType::Ferro ) :
                    {
                        ++f ;
                        break ;
                    }
                    case( DomainType::Air ) :
                    {
                        ++a ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "This should never happen!");
                    }
                }
            }

            // allocate ids
            Vector< id_t > tAir( a );
            Vector< id_t > tConductor( c ) ;
            Vector< id_t > tFerro( f ) ;

            // reset counters
            a = 0 ;
            f = 0 ;
            c = 0 ;

            // popolate ID vectors
            for( uint s=0; s<tNumSideSets; ++s )
            {
                switch (  static_cast< DomainType >( tMasterTypes( s ) ) )
                {
                    case( DomainType::Conductor ) :
                    {
                        tConductor( c++ ) = tSideSetIDs( s ) ;
                        break ;
                    }
                    case( DomainType::Ferro ) :
                    {
                        tFerro( f++ ) = tSideSetIDs( s ) ;
                        break ;
                    }
                    case( DomainType::Air ) :
                    {
                        tAir( a++ ) = tSideSetIDs( s ) ;
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "This should never happen!");
                    }
                }
            }

            // create sidesets
            if( a > 0 )
            {
                DomainGroup * tGroup = new DomainGroup( tAirType, tAirLabel, tAir );
                mSideSets.push( tGroup );
            }
            if( f > 0 )
            {
                DomainGroup * tGroup = new DomainGroup( tFerroType, tFerroLabel, tFerro );
                mSideSets.push( tGroup );
            }
            if( c > 0 )
            {
                DomainGroup * tGroup = new DomainGroup( tConductorType, tConductorLabel, tConductor );
                mSideSets.push( tGroup );
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::select_bcs_and_cuts()
        {
            if( is_maxwell_ha( mFormulation ) )
            {
                uint tCurrentCount = 0;

                for ( DomainGroup * tGroup: mBlocks )
                {
                    switch ( tGroup->type())
                    {
                        case ( DomainType::Conductor ) :
                        {
                            tCurrentCount += tGroup->groups().length();
                            break;
                        }
                        case ( DomainType::Coil ) :
                        {
                            // check if a current bc exists
                            if ( mCurrentBcMap.key_exists( tGroup->label()))
                            {

                                // link Bc with blocks
                                mCurrentBcMap( tGroup->label())->set_blocks( tGroup->groups());
                            }
                            break;
                        }
                        default :
                        {
                            // pass
                        }
                    }
                }

                mCurrentBlocks.set_size( tCurrentCount );
                tCurrentCount = 0;
                for ( DomainGroup * tGroup: mBlocks )
                {
                    if ( tGroup->type() == DomainType::Conductor )
                    {
                        for ( id_t tID: tGroup->groups())
                        {
                            mCurrentBlocks( tCurrentCount++ ) = tID;
                        }
                    }
                }
            }
            else if (  is_maxwell_hphi( mFormulation ) )
            {
                // count number of air blocks
                index_t tAirCount = 0 ;
                index_t tCurrentCount = 0 ;
                for( DomainGroup * tGroup : mBlocks )
                {
                    switch( tGroup->type() )
                    {
                        case( DomainType::Air ) :
                        {
                            tAirCount += tGroup->groups().length() ;
                            break ;
                        }
                        case( DomainType::Conductor ) :
                        {
                            tCurrentCount += tGroup->groups().length() ;
                            break ;
                        }
                        default :
                        {
                            // pass
                        }
                    }
                }

                // write number of air blocks into temporary vector
                Vector< id_t > tAirBlocks( tAirCount );
                mCurrentBlocks.set_size( tCurrentCount );

                tAirCount = 0 ;
                tCurrentCount = 0 ;

                for( DomainGroup * tGroup : mBlocks )
                {
                    switch( tGroup->type() )
                    {
                        case( DomainType::Air ) :
                        {
                            for ( id_t tID: tGroup->groups( ))
                            {
                                tAirBlocks( tAirCount++ ) = tID;
                            }
                            break;
                        }
                        case( DomainType::Conductor ) :
                        {
                            for ( id_t tID: tGroup->groups() )
                            {
                                mCurrentBlocks( tCurrentCount ++ ) = tID ;
                            }
                            break ;
                        }
                        default :
                        {
                            // pass
                        }
                    }
                }

                mesh::Scissors * tScissors = nullptr ;

                if( comm_rank() == 0 )
                {
                    // create scissors
                    tScissors = new mesh::Scissors( mMesh, tAirBlocks );
                    // add thin sheets
                    for ( DomainGroup * tTape: mTapes )
                    {
                        tScissors->cut( tTape->groups(), tTape->minus(), tTape->plus(), true );
                    }


                    for ( DomainGroup * tCut: mCuts )
                    {
                        tScissors->cut( tCut->groups(), tCut->minus(), tCut->plus());
                    }

                    // finalize cuts
                    tScissors->finalize();
                }

                // add tapes
                if( mTapes.size() > 0 )
                {
                    if( comm_rank() == 0 )
                    {
                        mesh::TapeRoller tTapeRoller( mMesh,
                                                      mTapeMaterialLabels.size(),
                                                      mMesh->max_element_order() );



                        for ( DomainGroup * tTape: mTapes )
                        {
                            tTapeRoller.add_sidesets( tTape->groups());
                            tTapeRoller.add_master_blocks( tTape->master() );
                        }


                        // change the element orientation so that the normals point into the right direction
                        tTapeRoller.flip_element_orientation();

                        mMaxBlockID = tTapeRoller.run();

                        comm_barrier() ;
                        send( mCommTable, mMaxBlockID );

                        // grab the sidesets and compute the normals
                        tTapeRoller.get_sidesets( mThinShellSideSets );
                        mesh::compute_surface_normals( mMesh, mThinShellSideSets, GroupType::SIDESET, false );

                        // revert the element orientation, otherwise the logic is messed up
                        // ( todo: why? I now think that we must not reverse or layers are backwards!)
                        // tTapeRoller.revert_element_orientation();

                        // correct the node positions
                        tTapeRoller.shift_nodes( mTapeThicknesses );

                        // compute the normals
                        if( mMesh->number_of_dimensions() == 2 )
                        {
                            tTapeRoller.compute_edge_signs_2d() ;
                        }

                        // set the labels of the blocks
                        this->set_layer_labels();

                        // remember ids of new sidesets
                        mGhostSideSets = tTapeRoller.ghost_sideset_ids();
                        mGhostBlocks   = tTapeRoller.ghost_block_ids() ;

                        comm_barrier() ;
                        send_same( mCommTable, mThinShellSideSets ) ;
                        send_same( mCommTable, mGhostSideSets );
                        send_same( mCommTable, mGhostBlocks );

                    }
                    else
                    {
                        comm_barrier() ;
                        receive( 0, mMaxBlockID );
                        comm_barrier() ;
                        receive( 0, mThinShellSideSets );
                        receive( 0, mGhostSideSets );
                        receive( 0, mGhostBlocks );
                    }

                }

                comm_barrier();
                // link cut BCs
                for( DomainGroup * tCut : mCuts )
                {
                    // check if a current bc is given for this cut
                    if ( mCurrentBcMap.key_exists( tCut->label()) )
                    {
                        // get bc
                        MaxwellBoundaryConditionCurrent * tBC = mCurrentBcMap( tCut->label() );

                        // check if an associated block exists
                        if( mBlockMap.key_exists( tCut->label() ) )
                        {
                            // link with groups
                            tBC->set_blocks( mBlockMap( tCut->label() )->groups() );
                        }

                        // get the cut IDS
                        Vector< id_t > tCutIDs ;
                        if( comm_rank() == 0 )
                        {
                            // grab IDs from scissors
                            tScissors->get_cut_ids( tCut->groups(), tCutIDs );

                            // distribute ids
                            send_same( mCommTable, tCutIDs );
                        }
                        else
                        {
                            receive( 0, tCutIDs );
                        }

                        // link IDs to BC
                        tBC->set_sidesets( tCutIDs );

                        // add bc to list
                        mSideSets.push( tCut );
                    }
                }

                // wait for all procs
                comm_barrier();

                // create the cut map that associates sidesets with cuts


                if( gComm.rank() == 0 )
                {
                    // refresh cut map
                    tScissors->collect_cut_data();

                    const Matrix< id_t > & tCutData = tScissors->cut_data() ;

                    send_same( mCommTable, tScissors->cut_data() );

                    mSideSetToCutMap.clear();

                    uint tN = tCutData.n_cols() ;
                    for( uint k=0; k<tN; ++k )
                    {
                        mSideSetToCutMap[ tCutData( 1, k ) ] = tCutData( 0, k );
                    }

                    delete tScissors ;
                }
                else
                {
                    Matrix< id_t > tCutData ;
                    receive( 0 , tCutData );

                    // crate the map
                    mSideSetToCutMap.clear();

                    uint tN = tCutData.n_cols() ;
                    for( uint k=0; k<tN; ++k )
                    {
                        mSideSetToCutMap[ tCutData( 1, k ) ] = tCutData( 0, k );
                    }

                }

                this->remove_coils();
            } // end h-phi formulation
            else
            {
                BELFEM_ERROR( false, "Unsupported IWG type") ;
            }

            // set boundaries
            for(  DomainGroup * tGroup : mBoundaries )
            {
                if( mMagfieldBcMap.key_exists( tGroup->label() ) )
                {
                    mMagfieldBcMap( tGroup->label() )->set_sidesets( tGroup->groups() );
                    mSideSets.push( tGroup );
                }
            }

            // add thin shells
            for( DomainGroup * tGroup : mTapes )
            {
                mSideSets.push( tGroup );
            }

            // create the comm table

            // wait for all procs
            comm_barrier();
        }


// -----------------------------------------------------------------------------

        void
        MaxwellFactory::select_sidesets()
        {
            // counter for sidesets
            index_t tCount = 0 ;

            // create typemap
            Map< id_t, DomainType > tSideSetTypeMap ;

            // loop over all groups
            for( DomainGroup * tGroup : mSideSets )
            {
                // add length to counter
                tCount += tGroup->groups().length() ;
            }

            // allocate memory
            mSideSetIDs.set_size( tCount );

            // reset counter
            tCount = 0 ;

            // collect group IDs
            for( DomainGroup * tGroup : mSideSets )
            {
                for( id_t tID : tGroup->groups() )
                {
                    if ( !tSideSetTypeMap.key_exists( tID ) )
                    {
                        tSideSetTypeMap[ tID ] = tGroup->type();

                    }
                    mSideSetIDs( tCount++ ) = tID ;
                }
            }

            // make sure that everything is unique
            unique( mSideSetIDs );

            BELFEM_ERROR( mSideSetIDs.length() == tCount, "List of SideSets or Cuts is not unique" );

            // populate type container
            mSideSetTypes.set_size( tCount, DomainType::UNDEFINED );
            tCount = 0 ;
            for( id_t tID : mSideSetIDs )
            {
                mSideSetTypes( tCount++ ) = tSideSetTypeMap( tID );
            }
        }

// -----------------------------------------------------------------------------

        Spline *
        MaxwellFactory::read_bhfile(
                const string & aPath,
                const string & aUnitB,
                const string & aUnitH,
                const value    aMaxB,
                      real   & aM )
        {

            if( comm_rank() == 0 )
            {
                // - - - - - - - - - - - - - - - - - - -
                // Step 0 : check unuts
                // - - - - - - - - - - - - - - - - - - -
                value tUnitB = unit_to_si( aUnitB );
                value tUnitH = unit_to_si( aUnitH );
                BELFEM_ERROR( check_unit( tUnitB, "T" ), "Wrong unit for b-field : %s", aUnitB.c_str());
                BELFEM_ERROR( check_unit( tUnitH, "A/m" ), "Wrong unit for h-field : %s", aUnitB.c_str());

                // - - - - - - - - - - - - - - - - - - -
                // Step 1 : read textfile
                // - - - - - - - - - - - - - - - - - - -

                // read file
                Ascii tFile( aPath, FileMode::OPEN_RDONLY );

                // check if first line contains an index
                Cell <string> tWords = string_to_words( clean_string( tFile.line( 0 )));

                // override input settings
                if( tWords.size() == 3 )
                {
                    tUnitB = string_to_bool( tWords( 1 ) ) ? unit_to_si("T") : unit_to_si("G");
                    tUnitH = string_to_bool( tWords( 2 ) ) ? unit_to_si("A/m") : unit_to_si("Oe");
                }

                uint l0 = tWords.size() == 2 ? 0 : 1;
                uint ln = tWords.size() == 2 ? tFile.length() : std::stoi( tWords( 0 )) + 1;

                // create vectors with values
                uint tn0 = ln - l0 ;
                Vector <real> tB0( tn0, 0.0 );
                Vector <real> tH0( tn0, 0.0 );

                // read data
                uint tCount = 0;
                for ( index_t l = l0; l < ln; ++l )
                {
                    Cell <string> tBH = string_to_words( clean_string( tFile.line( l )));
                    tB0( tCount ) = std::stod( tBH( 0 )) * tUnitB.first;
                    tH0( tCount++ ) = std::stod( tBH( 1 )) * tUnitH.first;
                }

                BELFEM_ERROR( tB0( 0 ) == 0.0, "b-datapoints must begin with zero" );
                BELFEM_ERROR( tn0 > 2, "need at least three datapoints!" );

                // offset
                aM = tH0( 0 );
                tH0 -= aM ;

                // check inputs

                real tBn = tB0( tn0 - 1 );
                real tBmax = aMaxB.first;

                BELFEM_ERROR( tBmax > tBn, "max value of b mus be bigger than in textfile!" );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2 : compute derivative at final point
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -


                real tdHdBn  ;
                real td2HdB2n  ;
                real tHn ;

                // number of samples for poly
                uint tn = tn0 < 5 ? tn0 : 5 ;

                Vector< real > tB1( tn );
                Vector< real > tH1( tn );
                Vector< real > tC( 3 );

                for( uint k=0; k<tn; ++k )
                {
                    tB1( k ) = tB0( tn0 - tn + k );
                    tH1( k ) = tH0( tn0 - tn + k );
                }
                polyfit( tB1, tH1, 2, tC );
                tdHdBn = dpolyval( tC, tBn );
                td2HdB2n = 2 * tC( 0 );
                tHn = polyval( tC, tBn );


                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3 : extrapolate curve
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // create the ode
                BhExtrapolationOde tOde( tBn, tdHdBn, td2HdB2n, tBmax );

                // create the integrator
                ode::Integrator tIntegrator( tOde, ode::Type::RK45 );

                real tX = tBn ;
                tIntegrator.timestep() = tB0( tn0-1 ) - tB0( tn0 -2 ) ;
                tIntegrator.maxtime() = tBmax ;
                Vector< real > tY( 1 , tHn );

                uint tn2 = 0 ;
                Cell< real > tB3 ;
                Cell< real > tH3 ;
                while ( ( uint ) tIntegrator.step( tX, tY ) < 2 )
                {
                    tB3.push( tX );
                    tH3.push( tY( 0 ) );
                    if( tX > tOde.bmax() )
                    {
                        break ;
                    }
                }
                tBmax = tX ;

                uint tn1 = tB3.size() ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4 : merge datasets
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                Vector< real > tB2( tn0 + tn1 );
                Vector< real > tH2( tn0 + tn1 );


                // first dataset
                for( uint k=1; k<tn0; ++k )
                {
                    tB2( k ) = tB0( k );
                    tH2( k ) = tH0( k );
                }

                tCount = tn0 ;

                // second dataset
                for( uint k=0; k<tn1; ++k )
                {
                    tB2( tCount ) = tB3( k );
                    tH2( tCount ) =tH3( k ) ;
                    ++tCount ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5 : project data to an equidistant grid
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                uint tNumPoints = 1001;
                Vector <real> tB( tNumPoints );
                Vector <real> tH( tNumPoints );

                linspace( 0.0, tBmax, tNumPoints, tB );
                OneDMapper tMapper( tB, 1 );
                tMapper.project( tB2, tH2, tH );
                tY.set_size( tNumPoints );
                for( uint k=1; k<tNumPoints; ++k )
                {
                    tY( k ) = std::log( tH( k ) / tB( k ) );
                }


                // fix peak
                real tBmin = 0.2 ;

                for( tCount = 0 ; tB( tCount ) < tBmin  ; ++tCount );

                // sign :
                real tS0 = tY( tCount + 1 ) - tY( tCount );

                bool tFoundPeak = false ;

                // look for peak ( a numerical bug )
                for( int k=tCount; k>0; --k )
                {
                    real tS = tY( k + 1 ) - tY( k ) ;
                    if( tS * tS0 < 0 )
                    {
                        tFoundPeak = true ;
                        tCount = k ;
                        break ;
                    }
                }

                if( ! tFoundPeak )
                {
                    tCount = 1 ;
                    tFoundPeak = true ;
                }

                // check if we have found a peek
                if( tFoundPeak )
                {
                    // derivative at reference point
                    real a = ( tY( tCount+2 ) - tY( tCount+1 ) ) / ( tB( tCount + 2 ) - tB( tCount + 1 ) );
                    real b = tY( tCount + 1 ) - a * tB( tCount + 1 );

                    // fix peak
                    for( uint k=0; k<=tCount; ++k )
                    {
                        tY( k ) = a * tB( k ) + b ;
                    }
                }

                // - - - - - - - - - -   - - - - - - - - - - - - - - - - - -
                // Step 7 : synchronize M0 value
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                proc_t tCommSize = gComm.size() ;

                if(tCommSize >  1 )
                {
                    // create commlist
                    Vector< proc_t > tCommList ( tCommSize );

                    for( proc_t p=0; p<tCommSize; ++p )
                    {
                        tCommList( p ) = p ;
                    }
                    Vector< real > tAllM( tCommSize, aM );

                    send( tCommList, tAllM );
                }


                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 7 : create the spline
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                tn2 = tn0 + tn1 ;

                // help matrix for spline
                SpMatrix tHelpMatrix;
                spline::create_helpmatrix( tNumPoints, tB( 1 ), tHelpMatrix );

                Cell< Spline * > aSplines( 2 , nullptr );
                return new Spline( tB, tY, tHelpMatrix, 0, 0, 0  );

            }
            else
            {
                receive( 0, aM );
                return new Spline( 0 );
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::remove_coils()
        {
            // remove coils from lists
            Cell< DomainType > tBlockTypes = mBlockTypes ;
            Vector< id_t > tBlockIDs = mBlockIDs ;

            // count non coil types
            uint tCount = 0 ;

            for( DomainType tType : mBlockTypes )
            {
                if( tType != DomainType::Coil )
                {
                    ++tCount ;
                }
            }

            mBlockIDs.set_size( tCount );
            mBlockTypes.set_size( tCount, DomainType::UNDEFINED );
            tCount = 0 ;


            for( uint b=0; b<tBlockIDs.length(); ++b )
            {
                if( tBlockTypes( b ) != DomainType::Coil )
                {
                    mBlockIDs( tCount ) = tBlockIDs( b );
                    mBlockTypes( tCount++ ) = tBlockTypes( b );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::set_layer_labels()
        {
            // set the names of the sidesets
            unsigned int tNumLayers = mTapeThicknesses.length();

            for( unsigned int k=0; k<tNumLayers; ++k )
            {
                // grab Block
                mesh::Block * tBlock = mMesh->block( k + mMaxBlockID + 1 );

                // create layer name
                tBlock->label() =
                sprint( "Layer_%u_%s_%uu",
                         k + 1,
                         mTapeMaterialLabels( k  ).c_str() ,
                         ( unsigned int ) std::ceil( mTapeThicknesses( k ) *1e6 ) );

            }
        }

// -----------------------------------------------------------------------------

        void
        MaxwellFactory::read_id_line( const string & aLine, Vector< id_t > & aIDs )
        {
            // tidy up and create words
            Cell< string > tWords = string_to_words(
                    search_and_replace(
                            search_and_replace( aLine, ",", " " ),
                            ";", " " ));

            // count numbers
            uint tCount = 0 ;
            uint tNumWords = tWords.size() ;

            // count numbers
            for ( uint w=0; w<tNumWords; ++w )
            {
                if( tWords( w ) == ":" )
                {
                    tCount += std::stoi( tWords( w+1 ) ) - std::stoi( tWords( w-1 ) ) - 1 ;
                }
                else
                {
                    ++tCount ;
                }
            }

            // allocate vector
            aIDs.set_size( tCount );
            tCount = 0 ;
            for ( uint w=0; w<tNumWords; ++w )
            {
                if( tWords( w ) == ":" )
                {
                    id_t tA = std::stoi( tWords( w-1 ) ) + 1 ;
                    id_t tB = std::stoi( tWords( w+1 ) ) ;
                    for( id_t tX = tA ; tX < tB; ++tX )
                    {
                        aIDs( tCount++ ) = tX ;
                    }
                }
                else
                {
                    aIDs( tCount++ ) = std::stoi( tWords( w ) );
                }
            }
        }
// -----------------------------------------------------------------------------

        void
        MaxwellFactory::reduce_order_of_ferro_blocks_tri6()
        {
            // first we need to flag all nodes that are meant to be deleted
            mMesh->unflag_all_nodes() ;
            mMesh->unflag_all_elements() ;

            uint tNumBlocks   = mBlockIDs.length() ;

            mesh::ElementFactory tFactory ;

            // count blocks with new elements that are to be created
            index_t tCount = 0 ;
            for( uint b=0; b<tNumBlocks; ++b )
            {
                // check if this is a ferro block
                if ( mBlockTypes( b ) == DomainType::Ferro )
                {
                    tCount += mMesh->block( mBlockIDs( b ) )->number_of_elements() ;
                }
            }

            // allocate container with elements that we want to delete
            Cell< mesh::Element * > tOldElements( tCount, nullptr );

            // reset the counter
            tCount = 0 ;

            // we also grab all elements form the mesh
            Cell< mesh::Element * > & tMeshElements = mMesh->elements() ;

            for( uint b=0; b<tNumBlocks; ++b )
            {
                // check if this is a ferro block
                if ( mBlockTypes( b ) == DomainType::Ferro )
                {
                    // get the elements
                    Cell< mesh::Element * > & tBlockElements
                            = mMesh->block( mBlockIDs( b ) )->elements() ;

                    // get the number of elements
                    index_t tNumElems = tBlockElements.size() ;

                    for( index_t e=0; e<tNumElems; ++e )
                    {
                        // grab the element
                        mesh::Element * tOldElement = tBlockElements( e );
                        tOldElement->flag();

                        // flag node for deletion
                        for( uint k=3; k<6;  ++k )
                        {
                            tOldElement->node( k )->flag() ;
                        }

                        // create a new element
                        mesh::Element * tNewElement = tFactory.create_lagrange_element( ElementType::TRI3, tOldElement->id() );

                        // copy some data
                        tNewElement->set_block_id( tOldElement->block_id() );
                        tNewElement->set_geometry_tag( tOldElement->geometry_tag() );
                        tNewElement->set_physical_tag( tOldElement->physical_tag() );

                        tNewElement->set_index( tOldElement->index() );

                        // link nodes
                        for( uint k=0; k<3;  ++k )
                        {
                            tNewElement->insert_node( tOldElement->node( k ), k );
                        }

                        // replace old element in containers
                        tMeshElements( tOldElement->index() ) = tNewElement ;
                        tBlockElements( e ) = tNewElement ;

                        // add old element to list for deletion
                        tOldElements( tCount++ ) = tOldElement ;

                    } // end loop over all elements

                } // end ferro block
            } // end loop over all blocks

            index_t tNumFacets = mMesh->number_of_facets() ;
            Cell< mesh::Facet * > & tMeshFacets = mMesh->facets() ;

            for( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                Cell< mesh::Facet * > & tSideSetFacets = tSideSet->facets();
                index_t tNumFacets = tSideSetFacets.size();

                for ( index_t f = 0; f < tNumFacets; ++f )
                {
                    mesh::Facet * tOldFacet = tSideSetFacets( f );

                    // switch if we delete the midnode
                    bool tDeleteMidnode = false;

                    // switch if we realign the midnode
                    bool tChangeMaster = false;

                    // switch if we realign the midnode
                    bool tChangeSlave = false;

                    // check if facet has a master
                    if ( tOldFacet->has_master() )
                    {
                        // check if master is ferro
                        if ( tOldFacet->master()->is_flagged())
                        {
                            // we change the master
                            tChangeMaster = true;

                            // check if facet has slave
                            if ( tOldFacet->has_slave() )
                            {
                                // check if facet is ferro
                                if ( tOldFacet->slave()->is_flagged())
                                {
                                    // we delete the midnode
                                    tDeleteMidnode = true;
                                }
                            }
                            else
                            {
                                // we delete the midnode
                                tDeleteMidnode = true;
                            }
                        }
                    }

                    if ( tOldFacet->has_slave() )
                    {
                        // check if facet is ferro
                        if ( tOldFacet->slave()->is_flagged() )
                        {
                            // we change  the slave
                            tChangeSlave = true;

                            if( ! tOldFacet->has_master() )
                            {
                                tDeleteMidnode = true ;
                            }
                        }
                    }

                    if ( tDeleteMidnode )
                    {
                        // create a new facet
                        mesh::Element * tElement = tFactory.create_lagrange_element( ElementType::LINE2,
                                                                                     tOldFacet->id());
                        tElement->insert_node( tOldFacet->element()->node( 0 ), 0 );
                        tElement->insert_node( tOldFacet->element()->node( 1 ), 1 );

                        mesh::Facet * tNewFacet = new mesh::Facet( tElement );
                        tNewFacet->set_index( tOldFacet->index());
                        if ( tOldFacet->has_master())
                        {
                            tNewFacet->set_master( tMeshElements( tOldFacet->master()->index()),
                                                   tOldFacet->master_index(), false );
                        }
                        if ( tOldFacet->has_slave())
                        {
                            tNewFacet->set_slave( tMeshElements( tOldFacet->slave()->index()), tOldFacet->slave_index());
                        }

                        // replace facet in containers
                        tSideSetFacets( f ) = tNewFacet ;
                        tMeshFacets( tOldFacet->index() ) = tNewFacet;
                        delete tOldFacet;
                    }
                    else
                    {
                        if ( tChangeMaster )
                        {
                            tOldFacet->set_master( tMeshElements( tOldFacet->master()->index() ), tOldFacet->master_index(),
                                                   false );
                        }
                        if ( tChangeSlave )
                        {
                            tOldFacet->set_slave( tMeshElements( tOldFacet->slave()->index()), tOldFacet->slave_index());
                        }

                        if ( tChangeMaster || tChangeSlave )
                        {
                            // center this node
                            real tX = 0.5 * ( tOldFacet->node( 0 )->x() + tOldFacet->node( 1 )->x());
                            real tY = 0.5 * ( tOldFacet->node( 0 )->y() + tOldFacet->node( 1 )->y());
                            tOldFacet->node( 2 )->set_coords( tX, tY );

                            // we don't want to delete this node anymore
                            tOldFacet->node( 2 )->unflag();
                        }
                    }
                } // end loop over all sidesets
            } // end loop over all facets

            // delete elements
            for( mesh::Element * tElement : tOldElements )
            {
                delete tElement ;
            }

            // delete nodes
            Cell< mesh::Node * > & tNewNodes = mMesh->nodes() ;
            Cell< mesh::Node * > tOldNodes;
            tOldNodes.vector_data() =  std::move( tNewNodes.vector_data() );

            tCount = 0 ;
            for( mesh::Node * tNode : tOldNodes )
            {
                if ( ! tNode->is_flagged() )
                {
                    ++tCount ;
                }
            }
            tNewNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for( mesh::Node * tNode : tOldNodes )
            {
                if ( tNode->is_flagged() )
                {
                    delete tNode ;
                }
                else
                {
                    tNewNodes( tCount++ ) = tNode ;
                }
            }


            // recreate mesh
            mMesh->unfinalize() ;
            mMesh->finalize() ;
        }

// -----------------------------------------------------------------------------
    }
}