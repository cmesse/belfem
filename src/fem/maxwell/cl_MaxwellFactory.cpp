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


#include "cl_IWG_Maxwell_HA_Tri3.hpp"
#include "cl_IWG_Maxwell_HA_Tri6.hpp"

#include "cl_IWG_Maxwell_HPhi_Tri3.hpp"
#include "cl_IWG_Maxwell_HPhi_Tri6.hpp"
#include "cl_IWG_Maxwell_HPhi_Tet4.hpp"
#include "cl_IWG_Maxwell_HPhi_Tet10.hpp"

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

            // create the edges for this mesh
            mMesh->create_edges( false, mNedelecBlocks, mNedelecSideSets ) ;


            // check if this is a higher order mesh
            if( mMesh->max_element_order() > 1 )
            {
                // create faces
                mMesh->create_faces( false, mNedelecBlocks, mNedelecSideSets ) ;
            }

            this->read_formulation();

            // read the timestepping information
            this->read_timestep() ;

            this->create_bcs() ;

            // define which blocks and sidesets are used and what they
            // represent
            this->create_topology();

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

            Mesh * aMesh = new Mesh( aInputFile.section("mesh")->get_string("meshFile")  );

            // get mesh unit
            string tUnit = aInputFile.section("mesh")->get_string("meshUnit");

            // we only need to scale the mesh if it is not in meters
            if( tUnit != "m" )
            {
                // grab value
                value tValue = unit_to_si( tUnit );

                BELFEM_ERROR( check_unit( tValue, "m" ), "Invalid length unit: %s", tUnit.c_str() );

                aMesh->scale_mesh( tValue.first );
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
                    }
                    else if ( tFlag == 4 )
                    {
                        mAlphaB = to_real( tWord );
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

                // check if resistance exists
                if ( aInput->key_exists( "resistance" ) )
                {
                    // check type of resistance law
                    const string tType = string_to_lower( aInput->get_string( "resistance" ));

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

                        BELFEM_ERROR(
                                ( aInput->key_exists( "jc" ) && !aInput->key_exists( "jc0" ))
                                || ( !aInput->key_exists( "jc" ) && aInput->key_exists( "jc0" )),
                                "variable must be called either jc or jc0, not both" );

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
                    else if ( tType == "powerlaw-ej" )
                    {
                        real tJc = aInput->key_exists( "jc0" ) ?
                                   aInput->get_value( "jc0", "A/m^2" ).first :
                                   aInput->get_value( "jc", "A/m^2" ).first;

                        real tN = aInput->key_exists( "n0" ) ?
                                  aInput->get_value( "n0", "-" ).first :
                                  aInput->get_value( "n", "-" ).first;

                        real tEc = aInput->key_exists( "rho_c" ) ?
                                   aInput->get_value( "rho_c", "Ohm*m" ).first * tJc :
                                   aInput->get_value( "ec", "V/m" ).first;

                        aMaterial->set_rho_el_ej( tEc, tJc, tN );

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
                        aMaterial->set_nu_s( this->read_bhfile( tPath, tUnitB, tUnitH, tMaxB, tM ));
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
                if(    tGroup->type() == DomainType::SuperConductor
                    || tGroup->type() == DomainType::FerroMagnetic )
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
            const string tTypeString = tString.substr( 0, tSplit );


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

            // tidy up and create words
            Cell< string > tWords = string_to_words(
                    clean_string(
                            search_and_replace(
                                    search_and_replace( tLine, ",", " " ),
                                    ";", " " )));


            // convert to vector
            Vector< id_t > tData( tWords.size());

            for ( index_t k = 0; k < tWords.size(); ++k )
            {
                tData( k ) = std::stoi( tWords( k ));
            }

            // end get raw data
            // -  - - - - - - - - - - - - - - -

            // create domain type
            DomainType tType = domain_type( tTypeString );

            switch ( tType )
            {
                case ( DomainType::SuperConductor ) :
                case ( DomainType::FerroMagnetic )  :
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
                    BELFEM_ERROR( tData.length() % 4 == 0, "Number of cut data must be divisible by four" );

                    index_t tN = tData.length() / 4 ;

                    // populate data
                    Vector< id_t > tGroups( tN );
                    Vector< id_t > tPlus( tN );
                    Vector< id_t > tMinus( tN );
                    Vector< id_t > tMaster( tN );

                    index_t tCount = 0;

                    for ( uint i = 0; i < tN; ++i )
                    {
                        tGroups( i ) = tData( tCount++ );
                        tPlus( i )   = tData( tCount++ );
                        tMinus( i )  = tData( tCount++ );
                        tMaster( i ) = tData( tCount++ );
                    }

                    // create the group object
                    return new DomainThinShell( tLabel, tGroups,
                                                tPlus, tMinus, tMaster, tType  );
                }
                case ( DomainType::Air ) :
                case ( DomainType::Boundary ) :
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
                if( tScalingType == "farfield"  )
                {
                    tMagfieldBcType = MagfieldBcType::Farfied ;
                }
                else if( tScalingType == "symmetry"  )
                {
                    tMagfieldBcType = MagfieldBcType::Symmetry ;
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
                        BoundaryConditionImposing tImposing =  tMagfieldBcType == MagfieldBcType::Farfied ?
                                                                 BoundaryConditionImposing::Dirichlet :
                                                                 BoundaryConditionImposing::Weak ;

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
                }
                else if ( comm_rank() == 0 )
                {
                    std::cout << "Warning: bc type " << tWords( 0 )  << " not implemented"  << std::endl ;
                }
            }


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

            //mParameters->set_sideset_integration_orders( 5 );
            //mParameters->set_block_integration_orders( 5 );
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
                case( IwgType::MAXWELL_HA_TRI3 ) :
                {
                    return new IWG_Maxwell_HA_Tri3() ;
                }
                case( IwgType::MAXWELL_HA_TRI6 ) :
                {
                    return new IWG_Maxwell_HA_Tri6() ;
                }
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
                case( IwgType::MAXWELL_HPHI_TET4 ) :
                {
                    return new IWG_Maxwell_HPhi_Tet4() ;
                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TET10 ) :
                {
                    return new IWG_Maxwell_HPhi_Tet10() ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "invalid type");
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
            // this is a hack, we will remove the old maxwell things
            IWG_Maxwell * tMaxwell = this->create_iwg( mFormulation );

            // Dirichlet bcs do not need elements. We remove those
            this->blacklist_dirichlet_sidesets();

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
            mMagneticField->create_fields( tMaxwell );

            // link field with materials
            this->link_materials( mMagneticField );

            // set the tape materials, if tapes exist
            this->set_tape_materials( mMagneticField );

            // set timestepping info
            tMaxwell->set_timestep( mTimeStep );

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
                if( mBlockTypes( b ) == DomainType::SuperConductor )
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
                    if ( mBlockTypes( b ) == DomainType::SuperConductor )
                    {
                        tBlocks( tCount++ ) = mBlockIDs( b );
                    }
                }

                Cell< DomainType > tTypes( tCount, DomainType::SuperConductor );

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
                    case( DomainType::FerroMagnetic ) :
                    {
                        ++tFerroCount ;
                        break ;
                    }
                    case( DomainType::SuperConductor ) :
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
                    case( DomainType::FerroMagnetic ) :
                    {
                        tFerroBlocks( tFerroCount++ ) =   mBlockIDs( b );
                        break ;
                    }
                    case( DomainType::SuperConductor ) :
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
                this->add_postprocessor_to_kernel( tIWG );

                // fix group types of IWG
                tIWG->update_group_types();
            }

            // crate ferro projector
            if( tFerroCount > 0 )
            {
                // create the iwg
                IWG_Maxwell * tIWG = new IWG_Maxwell_L2_Magfield(
                        this->element_type( tFerroBlocks( 0 )),
                        tMode,
                        mAlphaB );

                Cell< DomainType > tTypes( tFerroBlocks.length(), DomainType::FerroMagnetic );

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

                Cell< DomainType > tTypes( tConductorBlocks.length(), DomainType::SuperConductor );

                // set the blocks of the IWG
                tIWG->set_blocks( tConductorBlocks, tTypes );

                // add iwg to kernel
                this->add_postprocessor_to_kernel( tIWG );

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

        void
        MaxwellFactory::add_postprocessor_to_kernel( IWG_Maxwell * aEquation )
        {
            // add equation to kernel
            mKernel->add_equation( aEquation );

            // create dof projector
            DofManager * tProjector = mKernel->create_field( aEquation );

            // link field with materials
            this->link_materials( tProjector );

            // make sure that data arrays on mesh exist
            tProjector->create_fields( aEquation );

            // set the solver of this field based on the input
            this->set_solver( tProjector,
                              mInputFile.section("maxwell")
                                      ->section("solver")->section("linear") );

            // no BLR for L2, otherwise crash!
            tProjector->solver()->set_mumps_blr( BlockLowRanking::Off, 0.0 );

            tProjector->fix_node_dofs_on_ghost_sidesets() ;

            // add projector to list
            mMagneticField->postprocessors().push( tProjector );
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
                this->write_sideset_phystags_to_mesh();
            }
            else
            {
                this->create_interfaces();
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
                        if( tGroup->type() == DomainType::SuperConductor )
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
                if( tBlockTypeMap( tID ) == DomainType::SuperConductor )
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

            for( mesh::Block * tBlock : mMesh->blocks() )
            {
                // generate new default label for sideset
                std::ostringstream tLabel ;
                if( mMesh->number_of_blocks() > 9 && tBlock->id() < 10 )
                {
                    tLabel << "0" ;
                }
                tLabel << tBlock->id() << "_Block";
                tBlock->label() = tLabel.str() ;
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

                        // generate new label for block
                        std::ostringstream tLabel;
                        if ( mMesh->number_of_blocks() > 9 && tBlockID < 10 )
                        {
                            tLabel << "0";
                        }
                        tLabel << tBlock->id() << "_";
                        if ( string_to_lower( tGroup->label()) == string_to_lower( to_string( tGroup->type())))
                        {
                            tLabel << to_string( tGroup->type());
                        }
                        else
                        {
                            tLabel << to_string( tGroup->type()) << "_" << tGroup->label();
                        }
                        tBlock->label() = tLabel.str();
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
            // get the sidesets on the mesh
            Cell< mesh::Facet * > & tFacets = mMesh->facets() ;

            // loop over all facets
            for( mesh::Facet * tFacet : tFacets )
            {
                // check if facet has a slave
                if( tFacet->has_slave() )
                {
                    // grab master and slave elements
                    mesh::Element * tMaster = tFacet->master() ;
                    mesh::Element * tSlave  = tFacet->slave() ;

                    // if master tag is smaller than slave, we must swap
                    if( tSlave->physical_tag() > tMaster->physical_tag() )
                    {
                        // grab indices
                        index_t tMasterIndex = tFacet->master_index() ;
                        index_t tSlaveIndex  = tFacet->slave_index() ;

                        // swap elements
                        tFacet->set_master( tSlave, tSlaveIndex );
                        tFacet->set_slave( tMaster, tMasterIndex );

                    }
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

                            // check if this sideset is an interface
                            if ( ( tMasterType == DomainType::SuperConductor && tSlaveType == DomainType::Air )
                               ||( tMasterType == DomainType::SuperConductor && tSlaveType == DomainType::FerroMagnetic )
                               ||( tMasterType == DomainType::FerroMagnetic  && tSlaveType == DomainType::Air ) )
                            {
                                ++tCount ;
                            }
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
                            if ( tMasterType == DomainType::SuperConductor && tSlaveType == DomainType::Air )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceScAir );
                            }
                            else if ( tMasterType == DomainType::SuperConductor && tSlaveType == DomainType::FerroMagnetic )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceScFm );
                            }
                            else if ( tMasterType == DomainType::FerroMagnetic && tSlaveType == DomainType::Air )
                            {
                                tSideSetIDs( tCount )     = tSideSet->id();
                                tMasterBlockIDs( tCount ) = tFacet->master()->geometry_tag() ;
                                tSlaveBlockIDs( tCount )  = tFacet->slave()->geometry_tag() ;
                                tSideSetTypes( tCount++ ) = static_cast< uint >( DomainType::InterfaceFmAir );
                            }
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
        MaxwellFactory::select_bcs_and_cuts()
        {
            if( is_maxwell_ha( mFormulation ) )
            {
                uint tCurrentCount = 0;

                for ( DomainGroup * tGroup: mBlocks )
                {
                    switch ( tGroup->type())
                    {
                        case ( DomainType::SuperConductor ) :
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
                    if ( tGroup->type() == DomainType::SuperConductor )
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
                        case( DomainType::SuperConductor ) :
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
                        case( DomainType::SuperConductor ) :
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

                // create scissors
                mesh::Scissors tScissors( mMesh, tAirBlocks );
                // add thin sheets
                for( DomainGroup * tTape : mTapes )
                {
                    tScissors.cut( tTape->groups(), tTape->minus(), tTape->plus(), true );
                }



                for( DomainGroup * tCut : mCuts )
                {
                    tScissors.cut( tCut->groups(), tCut->minus(), tCut->plus() );
                }

                // finalize cuts
                tScissors.finalize() ;

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

                        mMaxBlockID = tTapeRoller.run();

                        comm_barrier() ;
                        send( mCommTable, mMaxBlockID );

                        // change the element orientation so that the normals point into the right direction
                        tTapeRoller.flip_element_orientation();

                        // grab the sidesets and compute the normals
                        Vector< id_t > tSideSets;
                        tTapeRoller.get_sidesets( tSideSets );
                        mesh::compute_surface_normals( mMesh, tSideSets, GroupType::SIDESET, false );

                        // revert the element orientation, otherwise the logic is messed up
                        tTapeRoller.revert_element_orientation();

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
                        send_same( mCommTable, mGhostSideSets );
                        send_same( mCommTable, mGhostBlocks );

                    }
                    else
                    {
                        comm_barrier() ;
                        receive( 0, mMaxBlockID );
                        comm_barrier() ;
                        receive( 0, mGhostSideSets );
                        receive( 0, mGhostBlocks );
                    }
                }

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

                        // grab IDs from scissors
                        tScissors.get_cut_ids( tCut->groups(), tCutIDs );


                        // link IDs to BC
                        tBC->set_sidesets( tCutIDs );

                        // add bc to list
                        mSideSets.push( tCut );
                    }
                }

                // wait for all procs
                comm_barrier();

                // create the cut map that associates sidesets with cuts
                Vector< id_t > tCutIDs ;
                Vector< id_t > tSideSetIDs ;
                if( gComm.rank() == 0 )
                {
                    Vector< proc_t > tCommTable( gComm.size() );
                    for( proc_t p=0; p<gComm.size(); ++p )
                    {
                        tCommTable( p ) = p ;
                    }

                    // refresh cut map
                    tScissors.collect_cut_data();

                    tCutIDs   = trans( tScissors.cut_data().row( 0 ) );
                    tSideSetIDs = trans( tScissors.cut_data().row( 1 ) );

                    send_same( tCommTable, tCutIDs );
                    send_same( tCommTable, tSideSetIDs );
                }
                else
                {
                    receive( 0 , tCutIDs );
                    receive( 0 , tSideSetIDs );
                }

                // crate the map
                mSideSetToCutMap.clear();
                for( uint k=0; k<tCutIDs.length(); ++k )
                {
                    mSideSetToCutMap[ tSideSetIDs( k ) ] = tCutIDs( k );
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
                Vector <real> tB0( ln - l0, 0.0 );
                Vector <real> tH0( ln - l0, 0.0 );

                // read data
                uint tCount = 0;
                for ( index_t l = l0; l < ln; ++l )
                {
                    Cell <string> tBH = string_to_words( clean_string( tFile.line( l )));
                    tB0( tCount ) = std::stod( tBH( 0 )) * tUnitB.first;
                    tH0( tCount++ ) = std::stod( tBH( 1 )) * tUnitH.first;
                }

                BELFEM_ERROR( tB0( 0 ) == 0.0, "b-datapoints must begin with zero" );

                // offset
                aM = tH0( 0 );
                tH0 -= aM ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2 : project data to finer and equidistant grid
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                uint tNumPoints0 = 100;
                Vector <real> tB1( tNumPoints0 );
                Vector <real> tH1( tNumPoints0 );

                // get last values in textfile
                real tBn = tB0( tCount - 1 );
                real tHn = tH0( tCount - 1 );

                linspace( 0.0, tBn, tNumPoints0, tB1 );
                OneDMapper tMapper( tB1, 1 );
                tMapper.project( tB0, tH0, tH1 );

                // fix number
                tH1( 0 ) = 0.0;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3 : compute derivative at final point
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                Vector <real> tV( 5 );
                Vector <real> tW( 5 );
                Vector <real> tC( 3 );

                tV( 0 ) = tB1( tNumPoints0 - 5 );
                tV( 1 ) = tB1( tNumPoints0 - 4 );
                tV( 2 ) = tB1( tNumPoints0 - 3 );
                tV( 3 ) = tB1( tNumPoints0 - 2 );
                tV( 4 ) = tB1( tNumPoints0 - 1 );

                tW( 0 ) = tH1( tNumPoints0 - 5 );
                tW( 1 ) = tH1( tNumPoints0 - 4 );
                tW( 2 ) = tH1( tNumPoints0 - 3 );
                tW( 3 ) = tH1( tNumPoints0 - 2 );
                tW( 4 ) = tH1( tNumPoints0 - 1 );

                polyfit( tV, tW, 2, tC );
                real tDHDBn = dpolyval( tC, tBn );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4 : create polynomial for extrapolation
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                real tBmax = aMaxB.first;

                BELFEM_ERROR( tBmax > tBn, "max value of b mus be bigger than in textfile!" );
                real tHmax = tHn + ( tBmax - tBn ) * constant::nu0;

                create_beam_poly( tBn, tHn, tDHDBn, tBmax, tHmax, constant::nu0, tC );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5 : populate extended table
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                uint tNumPoints1 = tNumPoints0 * tBmax / tBn + 1;
                Vector <real> tB2( tNumPoints1 );
                Vector <real> tH2( tNumPoints1 );
                Vector <real> tY2( tNumPoints1 );
                real tDeltaB = tB1( 1 );
                for ( uint k = 0; k < tNumPoints0; ++k )
                {
                    tB2( k ) = tB1( k );
                    tH2( k ) = tH1( k );
                }
                for ( uint k = tNumPoints0; k < tNumPoints1; ++k )
                {
                    tB2( k ) = tB2( k - 1 ) + tDeltaB;
                    tH2( k ) = polyval( tC, tB2( k ));
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6 : compute value for lnmu
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                for ( uint k = 1; k < tNumPoints1; ++k )
                {
                    tY2( k ) = std::log( tH2( k ) / tB2( k ) );
                }

                // extrapolate first value
                tY2( 0 ) = tY2( 1 ) - ( tY2( 2 ) - tY2( 1 )) / ( tB2( 2 ) - tB2( 1 )) * tB2( 1 );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
                // Step 8 : create the spline
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -


                // help matrix for spline
                SpMatrix tHelpMatrix;
                spline::create_helpmatrix( tNumPoints1, tB2( 1 ), tHelpMatrix );

                // create the spline and send its parameters to the other procs
                return new Spline( tB2, tY2, tHelpMatrix, 0, 0,0  );
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
    }
}