//
// Created by Christian Messe on 02.12.19.
//

#include "cl_IWG_TATCAD.hpp"

#include "commtools.hpp"

#include "cl_FEM_Group.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Kernel.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

#include "cl_Gas.hpp"
#include "geometrytools.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_TATCAD::IWG_TATCAD( const uint aNumberOfDimensions ) :
                IWG_TransientHeatConduction(
                        aNumberOfDimensions,
                        IwgType::TATCAD )
        {
            this->add_fields( { "Cp", "St", "r", "hw0", "cpw0" } );
        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::set_gas( Gas * aGas )
        {
            mGas = aGas;
        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::set_wetted_sidesets_hotcold(
                const Vector< id_t > & aHotgasSidesets,
                const Vector< id_t > & aColdgasSidesets )
        {
            // combine vectors
            Vector< id_t > tAllSidesets(
                     aHotgasSidesets.length()
                    +aColdgasSidesets.length() ) ;

            index_t tCount = 0 ;
            for( id_t tID : aHotgasSidesets )
            {
                tAllSidesets( tCount++ ) = tID ;
            }
            for( id_t tID : aColdgasSidesets )
            {
                tAllSidesets( tCount++ ) = tID ;
            }

            // call the function from the parent
            IWG::set_wetted_sidesets( tAllSidesets );

            BELFEM_ERROR( mField != nullptr, "IWG does not seem to be linked agaist a field" );

            // find nodes on hotgas side
            Mesh * tMesh = mField->mesh() ;

            BELFEM_ERROR( tMesh != nullptr, "the mesh of the field seems to be a null pointer" );

            // loop over all hotgas sidesets
            for ( id_t tID : aHotgasSidesets )
            {
                // check if sideset exists
                if( tMesh->sideset_exists( tID ) )
                {
                    // flag all nodes on that sideset
                    tMesh->sideset( tID )->flag_all_nodes() ;
                }
            }

            // get ref to node container
            Cell< mesh::Node * > & tNodes = tMesh->nodes() ;

            // reset counter
            tCount = 0 ;

            // count flagged nodes
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // allocate memory
            mNodesOnHotgasSiteSets.set_size( tCount );

            // reset counter
            tCount = 0 ;
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    mNodesOnHotgasSiteSets( tCount++ ) = tNode->index() ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::set_surface_emmissivity( const real & aEpsilon )
        {
            mEpsilon = aEpsilon ;
        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::set_reference_conditions(
                real & aT,
                real & aP,
                real & aU )
        {


            if( mField->parent()->master() == mField->rank() )
            {
                mT = aT;
                mP = aP;
                mU = aU;

                // create a cell with data to send
                Vector< real > tData( 3 );
                tData( 0 ) = aT ;
                tData( 1 ) = aP ;
                tData( 2 ) = aU ;

                // send data
                send_same( mField->parent()->comm_table(), tData );
            }
            else
            {
                // data to communicate
                Vector< real > tData;

                receive( mField->parent()->master(), tData );

                mT = tData( 0 );
                mP = tData( 1 );
                mU = tData( 2 );

                // make also available for input data
                // this is not necessary, but cleaner
                aT = mT ;
                aP = mP ;
                aU = mU ;
            }

            mRho = mGas->rho( aT, aP );
            mH   = mGas->h  ( aT, aP );
            mHt  = mH + 0.5 * mU * mU;

            mTpow4 = aT * aT * aT * aT;
            mRhoU  = mRho * mU;
            mRhoU2 = mRhoU * mU;

        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::compute_heatloads()
        {
            // get mesh
            Mesh * tMesh = mField->mesh() ;


            if( mField->rank() == mField->parent()->master() )
            {
                // - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -
                // get data from mesh
                // - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - -

                // dummy vector if this is a 2d problem
                Vector< real > tEmpty ;

                // surface temperature at this timestep
                Vector< real > & tTw = tMesh->field_data( "T" );

                // surface temperature at last timestep
                Vector< real > & tTw0 = tMesh->field_data( "T0" );

                // surface enthalpy at last timestep
                Vector< real > & tHw0 = tMesh->field_data("hw0");

                // specific heat at last timestep
                Vector< real > & tCpw0 = tMesh->field_data( "cpw0" );

                // Stanton number
                Vector< real > & tSt = tMesh->field_data("St" );

                // recovery factor
                Vector< real > & tR = tMesh->field_data("r");

                Vector< real > & tDotQ = tMesh->field_data("dotQ");

                // enthalpy at wall
                real tHw ;

                // recovery enthalpy
                real tHr ;

                // scalar heatflux
                real tdotQ ;


                // loop over all nodes on the sidesets
                for ( index_t k : mNodesOnHotgasSiteSets )
                {
                    // enthalpy at wall
                    tHw = tHw0( k ) + tCpw0( k ) * ( tTw( k ) - tTw0( k ) ) * this->theta() ;

                    // compute recovery enthalpy
                    tHr = mH + tR( k ) * ( mHt - mH );

                    // compute convective heatflux
                    tdotQ = tSt( k ) * mRhoU * ( tHr - tHw );

                    // compute radiative heat flux
                    tdotQ -= constant::sigma * mEpsilon * ( std::pow( tTw( k ), 4 ) - mTpow4 );

                    // write value into node
                    tDotQ( k ) = tdotQ  ;
                }
            }

            // wait for all procs to arrive here
            comm_barrier() ;

            mField->distribute_fields( { "dotQ" } );

        }

//------------------------------------------------------------------------------

        void
        IWG_TATCAD::compute_surface_enthalpy()
        {
            // the gas class is not fully parallelized
            // we compute the whole thing on the master proc
            if( mField->rank() == mField->parent()->master() )
            {
                // get mesh
                Mesh * tMesh = mField->mesh() ;

                // get fields

                // surface temperature
                Vector< real > & tTw0 = tMesh->field_data( "T0" );

                // pressure coefficient
                Vector< real > & tCp = tMesh->field_data( "Cp" );

                // surface enthalpy
                Vector< real > & tHw0 = tMesh->field_data( "hw0" );


                // specific enthalpy at surface
                Vector< real > & tCpW0 = tMesh->field_data( "cpw0" );

                // pressure at wall
                real tPw ;

                // loop over all nodes on the sideset
                for ( index_t k : mNodesOnHotgasSiteSets )
                {
                    // compute surface pressure
                    // if this is an ideal gas, it is OK if tCp is zero
                    tPw = mP + 0.5 * mRhoU2 * tCp( k );

                    // compute enthalpy
                    tHw0( k ) = mGas->h( tTw0( k ), tPw );

                    // compute specific heat capacity
                    tCpW0( k ) = mGas->cp( tTw0( k ), tPw );

                }
            }

            // wait for master
            comm_barrier() ;

            // send relevant data to other procs
            mField->distribute_fields( { "hw0", "cpw0", "St", "r" } );
        }

//------------------------------------------------------------------------------
    }
}
