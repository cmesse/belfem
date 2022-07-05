//
// Created by Christian Messe on 15.09.19.
//

#include "cl_GM_EoS_AlphaFunctionFactory.hpp"
#include "cl_GT_GasData.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_empty()
        {
            return new AlphaFunction_Empty();
        }
//----------------------------------------------------------------------------
        AlphaFunction *
        AlphaFunctionFactory::create_srk( const gastables::GasData * aData )
        {
            real tOmega = aData->acentric();

            // calculate M
            // 10.1016/0009-2509(72)80096-4
            real tM = 0.480 + ( 1.574 - 0.175*tOmega ) * tOmega;

            // help factor
            real tHelp = 1.0/std::sqrt( aData->T_crit() );

            return new AlphaFunction_Classic( aData->T_crit(), tM, tHelp );
        }

//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_pr78( const gastables::GasData * aData )
        {
            real tTcrit = aData->T_crit();
            real tOmega = aData->acentric();

            // calculate M
            real tM;

            if( tOmega < 0.48225818664150998277716914005452 )
            {
                tM = 0.37464 + ( 1.54226 - 0.26992 * tOmega ) * tOmega;
            }
            else
            {
                tM = 0.374642 + ( 1.487503 + ( 0.016666*tOmega - 0.164423 ) *tOmega ) * tOmega;
            }


            // help factor
            real tHelp = 1.0/std::sqrt( tTcrit );

            return new AlphaFunction_Classic( tTcrit, tM, tHelp );
        }

//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_ccr_pr( const gastables::GasData * aData )
        {
            // get number of gas
            const string & tCAS = aData->cas();

            real tC1 = BELFEM_QUIET_NAN;
            real tC2 = BELFEM_QUIET_NAN;
            real tC3 = BELFEM_QUIET_NAN;

            // 10.1023/B:IJOT.0000022331.46865.2f
            // Table V
            if( tCAS == "1333-74-0" )
            {
                // hydrogen
                tC1 =  0.09406;
                tC2 = -0.22429;
                tC3 = -0.02458;
            }
            else if ( tCAS == "74-82-8" )
            {
                // methane
                tC1 =  0.41667;
                tC2 = -0.05156;
                tC3 =  0.38954;
            }
            else if ( tCAS == "7782-44-7" )
            {
                // oxygen
                tC1 =  0.41325;
                tC2 =  0.10376;
                tC3 =  0.10971;
            }
            else if ( tCAS == "7727-37-9" )
            {
                // nitrogen
                tC1 =  0.44950;
                tC2 = -0.03278;
                tC3 =  0.49308;
            }
            else if ( tCAS == "74-85-1" )
            {
                // ethylene
                tC1 =  0.51014;
                tC2 =  0.06247;
                tC3 =  0.32052;
            }
            else if ( tCAS == "04.06.83" )
            {
                // hydrogen sulfide
                tC1 =  0.50694;
                tC2 =  0.14188;
                tC3 =  0.31438;
            }
            else if ( tCAS == "74-84-0" )
            {
                // ethane
                tC1 =  0.52539;
                tC2 =  0.11674;
                tC3 =  0.13968;
            }
            else if ( tCAS == "74-98-6" )
            {
                // propane
                tC1 =  0.59311;
                tC2 =  0.17042;
                tC3 =  0.10182;
            }
            else if ( tCAS == "75-28-50" )
            {
                // isobutane
                tC1 =  0.64121;
                tC2 =  0.07005;
                tC3 =  0.42647;
            }
            else if ( tCAS == "106-97-8" )
            {
                // n-Buthane
                tC1 =  0.67084;
                tC2 =  0.09474;
                tC3 =  0.23091;
            }
            else if ( tCAS == "110-82-7" )
            {
                // cyclohexane
                tC1 =  0.68259;
                tC2 =  0.04522;
                tC3 =  0.53089;
            }
            else if ( tCAS == "71-43-2" )
            {
                // benzene
                tC1 =  0.69709;
                tC2 = -0.07749;
                tC3 =  0.86396;
            }
            else if ( tCAS == "124-38-9" )
            {
                // carbon dioxide
                tC1 =  0.68583;
                tC2 =  0.17408;
                tC3 =  0.18239;
            }
            else if ( tCAS == "78-78-4" )
            {
                // isopentane
                tC1 =  0.71103;
                tC2 =  0.06958;
                tC3 =  0.29784;
            }
            else if ( tCAS == "109-66-0" )
            {
                // pentane
                tC1 =  0.74373;
                tC2 =  0.05868;
                tC3 =  0.35254;
            }
            else if ( tCAS == "7664-41-7" )
            {
                // ammonia
                tC1 =  0.74852;
                tC2 =  0.07849;
                tC3 =  0.10073;
            }
            else if ( tCAS == "108-88-3" )
            {
                // toulene
                tC1 =  0.75554;
                tC2 =  0.11290;
                tC3 =  0.22419;
            }
            else if ( tCAS == "110-54-3" )
            {
                // hexane
                tC1 =  0.83968;
                tC2 = -0.19125;
                tC3 =  0.93864;
            }
            else if ( tCAS == "67-64-1" )
            {
                // acetone
                tC1 =  0.82577;
                tC2 =  0.04252;
                tC3 =  0.15901;
            }
            else if ( tCAS == "7732-18-5" )
            {
                // water
                tC1 =  0.91402;
                tC2 = -0.23571;
                tC3 =  0.54115;
            }
            else if ( tCAS == "142-82-5" )
            {
                // heptane
                tC1 =  0.87206;
                tC2 =  0.08945;
                tC3 =  0.28459;
            }
            else if ( tCAS == "111-65-9" )
            {
                // octane
                tC1 =  0.94934;
                tC2 = -0.00379;
                tC3 =  0.43788;
            }
            else
            {
                real tOmega = aData->acentric();

                // Equation ( 26 )
                tC1 = ( 1.3569 * tOmega + 0.9957 ) * tOmega + 0.4077;

                // Equation ( 27 )
                tC2 = ( 3.5590 -11.2986 * tOmega ) * tOmega - 0.1146;

                // Equation ( 28 )
                tC3 = ( 11.7802 * tOmega - 3.8901 ) * tOmega + 0.5033;
            }

            return new AlphaFunction_CCR( aData->T_crit(), tC1, tC2, tC3 );
        }

//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_ccr_mc_srk( const gastables::GasData * aData )
        {
            // get number of gas
            const string & tCAS = aData->cas();

            real tC1 = BELFEM_QUIET_NAN;
            real tC2 = BELFEM_QUIET_NAN;
            real tC3 = BELFEM_QUIET_NAN;

            // 10.1023/B:IJOT.0000022331.46865.2f
            // Table II
            if( tCAS == "1333-74-0" )
            {
                // hydrogen
                tC1 =  0.161;
                tC2 = -0.225;
                tC3 = -0.232;
            }
            else if ( tCAS == "74-82-8" )
            {
                // methane
                tC1 =  0.549;
                tC2 = -0.409;
                tC3 =  0.603;
            }
            else if ( tCAS == "7782-44-7" )
            {
                // oxygen
                tC1 =  0.545;
                tC2 = -0.235;
                tC3 =  0.292;
            }
            else if ( tCAS == "7727-37-9" )
            {
                // nitrogen
                tC1 =  0.584;
                tC2 = -0.396;
                tC3 =  0.736;
            }
            else if ( tCAS == "74-85-1" )
            {
                // ethylene
                tC1 =  0.652;
                tC2 = -0.315;
                tC3 =  0.563;
            }
            else if ( tCAS == "04.06.83" )
            {
                // hydrogen sulfide
                tC1 =  0.641;
                tC2 = -0.183;
                tC3 =  0.513;
            }
            else if ( tCAS == "74-84-0" )
            {
                // ethane
                tC1 =  0.711;
                tC2 = -0.573;
                tC3 =  0.894;
            }
            else if ( tCAS == "74-98-6" )
            {
                // propane
                tC1 =  0.775;
                tC2 = -0.476;
                tC3 =  0.815;
            }
            else if ( tCAS == "75-28-50" )
            {
                // isobutane
                tC1 =  0.807;
                tC2 = -0.432;
                tC3 =  0.91;
            }
            else if ( tCAS == "106-97-8" )
            {
                // n-Buthane
                tC1 =  0.823;
                tC2 = -0.267;
                tC3 =  0.402;
            }
            else if ( tCAS == "110-82-7" )
            {
                // cyclohexane
                tC1 =  0.86;
                tC2 = -0.566;
                tC3 =  1.375;
            }
            else if ( tCAS == "71-43-2" )
            {
                // benzene
                tC1 =  0.84;
                tC2 = -0.389;
                tC3 =  0.917;
            }
            else if ( tCAS == "124-38-9" )
            {
                // carbon dioxide
                tC1 =  0.867;
                tC2 = -0.674;
                tC3 =  2.471;
            }
            else if ( tCAS == "78-78-4" )
            {
                // isopentane
                tC1 =  0.876;
                tC2 = -0.386;
                tC3 =  0.66;
            }
            else if ( tCAS == "109-66-0" )
            {
                // pentane
                tC1 =  0.901;
                tC2 = -0.305;
                tC3 =  0.542;
            }
            else if ( tCAS == "7664-41-7" )
            {
                // ammonia
                tC1 =  0.916;
                tC2 = -0.369;
                tC3 =  0.417;
            }
            else if ( tCAS == "108-88-3" )
            {
                // toulene
                tC1 =  0.923;
                tC2 = -0.301;
                tC3 =  0.494;
            }
            else if ( tCAS == "110-54-3" )
            {
                // hexane
                tC1 =  1.005;
                tC2 = -0.591;
                tC3 =  1.203;
            }
            else if ( tCAS == "67-64-1" )
            {
                // acetone
                tC1 =  0.993;
                tC2 = -0.322;
                tC3 =  0.265;
            }
            else if ( tCAS == "7732-18-5" )
            {
                // water
                tC1 =  1.095;
                tC2 = -0.678;
                tC3 =  0.7;
            }
            else if ( tCAS == "142-82-5" )
            {
                // heptane
                tC1 =  1.036;
                tC2 = -0.258;
                tC3 =  0.488;
            }
            else if ( tCAS == "111-65-9" )
            {
                // octane
                tC1 =  1.15;
                tC2 = -0.587;
                tC3 = 1.096;
            }
            else
            {
                real tOmega = aData->acentric();

                // Equation ( 14 )
                tC1 = ( 0.16054 - 0.1094 * tOmega  ) * tOmega + 0.5178;

                // Equation ( 15 )
                tC2 = 0.3279 - 0.4291 * tOmega;

                // Equation ( 16 )
                tC3 = 1.3506 * tOmega + 0.4866;
            }

            return new AlphaFunction_MC( aData->T_crit(), tC1, tC2, tC3 );
        }

//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_pm_srk(  const gastables::GasData * aData )
        {
            BELFEM_ASSERT( aData->has_cubic(),
                    "The gas %s has no cubic data",
                    aData->label().c_str() );

            return new AlphaFunction_PM( aData->T_crit(), aData->srk() );
        }

//----------------------------------------------------------------------------

        AlphaFunction *
        AlphaFunctionFactory::create_pm_pr(  const gastables::GasData * aData )
        {
            BELFEM_ASSERT( aData->has_cubic(),
                          "The gas %s has no cubic data",
                          aData->label().c_str() );

            return new AlphaFunction_PM( aData->T_crit(), aData->pr() );
        }

//----------------------------------------------------------------------------
    }
}
