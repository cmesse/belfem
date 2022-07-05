//
// Created by Christian Messe on 30.04.20.
//

#ifndef BELFEM_CL_TABLEGAS_HPP
#define BELFEM_CL_TABLEGAS_HPP

#include "cl_Gas.hpp"
#include "cl_BS_LookupTable.hpp"
#include "fn_GT_data_path.hpp"


namespace belfem
{
    class TableGas : public Gas
    {

        bspline::LookupTable * mTable;

        // index for mass table
        index_t mIndexM;

        // index for enthalpy table
        index_t mIndexH;

        // index for entropy table
        index_t mIndexS;

        // index for viscisity table
        index_t mIndexMu;

        // index for frozen conductivity table
        index_t mIndexLambda;

        // index for dissociation enthalpy table
        index_t mIndexHd ;

        // constant for scaling derivative to p
        const real mdpscale = 1000.0 / std::log( 10.0 );

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        TableGas( const string aPath = "" );

//------------------------------------------------------------------------------

        ~TableGas();

//------------------------------------------------------------------------------

        void
        remix( const Vector<real> & aMolarFractions,
               bool aRemixHeat=true,
               bool aRemixTransport=true );

//------------------------------------------------------------------------------

        void
        remix_mass( const Vector<real> & aMassFractions,
                    bool aRemixHeat=true,
                    bool aRemixTransport=true );

//------------------------------------------------------------------------------

        void
        reset_mixture();

//------------------------------------------------------------------------------

        // expose the lookup table
        inline bspline::LookupTable &
        lookup_table()
        {
            return *mTable;
        };

//------------------------------------------------------------------------------
// MIXTURE
//------------------------------------------------------------------------------

        // Molar weight
        const real &
        M( const real & aT, const real & aP );

        // specific gas constant
        const real &
        R( const real & aT, const real & aP );

//------------------------------------------------------------------------------
// Caloric Properties
//------------------------------------------------------------------------------

        real
        cp( const real & aT, const real & aP );

        real
        cv( const real & aT, const real & aP );

        real
        gamma( const real & aT, const real & aP );

        real
        c( const real & aT, const real & aP );

        real
        h( const real & aT, const real & aP );

        real
        s( const real & aT, const real & aP );

        real
        dsdT( const real & aT, const real & aP );

        real
        dsdp( const real & aT, const real & aP );

        real
        dcpdT( const real & aT, const real & aP );

        // enthalpy derivative to pressure ( needed for total temperature )
        real
        dhdp( const real & aT, const real & aP );

        real
        hd( const real & aT, const real & aP );

//------------------------------------------------------------------------------
// Transport Properties
//------------------------------------------------------------------------------

        real
        mu( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        lambda( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        Pr( const real & aT, const real & aP );



//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        // special parameter for table
        real
        pi( const real & aT, const real & aP );

        void
        create_table( const string aPath );

//------------------------------------------------------------------------------

        void
        init_parameters();

//------------------------------------------------------------------------------

        void
        create_eos();

//------------------------------------------------------------------------------




    };
}
#endif //BELFEM_CL_TABLEGAS_HPP
