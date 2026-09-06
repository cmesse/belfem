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


#include "assert.hpp"
#include "fn_sum.hpp"
#include "fn_dot.hpp"
#include "cl_Cell.hpp"
#include "cl_Material_Abundance.hpp"
namespace belfem
{
    namespace material
    {
        Abundance::Abundance()
        {
            // taken from CRC Handbook of Chemistry and Physics

            // 1
            mDatabase[ "H" ] = AbundanceEntry { { 1.0078250321, 2.0141017780, 3.0160492675 } ,
                                               { 99.9850, 0.0115, 0.000149 }};

            // 2
            mDatabase[ "He" ] = AbundanceEntry { { 3.0160293097, 4.0026032497} ,
                                               { 0.000137, 99.999863}};

            // 3
            mDatabase[ "Li" ] = AbundanceEntry { { 6.0151223, 7.0160040} ,
                                               { 7.59, 92.41}};

            // 4
            mDatabase[ "Be" ] = AbundanceEntry { { 9.0121821 } ,
                                               { 100}};

            // 5
            mDatabase[ "B" ] = AbundanceEntry { { 10.0129370, 11.0093055} ,
                                               { 19.9, 80.1}};

            // 6
            mDatabase[ "C" ] = AbundanceEntry { { 12.0000000, 13.0033548378} ,
                                               { 98.93, 1.07}};

            // 7
            mDatabase[ "N" ] = AbundanceEntry { { 14.0030740052, 15.0001088984} ,
                                               { 99.632, 0.368}};

            // 8
            mDatabase[ "O" ] = AbundanceEntry { { 15.9949146221, 16.99913150, 17.9991604} ,
                                               { 99.757, 0.038, 0.205}};

            // 9
            mDatabase[ "F" ] = AbundanceEntry { { 18.99840320} , {100}};

            // 10
            mDatabase[ "Ne" ] = AbundanceEntry { { 19.9924401759, 20.99384674, 21.99138551} ,
                                               { 90.48, 0.27, 9.25}};

            // 11
            mDatabase[ "Na "] = AbundanceEntry { { 22.98976967 }, { 100 } };

            // 12
            mDatabase[ "Mg" ] = AbundanceEntry { { 23.98504190, 24.98583702, 25.98259304} ,
                                               { 78.99, 10.00, 11.01}};

            // 13
            mDatabase[ "Al" ] = AbundanceEntry { { 26.98153844 }, {100}};

            // 14
            mDatabase[ "Si" ] = AbundanceEntry { { 27.9769265327, 28.97649472, 29.97377022} ,
                                               { 92.2297, 4.6832, 3.0872}};

            // 15
            mDatabase[ "P" ] = AbundanceEntry { { 30.97376151 }, {100}};

            // 16
            mDatabase[ "S" ] = AbundanceEntry { { 31.97207069, 32.97145850, 33.96786683, 35.96708088} ,
                                               { 94.93, 0.76, 4.29, 0.02}};

            // 17
            mDatabase[ "Cl" ] = AbundanceEntry { { 34.96885271, 36.96590260} ,
                                               { 75.78, 24.22}};

            // 18
            mDatabase[ "Ar" ] = AbundanceEntry { { 35.96754628, 37.9627322, 39.962383123} ,
                                               { 0.3365, 0.0632, 99.6003}};

            // 19
            mDatabase[ "K" ] = AbundanceEntry { { 38.9637069, 39.96399867, 40.96182597} ,
                                               { 93.2581, 0.0117, 6.7302}};

            // 20
            mDatabase[ "Ca" ] = AbundanceEntry { { 39.9625912, 41.9586183, 42.9587668, 43.9554811, 45.9536928, 47.952534} ,
                                               { 96.941, 0.647, 0.135, 2.086, 0.004, 0.187}};

            // 21
            mDatabase[ "Sc" ] = AbundanceEntry { { 44.9559102}, {100}};

            // 22
            mDatabase[ "Ti" ] = AbundanceEntry { { 45.9526295, 46.9517638, 47.9479471, 48.9478708, 49.9447921} ,
                                               { 8.25, 7.44, 73.72, 5.41, 5.18}};

            // 23
            mDatabase[ "V" ] = AbundanceEntry { { 49.9471628, 50.9439637} ,
                                               { 0.250, 99.750}};

            // 24
            mDatabase[ "Cr" ] = AbundanceEntry { { 49.9460496, 51.9405119, 52.9406538, 53.9388849} ,
                                               { 4.345, 83.789, 9.501, 2.365}};

            // 25
            mDatabase[ "Mn" ] = AbundanceEntry { { 54.9380496} , { 100}} ;

            // 26
            mDatabase[ "Fe" ] = AbundanceEntry { { 53.9396148, 55.9349421, 56.9353987, 57.9332805} ,
                                               { 5.845, 91.754, 2.119, 0.282}};

            // 27
            mDatabase[ "Co" ] = AbundanceEntry { { 58.9331950} , { 100}} ;

            // 28
            mDatabase[ "Ni" ] = AbundanceEntry { { 57.9353479, 59.9307906, 60.9310604, 61.9283488, 63.9279696} ,
                                               { 68.0769, 26.2231, 1.1399, 3.6345, 0.9256}};

            // 29
            mDatabase[ "Cu" ] = AbundanceEntry { { 62.9296011, 64.9277937} ,
                                               { 69.17, 30.83}};

            // 30
            mDatabase[ "Zn" ] = AbundanceEntry { { 63.9291466, 65.9260368, 66.9271309, 67.9248476, 69.925325} ,
                                               { 48.63, 27.90, 4.10, 18.75, 0.62}};

            // 31
            mDatabase[ "Ga" ] = AbundanceEntry { { 68.925581, 70.9247050} ,
                                               { 60.108, 39.892}};

            // 32
            mDatabase[ "Ge" ] = AbundanceEntry { { 69.9242504, 71.9220762, 72.9234594, 73.9211782, 75.9214027} ,
                                               { 20.84, 27.54, 7.73, 36.28, 7.61}};

            // 33
            mDatabase[ "As" ] = AbundanceEntry { { 74.9215964} , {100} };

            // 34
            mDatabase[ "Se" ] = AbundanceEntry { { 73.9224766, 75.9192141, 76.9199146, 77.9173095, 79.9165218, 81.9167000} ,
                                               { 0.89, 9.37, 7.63, 23.77, 49.61, 8.73}};

            // 35
            mDatabase[ "Br" ] = AbundanceEntry { { 78.9183376, 80.916291} ,
                                               { 50.69, 49.31}};

            // 36
            mDatabase[ "Kr" ] = AbundanceEntry { { 77.920386, 79.916378, 81.9134846, 82.914136, 83.911507, 85.9106103} ,
                                               { 0.35, 2.28, 11.58, 11.49, 57.00, 17.30}};

            // 37
            mDatabase[ "Rb" ] = AbundanceEntry { { 84.9117893, 86.9091835} ,
                                               { 72.17, 27.83}};

            // 38
            mDatabase[ "Sr" ] = AbundanceEntry { { 83.913425, 85.9092624, 86.9088793, 87.9056143} ,
                                               { 0.56, 9.86, 7.00, 82.58}};

            // 39
            mDatabase[ "Y" ] = AbundanceEntry { { 88.9058479} , {100.} };

            // 40
            mDatabase[ "Zr" ] = AbundanceEntry { { 89.9047037, 90.9056450, 91.9050401, 93.9063158, 95.908276} ,
                                               { 51.45, 11.22, 17.15, 17.38, 2.80}};

            // 41
            mDatabase[ "Nb"] = AbundanceEntry { { 92.906377}, {100}};

            // 42
            mDatabase[ "Mo"] = AbundanceEntry { { 91.906810, 93.9050876, 94.9058415, 95.9046789, 96.9060210, 97.9054078, 99.907477 },
                                                { 14.84, 9.25, 15.92, 16.68, 9.55, 24.13, 9.63}};

            // 43
            mDatabase[ "Tc" ] = AbundanceEntry{ { 96.906365,97.907216, 98.9062546 },{}};

            // 44
            mDatabase[ "Ru" ] = AbundanceEntry { { 95.907598, 97.905287, 98.9059393, 99.9042197, 100.9055822, 101.9043495, 103.905430} ,
                                               { 5.54, 1.87, 12.76, 12.60, 17.06, 31.55, 18.62}};

            // 45
            mDatabase[ "Rh" ] = { { 102.905504 }, {100} };

            // 46
            mDatabase[ "Pd" ] = AbundanceEntry { { 101.905608, 103.904035, 104.905084, 105.903483, 107.903894, 109.905152} ,
                                               { 1.02, 11.14, 22.33, 27.33, 26.46, 11.72}};

            // 47
            mDatabase[ "Ag" ] = AbundanceEntry { { 106.905093, 108.904756} ,
                                               { 51.839, 48.161}};

            // 48
            mDatabase[ "Cd" ] = AbundanceEntry { { 105.906458, 107.904183, 109.903006, 110.904182, 111.9027572, 112.9044009, 113.9033581, 115.904755} ,
                                               { 1.25, 0.89, 12.49, 12.80, 24.13, 12.22, 28.73, 7.49}};

            // 49
            mDatabase[ "In" ] = AbundanceEntry { { 112.904061, 114.903878} ,
                                               { 4.29, 95.71}};

            // 50
            mDatabase[ "Sn" ] = AbundanceEntry { { 111.904821, 113.902782, 114.903346, 115.901744, 116.902954, 117.901606, 118.903309, 119.9021966, 121.9034401, 123.9052746} ,
                                               { 0.97, 0.66, 0.34, 14.54, 7.68, 24.22, 8.59, 32.58, 4.63, 5.79}};

            // 51
            mDatabase[ "Sb" ] = AbundanceEntry { { 120.9038180, 122.9042157} ,
                                               { 57.21, 42.79}};

            // 52
            mDatabase[ "Te" ] = AbundanceEntry { { 119.904020, 121.9030471, 122.9042730, 123.9028195, 124.9044247, 125.9033055, 127.9044614, 129.9062228} ,
                                               { 0.09, 2.55, 0.89, 4.74, 7.07, 18.84, 31.74, 34.08}};

            // 53
            mDatabase[ "I" ] = AbundanceEntry { { 126.904473}, {100}};

            // 54
            mDatabase[ "Xe" ] = AbundanceEntry { { 123.9058958, 125.904269, 127.9035304, 128.9047795, 129.9035079, 130.9050819, 131.9041545, 133.9053945, 135.907220} ,
                                               { 0.09, 0.09, 1.92, 26.44, 4.08, 21.18, 26.89, 10.44, 8.87}};

            // 55
            mDatabase[ "Cs" ] = AbundanceEntry { { 132.905447}, {100}};

            // 56
            mDatabase[ "Ba" ] = AbundanceEntry { { 129.906310, 131.905056, 133.904503, 134.905683, 135.904570, 136.905821, 137.905241} ,
                                               { 0.106, 0.101, 2.417, 6.592, 7.854, 11.232, 71.698}};

            // 57
            mDatabase[ "La" ] = AbundanceEntry { { 137.907107, 138.906348} ,
                                               { 0.090, 99.910}};

            // 58
            mDatabase[ "Ce" ] = AbundanceEntry { { 135.907140, 137.905986, 139.905434, 141.909240} ,
                                               { 0.185, 0.251, 88.450, 11.114}};

            // 59
            mDatabase[ "Pr" ] = AbundanceEntry{ {140.907648}, {100}};

            // 60
            mDatabase[ "Nd" ] = AbundanceEntry { { 141.907719, 142.909810, 143.910083, 144.912569, 145.913112, 147.916889, 149.920887 } ,
                                            { 27.2, 12.2, 23.8, 8.3, 17.2, 5.7, 5.6}};

            // 61
            mDatabase[ "Pm"] = AbundanceEntry { { 144.912750, 146.914893 }, {}};

            // 62
            mDatabase[ "Sm" ] = AbundanceEntry { {  143.911995, 147.914818, 148.917180, 149.917271, 151.919728, 153.922205} ,
                                               { 3.07, 14.99, 11.24, 13.82, 7.38, 26.75, 22.75}};

            // 63
            mDatabase[ "Eu" ] = AbundanceEntry { { 150.919846, 152.921226} ,
                                               { 47.81, 52.19}};

            // 64
            mDatabase[ "Gd" ] = AbundanceEntry { { 151.919788, 153.920862, 154.922619, 155.922120, 156.923957, 157.924101, 159.927051} ,
                                               { 0.20, 2.18, 14.80, 20.47, 15.65, 24.84, 21.86}};

            // 65
            mDatabase[ "Tb" ] = AbundanceEntry {{ 158.925343 }, {100}};

            // 66
            mDatabase[ "Dy" ] = AbundanceEntry { { 157.924405, 159.925194, 160.926930, 161.926795, 162.928728, 163.929171} ,
                                               { 0.06, 0.10, 2.34, 18.91, 25.51, 24.90, 28.18}};

            // 67
            mDatabase[ "Ho" ] = AbundanceEntry {{ 164.930319 }, {100}};

            // 68
            mDatabase[ "Er" ] = AbundanceEntry { { 161.928775, 163.929197, 165.930290, 166.932045, 167.932368, 169.935460} ,
                                               { 0.14, 1.61, 33.61, 22.93, 26.78, 14.93}};
            // 69
            mDatabase[ "Tm"] = AbundanceEntry {{ 168.934211 }, {100}};
            // 70
            mDatabase[ "Yb" ] = AbundanceEntry { { 167.933894, 169.934759, 170.936322, 171.9363777, 172.9382068, 173.9388581, 175.942568} ,
                                               { 0.13, 3.04, 14.28, 21.83, 16.13, 31.83, 12.76}};

            // 71
            mDatabase[ "Lu" ] = AbundanceEntry { { 174.9407679, 175.9426824} ,
                                               { 97.41, 2.59}};

            // 72
            mDatabase[ "Hf" ] = AbundanceEntry { { 173.940040, 175.9414018, 176.9432200, 177.9436977, 178.9458151, 179.9465488} ,
                                               { 0.16, 5.26, 18.60, 27.28, 13.62, 35.08}};

            // 73
            mDatabase[ "Ta" ] = AbundanceEntry { { 179.947466, 180.947996} ,
                                               { 0.012, 99.988}};

            // 74
            mDatabase[ "W" ] = AbundanceEntry { { 179.946706, 181.948206, 182.9502245, 183.9509326, 185.954362} ,
                                               { 0.12, 26.50, 14.31, 30.64, 28.43}};

            // 75
            mDatabase[ "Re" ] = AbundanceEntry { { 184.9529557, 186.9557508} ,
                                               { 37.40, 62.60}};

            // 76
            mDatabase[ "Os" ] = AbundanceEntry { { 183.952491, 185.953838, 186.9557479, 187.9558360, 188.9581449, 189.958445, 191.961479} ,
                                               { 0.02, 1.59, 1.96, 13.24, 16.15, 26.26, 40.78}};

            // 77
            mDatabase[ "Ir" ] = AbundanceEntry { { 190.960591, 192.962924} ,
                                               { 37.3, 62.7}};

            // 78
            mDatabase[ "Pt" ] = AbundanceEntry { { 189.959930, 191.961035, 193.962664, 194.964774, 195.964935, 197.967876} ,
                                               { 0.014, 0.782, 32.967, 33.832, 25.242, 7.163}};

            // 79
            mDatabase[ "Au"] = AbundanceEntry { { 196.966552} ,
                                               { 100 }};

            // 80
            mDatabase[ "Hg" ] = AbundanceEntry { { 195.965815, 197.966752, 198.968262, 199.968309, 200.970285, 201.970626, 203.973476} ,
                                               { 0.15, 9.97, 16.87, 23.10, 13.18, 29.86, 6.87}};

            // 81
            mDatabase[ "Tl" ] = AbundanceEntry { { 202.972329, 204.974412} ,
                                               { 29.524, 70.476}};

            // 82
            mDatabase[ "Pb" ] = AbundanceEntry { { 203.973029, 205.974449, 206.975881, 207.976636} ,
                                               { 1.4, 24.1, 22.1, 52.4}};

            // 83
            mDatabase[ "Bi" ] = AbundanceEntry { { 208.980383} ,
                                               { 100}};
        }

        void
        Abundance::compute_molar_mass_and_impurity_from_volumes(
                        const Cell< string > & aElements,
                        const Vector< real > & aVolumeFractions,
                          real & aMolarMass,
                          real & aImpurity ) const
        {
            const Vector< real > & X = aVolumeFractions ;
            index_t n = aVolumeFractions.length() ;

            Vector< real > M( n, 0.0 );

            for ( index_t i = 0 ; i < n ; ++i )
            {
                const AbundanceEntry & tEntry = mDatabase( aElements( i ) );

                M( i ) = 0.01 * dot( tEntry.first, tEntry.second ) ;
            }

            aImpurity = 0.0 ;


            aMolarMass= dot( M, X ) ;

            for ( index_t i = 0 ; i < n ; ++i )
            {
                const AbundanceEntry & tEntry = mDatabase( aElements( i ) );

                const Vector< real > & Mj = tEntry.first ;
                const Vector< real > & Aj = tEntry.second ;
                uint m = Mj.length() ;
                for ( uint j=0; j<m; ++j )
                {
                    aImpurity += X( i ) * Aj( j ) * std::pow( ( Mj( j ) - M( i ) ) / aMolarMass , 2 );
                }
            }
            aImpurity *= 0.01 ;
            aMolarMass *= 0.001 ;
        }

        void
        Abundance::compute_molar_mass_and_impurity_from_masses(
                       const Cell< string > & aElements,
                       const Vector< real > & aMassFractions,
                         real & aMolarMass,
                         real & aImpurity ) const
        {
            // find balance
            Vector< real > Y( aMassFractions );
            real s = sum( Y );

            index_t n = aMassFractions.length() ;

            // check if data are given in percent
            if ( std::abs( 100. - s ) < 1e-6 )
            {
                Y *= 0.01 ;
            }

            // check if balance is correct
            else if ( std::abs(s - 1.0) > 1e-6 )
            {
                // check if one and only one entry is zero
                index_t b = 0 ;
                index_t tCount = 0 ;
                for ( index_t i = 0 ; i < n ; ++i )
                {
                    if ( Y( i ) == 0.0 )
                    {
                        b = i ;
                        ++tCount ;
                    }
                }
                BELFEM_ERROR( tCount == 1, "Mass Fractions are not balanced" );

                // create balance
                Y( b ) = 1.0 - s;
            }

            // compute the molar masses for the individual components
            Vector< real > M( n, 0.0 );

            for ( index_t i = 0 ; i < n ; ++i )
            {
                const AbundanceEntry & tEntry = mDatabase( aElements( i ) );

                M( i ) = 0.01 * dot( tEntry.first, tEntry.second ) ;
            }

            // compute the molar fractions
            Vector< real > X( n, 0.0 );

            for ( index_t i = 0 ; i < n ; ++i )
            {
                X( i ) = Y( i ) / M( i ) ;
            }
            X /= sum( X );

            aMolarMass = dot( M, X ) ;

            for ( index_t i = 0 ; i < n ; ++i )
            {
                const AbundanceEntry & tEntry = mDatabase( aElements( i ) );

                const Vector< real > & Mj = tEntry.first ;
                const Vector< real > & Aj = tEntry.second ;
                uint m = Mj.length() ;
                for ( uint j=0; j<m; ++j )
                {
                    aImpurity += X( i ) * Aj( j ) * std::pow( ( Mj( j ) - M( i ) ) / aMolarMass , 2 );
                }
            }
            aImpurity *= 0.01 ;
            aMolarMass *= 0.001 ;
        }

        real
        Abundance::compute_molar_mass(
            const string & aElement   ) const
        {
            const AbundanceEntry & tEntry = mDatabase( aElement );
            return 0.001 * 0.01 * dot( tEntry.first, tEntry.second ) ;
        }

    }
}
