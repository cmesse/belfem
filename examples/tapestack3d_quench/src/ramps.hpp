/*
 * The time programme of the tapestack3d quench run, in ONE place: the
 * transport current and the Ic defect are both driven from here, so the two
 * cannot drift apart between the deck and the plugin.
 *
 * Story of the run: a linear current ramp, as in an Ic test. Three seconds
 * in, the topmost tape loses 90 % of its jc over a 1.5 mm disk; the solder
 * diverts that tape's share into the seven healthy tapes, which have ample
 * margin at that current. The ramp does not stop, so the stack runs out of
 * margin a few seconds later and quenches, starting at the defect where the
 * diverted current crosses the resistive solder.
 *
 * The defect switch is the logistic the deck's `type : sigmoid` source uses
 * ( cl_SourceFunction.hpp, function_sigmoid ):
 *
 *     s( t ) = 1 / ( 1 + exp( beta * ( d - t ) ) ),
 *     d      = offset + period / 2,
 *     beta   = -2 ln( fuzzyness / ( 1 - fuzzyness ) ) / period,
 *
 * so s( offset ) = fuzzyness and s( offset + period ) = 1 - fuzzyness.
 * Included by defect.cpp and current.cpp; both compile into userdefect.so.
 */

#ifndef TAPESTACK3D_RAMPS_HPP
#define TAPESTACK3D_RAMPS_HPP

#include <algorithm>
#include <cmath>

namespace ramps
{
    // ---- transport current ( boundary condition MyCurrent ) ----------------
    // Linear ramp of the TOTAL stack current ( eight tapes in parallel ).
    // 200 A/s reaches 2000 A at the 10 s end of the run, above the stack's
    // self-field Ic at 77 K, so the quench happens DURING the ramp rather
    // than depending on a plateau value guessed right. sp-ap at 77.5 K,
    // self-field: 380 A/cm, i.e. 152 A per 4 mm tape with the table's 1 um
    // layer convention, 243 A with the deck's 1.6 um ybco layer.
    constexpr double gCurrentRate      = 200.0 ;   // A/s
    constexpr double gCurrentMax       = 2500.0 ;  // A, safety cap, not reached in 10 s

    // ---- Ic defect switch-on ( MyDefect ) -----------------------------------
    // Fully on at 3.5 s, when the stack carries 700 A ( ~88 A per tape ).
    // The top tape sheds ~80 A into its neighbours, which then sit near
    // 45 % of Ic: the defect alone must NOT quench the stack, the ramp does.
    constexpr double gDefectOffset     =  3.0 ;    // s
    constexpr double gDefectPeriod     =  0.5 ;    // s
    constexpr double gDefectFuzzyness  =  1.0e-3 ;

    inline double
    sigmoid( const double t, const double aOffset, const double aPeriod, const double aFuzzyness )
    {
        const double d    = aOffset + 0.5 * aPeriod ;
        const double beta = -2.0 * std::log( aFuzzyness / ( 1.0 - aFuzzyness ) ) / aPeriod ;
        return 1.0 / ( 1.0 + std::exp( beta * ( d - t ) ) );
    }

    inline double
    current( const double t )
    {
        return std::min( gCurrentRate * std::max( t, 0.0 ), gCurrentMax );
    }

    inline double
    defect_switch( const double t )
    {
        return sigmoid( t, gDefectOffset, gDefectPeriod, gDefectFuzzyness );
    }
}

#endif // TAPESTACK3D_RAMPS_HPP
