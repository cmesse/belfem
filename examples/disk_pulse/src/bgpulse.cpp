/*
 * Background-field pulse for the disk magnetization experiment.
 *
 * Ramp the applied flux density from zero to B0, hold it, then cut it to zero
 * in one step. The screening currents induced in the disk have nowhere to go
 * and keep circulating, so the disk stays magnetized after the applied field
 * is gone -- that persistence is what the run is meant to show.
 *
 * The function returns a FLUX DENSITY in tesla, which is what the deck's
 * `units : T ;` declares. BELFEM converts it to the field strength the weak
 * form imposes, H = B/mu0, and it is that H which appears in the solver.
 *
 * Build:
 *   mkdir build && cd build
 *   cmake -DBELFEM_DIR=/path/to/belfem .. && make
 */

#include <belfem_user_api>

#include <cmath>

using namespace belfem;

namespace
{
    const real gB0    = 1.0 ;      // plateau flux density          [ T ]
    const real gTUp   = 0.100 ;    // centre of the rise            [ s ]
    const real gTDown = 0.200 ;    // centre of the fall            [ s ]
    const real gTau   = 0.010 ;    // 10-90% transition ~ 4.4*gTau  [ s ]

    /**
     * logistic step, 0 -> 1, centred on t0
     */
    inline real
    sigmoid( const real t, const real t0, const real tau )
    {
        return 1.0 / ( 1.0 + std::exp( - ( t - t0 ) / tau ) );
    }
}

/**
 * @brief applied background flux density
 * @param t time (s)
 * @return  flux density (T)
 *
 * Two logistic steps: up at gTUp, down at gTDown. The difference gives a
 * plateau between them and returns to zero afterwards.
 *
 * WHY NOT A STEP. The first version of this file cut the field to zero in one
 * instant, and the run answered with J/Jc > 20. That is not stiffness, it is
 * the power law being handed an unphysical excitation: with n = 25, J/Jc = 20
 * means E/Ec = 20^25, about 3e32. The critical-state model caps |J| at Jc
 * exactly because the textbook field sweep is quasi-static -- Iwasa 2009
 * Sec. 5.2 sweeps 0 -> Hm -> 0 as a continuous sweep, and Russenschuck 2010
 * Sec. 16.1 says a time-transient field drives the current density only
 * "slightly above Jc", relaxing back to Jc once the sweep stops.
 *
 * A finite transition is therefore not a numerical convenience but the
 * physical statement. A logistic is steepest at its centre, where the slope is
 * B0/( 4*gTau ) = 25 T/s, and its 10-90% transition takes 4.4*gTau = 44 ms.
 * By E ~ ( dB/dt )*r/2 and J/Jc = ( E/Ec )^( 1/n ) that lands near J/Jc ~ 1.4 --
 * "slightly above", and still a fast cut against the 100 ms plateau.
 *
 * The two steps sit 10*gTau apart, so they overlap by exp( -5 ): the plateau
 * reaches 0.9866*B0 rather than B0 exactly. Widen the separation or shrink
 * gTau if a flat top matters more than a sharp cut.
 *
 * The physics under test is unchanged: after the field returns to zero the
 * screening currents have nowhere to go, and the disk keeps a remanent
 * magnetization ( Iwasa Eqs. 5.6-5.7, -M = Hp/2 for a Bean slab ) which no
 * external field can remove -- only warming past Tc does that.
 */
real background_pulse( const real t )
{
    return gB0 * (   sigmoid( t, gTUp,   gTau )
                   - sigmoid( t, gTDown, gTau ) );
}

extern "C" void BgPulse_init( SourceFunction * source )
{
    source->set_user_defined( &background_pulse );
}
