// Defines the entry point for the console application.
//Pleak is included in Pcool

//for bolometers (absorber + 2SINs) w/ Andreev Current
#ifndef CEB_2EQ_PARALLEL_LITE_H
#define CEB_2EQ_PARALLEL_LITE_H

#include <cmath>
#include <tuple>

#include "constants.h"

/*
 * Never used
 */
template<class ScalarLike>
double current(const ScalarLike& v, const double tau) {
    /*  Compute the current using approximation [delta(0) / eRn]
     *
     *  `v` is voltage,
     *  `tau` is for temperature
     */

    const double a0 = 1.0 + 0.375 * tau - 0.1171875 * std::pow(tau, 2);
    const double a1 = tau * std::pow(a0, 2);
    ScalarLike a2 = 1.0 + std::exp((std::abs(v) - 1.0) / tau - (1.15 + tau));
    ScalarLike a3 = (std::pow(v, 2) - 1.0) / (1.0 + std::exp(-(std::pow(v, 2) - 1.0) / tau));

    ScalarLike s = (v >= 0.0) * 2.0 - 1.0;

    // [delta(0) / eRn]
    return s * std::sqrt(2.0 * M_PI * a1 / a2 + a3)
           * (1.0 / (2.0 * std::exp((1.0 - v) / tau) + 1.0)
              - 1.0 / (2.0 * std::exp((1.0 + v) / tau) + 1.0));
}

template<class ScalarLike>
double currentIntegral(const double DT, const ScalarLike& v, const double tau, const double tauE) {
    /*
     * Compute the current using exact integral, rectangle method
     *
     * @vshampor: This is probably a function to calculate the integral value in dimensionless units
     *  and parameters should be renamed accordingly to move the physical meaning out of this function's params to
     *  a position higher in code
     *
     * `DT` is delta (energy gap),
     * `v` is voltage,
     * `tau` is for temperature,
     * `tauE` is for electron temperature
    */
    const double dx = INTEGRATION_SCALE * DT;
    const double inf = ESSENTIALLY_INFINITY_SCALE * DT;
    const double DT2 = std::pow(DT, 2);
    const ScalarLike v1 = std::abs(v);

    ScalarLike i = 0.0;

    double x = DT + dx; // energy
    while (x <= inf) {
        constexpr double accuracy = 1e-15;
        const double a2 = std::sqrt(std::pow(x, 2) - DT2);
        ScalarLike a0;
        double a1;
        ScalarLike i_step;

        a0 = 1.0 / (std::exp((x - v1) / tauE) + 1.0);
        a1 = 1.0 / (std::exp(x / tau) + 1.0);
        i_step = x / a2 * (a0 - a1);
        // accuracy check at every iteration
        if (std::abs(i_step) < accuracy)
            break;
        i += i_step;

        a0 = 1.0 / (std::exp((-x - v1) / tauE) + 1.0);
        a1 = 1.0 / (std::exp(-x / tau) + 1.0);
        i_step = x / a2 * (a0 - a1);
        // accuracy check at every iteration
        if (std::abs(i_step) < accuracy)
            break;
        i += i_step;

        x += dx;
    }

    ScalarLike s = (v >= 0.0) * 2.0 - 1.0;

    return i * dx * s;
}

template<class ScalarLike>
double PowerCoolPiece(const ScalarLike& v, const double tau) {
    /*
     *  A part of `PowerCool`
     */
    const ScalarLike a0 = (1.0 - v) / tau;
    const ScalarLike tmp1 = 2.0 * std::exp(a0);
    const ScalarLike tmp2 = std::exp(-2.5 * (a0 + 2.0));
    const double PItau = M_PI * tau;
    const ScalarLike a1 = std::sqrt(2.0 * PItau) * ((1.0 - v) / (tmp1 + 1.28) + 0.5 * tau / (tmp1 + 0.64))
                          / (tmp2 + 1.0);

    if (v - 1.0 > tau) {
        const ScalarLike a2 = std::sqrt(std::pow(v, 2) - 1.0);
        return a1 + 0.5 * (-v * a2 + std::log(v + a2) + std::pow(PItau, 2) / 3.0 * v / a2) / (1.0 / tmp2 + 1.0);
    }
    return a1;
}

template<class ScalarLike>
double PowerCool(const ScalarLike& v, const double tau, const double tauE) {
    /*
     *  Compute the cooling power of a SIN junction
     *
     *  `v` is voltage,
     *  `tau` is for temperature,
     *  `tauE` is for electron temperature
     */

    const ScalarLike a1 = PowerCoolPiece(v, tauE);
    const ScalarLike b1 = PowerCoolPiece(-v, tauE);
    const ScalarLike c1 = PowerCoolPiece(0.0, tau);

    return a1 + b1 - c1;
}

template<class ScalarLike>
ScalarLike AndCurrent(const double DT, const ScalarLike& v, const double tauE, const double Wt, const double tm) {
    /*
     *  Compute the Andreev current using exact integral, rectangle method
     *
     *  `DT` is delta (energy gap),
     *  `v` is voltage,
     *  `tauE` is for electron temperature,
     *  `Wt` is omega with ~ from Vasenko 2010, for the transparency of the barrier,
     *  `tm` is for the depairing energy
     */

    const double dx = INTEGRATION_SCALE * DT;
    const double DT2 = std::pow(DT, 2);
    const double twoWt = 2.0 * Wt;
    const double twoTauE = 2.0 * tauE;

    ScalarLike i = 0.0;

    double x = 0.0; // energy
    while (x <= DT) {
        const double DT2x2 = DT2 - std::pow(x, 2);
        const double sqrtDT2x2 = std::sqrt(DT2x2);
        const ScalarLike da0a1 = std::tanh((x + v) / twoTauE) - std::tanh((x - v) / twoTauE);
        const double a2 = twoWt * sqrtDT2x2 / tm
                          / (std::pow(x * (twoWt - sqrtDT2x2 / DT), 2) + DT2x2 / std::pow(DT * tm, 2));
        i += DT / sqrtDT2x2 * da0a1 * a2;
        x += dx;
    }
    i *= dx;
    return i;
}

template<class ScalarLike>
std::tuple<ScalarLike, ScalarLike> PowerCoolInt(const double DT, const ScalarLike& v,
                                                const double tau, const double tauE) {
    /*
     *  Compute the integral of the cooling power of SIN junction, rectangle method
     *
     *  `DT` is delta (energy gap),
     *  `v` is voltage,
     *  `tau` is for temperature,
     *  `tauE` is for electron temperature
     */

    // E = x, V = v, T = tau, [x] = [v] = [tau]

    const double dx = INTEGRATION_SCALE * DT;
    const double inf = ESSENTIALLY_INFINITY_SCALE * DT;
    const ScalarLike v1 = std::abs(v);
    const double DT2 = std::pow(DT, 2);

    ScalarLike po = 0.0; // power
    ScalarLike ps = 0.0; // SIN cooling power

    double x = DT + dx; // energy
    while (x <= inf) {
        double x2 = std::pow(x, 2);
        double a = std::sqrt(x2 - DT2);
        ScalarLike da2a1;

        da2a1 = (1.0 / (std::exp((x - v1) / tauE) + 1.0) - 1.0 / (std::exp(x / tau) + 1.0)) / a;
        po += x * (x - v1) * da2a1;
        ps += x2 * da2a1;

        da2a1 = (1.0 / (std::exp((-x - v1) / tauE) + 1.0) - 1.0 / (std::exp(-x / tau) + 1.0)) / a;
        po += x * (-x - v1) * da2a1;
        ps += -x2 * da2a1;

        x += dx;
    }
    po *= dx;
    ps *= dx;
    return std::make_tuple(po, ps);
}
#endif // CEB_2EQ_PARALLEL_LITE_H
