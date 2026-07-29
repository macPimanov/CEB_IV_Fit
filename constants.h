#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

/*
 *  Default filename here.
 *  Format: two columns (voltage (V), current(A)).
 *  Decimal delimiter: .
 */
const std::string DATA_FILE_NAME = "SPC-CEB_300mK_Triton11-2026.txt";

// CSV file delimiter
constexpr char SEP = '\t';

/*
 *  Ensure the math constants are defined one way or another
 */
#ifndef M_PI
const double M_PI = std::atan2(1.0, 1.0) * 4.0;
#endif

#ifndef M_SQRT2
const double M_SQRT2 = std::sqrt(2.0);
#endif

#ifndef M_SQRT1_2
const double M_SQRT1_2 = 1.0 / M_SQRT2;
#endif

/*
 *  2018 CODATA recommended values
 *
 *  see https://physics.nist.gov/cuu/Constants/index.html
 */
constexpr double ME = 9.1093837015e-31; // ± 0.0000000028e-31, [kg]
constexpr double E = 1.602176634e-19; // exact,                [C]
constexpr double H = 6.62607015e-34; // exact,                 [J / Hz]
constexpr double HBAR = H / (2.0 * M_PI); //                   [J / Hz]
constexpr double K = 1.380649e-23 / E; // exact,               [eV / K]
// DELTA = 2.78e-23 -> Tc = 1.14
// Rn = 2.24e3 // per 1 SIN
// coefficient follows from BCS theory (TODO: check for which metals)
constexpr double BCS_INTEGRAL = 1.764; // TODO: source? BASE VALUE: 1.764

#define THREADS 56

/*
 *  The scales to use when calculating integrals
 *
 *  The variable to integrate over makes the steps of this size
 *  relative to a characteristic value.
 */
// the step of an integral calculation
constexpr double INTEGRATION_SCALE = 0.2e-3;
// the large enough value
constexpr double ESSENTIALLY_INFINITY_SCALE = 20.0;

constexpr auto NUMBER_OF_SINS_IN_CEB = 2.0;

/*
 *  Voltage and current noise for some amplifiers
 */
#define AD745 745
#define OPA111 111
#define AD797 797
#define IFN146 146
#define OPA1641 1641
#define AMPLIFIER AD745
#ifndef AMPLIFIER
#   error "Undefined amplifier"
#elif AMPLIFIER == AD745
constexpr double VOLTAGE_NOISE = 1.0e-8; // [V/sqrt(Hz)]
constexpr double CURRENT_NOISE = 6.9e-15; // [A/sqrt(Hz)]
#elif AMPLIFIER == OPA111
constexpr double VOLTAGE_NOISE = 8.0e-9; // [V/sqrt(Hz)]
constexpr double CURRENT_NOISE = 0.8e-15; // [A/sqrt(Hz)]
#elif AMPLIFIER == AD797
constexpr double VOLTAGE_NOISE = 0.9e-9; // [V/sqrt(Hz)]
constexpr double CURRENT_NOISE = 2000.0e-15; // [A/sqrt(Hz)]
#elif AMPLIFIER == IFN146
constexpr double VOLTAGE_NOISE = 1.1e-9; // [V/sqrt(Hz)]
constexpr double CURRENT_NOISE = 0.3e-15; // [A/sqrt(Hz)]
#elif AMPLIFIER == OPA1641
constexpr double VOLTAGE_NOISE = 5.1e-9; // [V/sqrt(Hz)]
constexpr double CURRENT_NOISE = 0.8e-15; // [A/sqrt(Hz)]
#else
#   error "Unknown amplifier"
#endif
#undef AD745
#undef OPA111
#undef AD797
#undef IFN146
#undef OPA1641

constexpr double VOLTAGE_NOISE_2_AMPS = VOLTAGE_NOISE * M_SQRT2; // for 2 amp, × sqrt(2)
constexpr double CURRENT_NOISE_2_AMPS = CURRENT_NOISE * M_SQRT1_2; // for 2 amp, / sqrt(2)

#endif
