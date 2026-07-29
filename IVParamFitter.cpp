#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "constants.h"
#include "utils.h"
#include "CEBNumericModel.h"
#include "MinimizationAlgorithms.h"

#include "IVParamFitter.h"

IVParamFitter::IVParamFitter() {
    /*
     *  Set the base parameters for fitting
     */
    if (std::ifstream parfile("startparams.txt"); parfile) {
        std::string parname;
        double parvalue;
        bool parfit;
        while (!parfile.eof()) {
            parfile >> std::skipws >> parname >> parvalue >> parfit;

            par[parname] = parvalue;
            ToFit[parname] = parfit;

            std::clog << std::left << std::setw(18) << std::format("{} = {},", parname, parvalue)
                    << std::internal << "to fit = " << std::boolalpha << parfit
                    << std::endl;
        }
        parfile.close();
    } else {
        throw std::runtime_error("Can't read \"startparams.txt\"");
    }
}

double IVParamFitter::operator()(const double dParam) {
    par[parameterName] = dParam;

    computeCEBProperties();
    auto [Irex, Vrex] = resample();

    return ChiSqDer(Vnum, Inum, Irex);
}

void IVParamFitter::SeqFit(const size_t runCount, const std::valarray<double>& Irex) {
    /*
     *  Fit using Golden method
     */

    std::random_device r;
    std::default_random_engine generator(r());

    writeConverg(par[parameterName], ChiSq(Inum, Irex), std::chrono::steady_clock::now());

    for (size_t run = 0; run < runCount; ++run) {
        std::clog << "SeqFit run " << run << std::endl;

        std::vector<std::string> ParSeq;
        for (const auto& [key, value]: ToFit) {
            if (value) {
                ParSeq.push_back(key);
            }
        }
        std::ranges::shuffle(ParSeq, generator);

        double fmin = NAN;
        for (auto& parname: ParSeq) {
            parameterName = parname;
            std::tie(par[parameterName], fmin) = GoldenMinimize(
                *this,
                0.5 * par[parameterName], 2.0 * par[parameterName],
                1.0 * par[parameterName],
                1e-3
            );
        }

        // store the parameters after minimization and the result
        bool appendNewLine = std::filesystem::exists("fitparameters_new.txt")
                             && std::filesystem::file_size("fitparameters_new.txt");
        if (std::fstream params("fitparameters_new.txt", std::ios::app); params) {
            if (appendNewLine) {
                params << std::endl;
            }
            params << std::format("time = {}", std::chrono::system_clock::now()) << std::endl;
            for (const auto& [parname, parvalue]: par)
                params << std::format("{} = {} ({})", parname, parvalue,
                                      ToFit[parname] ? std::string("fit") : std::string("skip")) << std::endl;
            params << std::format("fmin = {}", fmin) << std::endl;
            params.close();
        } else {
            throw std::runtime_error("Unable to append to \"fitparameters_new.txt\"");
        }
    }
}

size_t IVParamFitter::loadExperimentData(const std::string& filename, const bool removeOffset) {
    std::tie(Iexp, Vexp) = getExperimentalData(filename, removeOffset);
    if (Iexp.size() != Vexp.size()) {
        throw std::length_error("Experimental I and V must be of the same size");
    }
    return Iexp.size();
}

std::tuple<std::valarray<double>, std::valarray<double>> IVParamFitter::resample() const {
    return Resample(Iexp, Vexp, Inum, Vnum);
}

size_t IVParamFitter::computeCEBProperties() {
    /*
     *  The function has equations from DOI: 10.1063/1.1351002
     */
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

    // --------- known/guessed physical parameters

    // `bolometersInSeries` and `bolometersInParallel` may be of an integer type,
    // but as they're used in floating-point operations, let them be `double`
    const auto bolometersInSeries = par["M"]; // number of bolometers in series
    const auto bolometersInParallel = par["MP"]; // number of bolometers in parallel
    const auto totalBolometersNumber = bolometersInSeries * bolometersInParallel;

    // incoming power for all structure [pW]
    const double Pbg = par["Pbg"];
    // returning power ratio, <1
    const double beta = par["beta"];
    // exponent for Te-ph, 7, 6, or 5
    const double TephPOW = par["TephPOW"];
    // volume of the absorber [um³]
    const double Vol = par["Vol"];
    // heat exchange in normal metal [nW/(K⁵×um³)]
    const double Sigma = par["Z"];
    // critical temperature [K]
    const double Tc = par["Tc"];
    // normal resistance for 1 bolometer [Ohm]
    const double Rn = par["Rn"] * bolometersInParallel / bolometersInSeries;
    // leakage resistance per 1 bolometer [Ohm]
    const double Rleak = par["Rleak"] * bolometersInParallel / bolometersInSeries;
    // transparency of the barrier
    const double Wt = par["Wt"];
    // depairing energy
    const double tm = par["tm"];
    // coefficient for Andreev current
    const double ii = par["ii"];
    // normal resistance of 1 absorber [Ohm]
    const double Rabs = par["Ra"];
    // phonon temperature [K]
    const double Tph = par["Tp"];
    // main frequency [GHz]
    const double FREQUENCY = par["F"];
    // bandwidth [GHz]
    const double BANDWIDTH = par["dF"];
    // voltage range end [V]
    const double dVFinVg = par["dVFinVg"];
    // voltage range start [V]
    const double dVStartVg = par["dVStartVg"];
    // voltage step [V]
    const double dV = par["dV"];

    // electron temperature to be found [K]
    double Te = Tph;
    // electron temperature in superconductor [K]
    const double Tsin = Tph;

    const double DeltaT = std::sqrt(1.0 - std::pow(Tsin / Tc, 3.2));
    // incoming power per 1 bolometer [pW]
    const double dPbg = Pbg /*/ totalBolometersNumber*/;
    // energy gap [K], Vg[eV] = Tc * BCS_INTEGRAL * 86.25e-6
    const double Delta = (BCS_INTEGRAL * Tc); // [K]

    // if there is a file named “Te.txt”, backup its content into “Te_old.txt”
    if (std::filesystem::exists("Te.txt")) {
        std::filesystem::rename("Te.txt", "Te_old.txt");
    }

    //---------- normalized constants
    // leave only the resistance of SIN junctions
    const double Rsin = (Rn - Rabs) / NUMBER_OF_SINS_IN_CEB;

    const double I0 = 1e9 * (Delta / Rsin * K); // [nA], units of current

    const double Vg = Delta * K; // [eV]

    const double tauSin = Tsin / Delta; // dimensionless
    double tauE = Te / Delta; // dimensionless

    //---------- calculation parameters

    // initial voltage
    const double Vstr = dVStartVg * Vg;
    // final voltage
    const double Vfin = dVFinVg * Vg;

    const auto voltageStepsCount = static_cast<size_t>(std::round((Vfin - Vstr) / dV)); // the number of voltage steps
    if (!voltageStepsCount) {
        throw std::length_error("No voltage steps to do");
    }

    Inum.resize(voltageStepsCount - 1);
    Vnum.resize(voltageStepsCount - 1);

    std::valarray<double> V(voltageStepsCount + 1); // [V]
    std::ranges::iota(V, 0);
    V = Vstr + (V * dV);

    std::ofstream file_Noise("Noise.txt");
    if (!file_Noise) {
        throw std::runtime_error("Unable to write \"file_Noise.txt\"");
    }
    std::ofstream file_Te("Te.txt");
    if (!file_Te) {
        throw std::runtime_error("Unable to write \"file_Te.txt\"");
    }
    std::ofstream file_NEP("NEP.txt");
    if (!file_NEP) {
        throw std::runtime_error("Unable to write \"NEP.txt\"");
    }
    std::ofstream file_G("G.txt");
    if (!file_G) {
        throw std::runtime_error("Unable to write \"G.txt\"");
    }

    file_Noise
            << "Voltage" << SEP
            << "NOISEep" << SEP
            << "NOISEs" << SEP
            << "NOISEa" << SEP
            << "NOISE" << SEP
            << "NOISEph" << SEP
            << "NOISE^2-NOISEph^2" << std::endl;
    file_Te
            << "Voltage" << SEP
            << "Current" << SEP
            << "Iqp" << SEP
            << "Iand" << SEP
            << "V/Rleak" << SEP
            << "Te" << SEP
            << "Ts" << SEP
            << "DeltaT" << SEP
            << "Peph" << SEP
            << "Pand" << SEP
            << "Pleak" << SEP
            << "Pabs" << SEP
            << "Pcool" << std::endl;
    file_NEP
            << "Voltage" << SEP
            << "Current" << SEP
            << "NEPeph" << SEP
            << "NEPs" << SEP
            << "NEPa" << SEP
            << "NEP" << SEP
            << "NEPph" << SEP
            << "Sv" << SEP
            << "NEP^2-NEPph^2" << std::endl;
    file_G
            << "Voltage" << SEP
            << "Ge" << SEP
            << "Gnis" << std::endl;

    std::valarray<double> I(voltageStepsCount + 1);
    std::valarray<double> I_A(voltageStepsCount + 1);

    for (size_t voltageStep = 1; voltageStep < voltageStepsCount; ++voltageStep) // next voltage; exclude edges
    {
        constexpr double dT = 0.005; // temperature step for derivative calculations

        double Pe_ph, Pabs, Pleak, Pcool, Ps, Pand; // for power

        // for (size_t n = 0; n < 5; ++n) // next interation
        // {
        double tauELower = 0.0; // dimentionless
        double tauEUpper = 3.0 / BCS_INTEGRAL; // dimentionless

        // find tauE so that Pheat == NUMBER_OF_SINS_IN_CEB * Pcool
        for (size_t l = 0; l < 15; ++l)
        // next iteration; TODO: check with Leonid's old version to see how this is different from `n` above
        {
            tauE = (tauELower + tauEUpper) / 2.0;

            I[voltageStep] = currentIntegral(DeltaT, V[voltageStep] / Vg, tauSin, tauE) * I0 + 1e9 * (V[voltageStep] / Rleak); // [nA]

//            I[voltageStep] = current( V[voltageStep] / Vg, tauE) * I0 + 1e9 * (V[voltageStep] / Rleak); // [nA]



            I_A[voltageStep] = ii * AndCurrent(DeltaT, V[voltageStep] / Vg, tauE, Wt, tm) * I0; // [nA]

            Pe_ph = Sigma * Vol
                                 * (std::pow(Tph, TephPOW)
                                    - std::pow(tauE * Delta, TephPOW))
                                 * 1e3; // [pW]

            Pabs = std::pow(I[voltageStep], 2) * Rabs * 1e-6; // [pW]

            Pleak = NUMBER_OF_SINS_IN_CEB * std::pow(V[voltageStep], 2) / Rleak * 1e12; // [pW]

            Pand = std::pow(I_A[voltageStep] * 1e-3 /*[uA]*/, 2) * Rabs /*[Ohm]*/
                   + 2.0/*TODO*/ * (I_A[voltageStep] * 1e3) /*[pA]*/ * V[voltageStep] /*[V]*/;
            // [pW], absorber + Andreev

            std::tie(Pcool, Ps) = PowerCoolInt(DeltaT, V[voltageStep] / Vg, tauSin, tauE);

            Pcool *= std::pow(Vg, 2) / Rsin * 1e12; // [pW]
            Ps *= std::pow(Vg, 2) / Rsin * 1e12; // [pW], returning power from S to N

            if (const double Pheat = Pe_ph + Pabs + Pand + dPbg + 2.0/*TODO*/ * beta * Ps + Pleak;
                Pheat < NUMBER_OF_SINS_IN_CEB * Pcool) {
                tauEUpper = tauE;
            } else {
                tauELower = tauE;
            }
        }
        // }

        Te = tauE * Delta;

        Inum[voltageStep - 1] = 1e-9 * (I[voltageStep] + I_A[voltageStep]) * bolometersInParallel;
        Vnum[voltageStep - 1] = (NUMBER_OF_SINS_IN_CEB * V[voltageStep] + 1e-9 * (I[voltageStep] + I_A[voltageStep]) *
                                 Rabs) * bolometersInSeries;

        file_Te
                << Vnum[voltageStep - 1] << SEP
                << Inum[voltageStep - 1] << SEP
                << 1e-9 * I[voltageStep] * bolometersInParallel << SEP
                << 1e-9 * I_A[voltageStep] * bolometersInParallel << SEP
                << 1e9 * (V[voltageStep] / Rleak) * bolometersInParallel << SEP
                << Te << SEP
                << Tsin << SEP
                << DeltaT << SEP
                << Pe_ph << SEP
                << Pand << SEP
                << Pleak << SEP
                << Pabs << SEP
                << Pcool << std::endl;

        //----- NEP ----------------------------------------------

        const double dPT = std::get<0>(PowerCoolInt(DeltaT, V[voltageStep] / Vg, tauSin, tauE + dT / Delta))
                           - std::get<0>(PowerCoolInt(DeltaT, V[voltageStep] / Vg, tauSin, tauE - dT / Delta));

        const double dPdT = 1e12 * (std::pow(Vg, 2) / Rsin) * dPT / (2.0 * dT); // [pW/K]

        const double dIdT = I0
                            * (currentIntegral(DeltaT, V[voltageStep] / Vg, tauSin, tauE + dT / Delta)
                               - currentIntegral(DeltaT, V[voltageStep] / Vg, tauSin, tauE - dT / Delta))
                            / (2.0 * dT); // [nA/K]

//        const double dIdT = I0
//                            * (current( V[voltageStep] / Vg, tauE + dT / Delta)
//                               - current(V[voltageStep] / Vg, tauE - dT / Delta))
//                            / (2.0 * dT); // [nA/K]

        const double dIdV = I0
                            * (currentIntegral(DeltaT, V[voltageStep + 1] / Vg, tauSin, tauE)
                               + ii * AndCurrent(DeltaT, V[voltageStep + 1] / Vg, tauE, Wt, tm)
                               - currentIntegral(DeltaT, V[voltageStep - 1] / Vg, tauSin, tauE)
                               - ii * AndCurrent(DeltaT, V[voltageStep - 1] / Vg, tauE, Wt, tm))
                            / (2.0 * dV); // [nA/V]

//        const double dIdV = I0
//                            * (current(V[voltageStep + 1] / Vg, tauE)
//                               + ii * AndCurrent(DeltaT, V[voltageStep + 1] / Vg, tauE, Wt, tm)
//                               - current(V[voltageStep - 1] / Vg, tauE)
//                               - ii * AndCurrent(DeltaT, V[voltageStep - 1] / Vg, tauE, Wt, tm))
//                            / (2.0 * dV); // [nA/V]

        const double dPdV = std::pow(Vg, 2) / Rsin * 1e12
                            * (std::get<0>(PowerCoolInt(DeltaT, V[voltageStep + 1] / Vg, tauSin, tauE))
                               - std::get<0>(PowerCoolInt(DeltaT, V[voltageStep - 1] / Vg, tauSin, tauE)))
                            / (2.0 * dV); // [pW/V]

        // heat conductance
        const double G_NIS = dPdT; // after eq. (10)
        const double G_e = 5.0 * Sigma * Vol * std::pow(Te, 4) * 1e3; // [pW/K], after eq. (10)

        const double G = G_e + NUMBER_OF_SINS_IN_CEB * (G_NIS - dIdT / dIdV * dPdV); // [pW/K]

        // the responsivity in the current biased regime, eq. (30)
        const double Sv = -2.0 * dIdT / dIdV / G / bolometersInParallel; // [V/pW], for 1 bolo

        // NEPe_ph squared, eq. (24)
        const double NEPe_ph2 = 10.0 * (E * K) * Sigma * Vol * (std::pow(Tph, TephPOW) + std::pow(Te, TephPOW)) * 1e3
                                * 1e12; // [pW²/Hz]

        // amplifier noise
        const double NoiA = std::pow(VOLTAGE_NOISE_2_AMPS, 2)
                            + std::pow(
                                CURRENT_NOISE_2_AMPS * (2.0 * 1e9 / dIdV + Rabs) * bolometersInSeries /
                                bolometersInParallel,
                                2); // [V²/Hz]

        const double NEPa = NoiA / std::pow(Sv, 2); // [pW²/Hz]

        //----- NEP SIN approximation ----------------------------

        const double dI = 1e9 * (2.0 * E * std::abs(I[voltageStep]) / std::pow(dIdV * Sv, 2)); // [pW²/Hz]

        const double dPdI = 1e9 * (2.0 * 2.0 * E * Pcool / (dIdV * Sv));
        // [pW²/Hz], second '2' is from comparison with integral

        const double mm =
                std::log(std::sqrt(2.0 * M_PI * K * Te * Vg) / (2.0 * std::abs(I[voltageStep]) * Rsin * 1e-9));

        const double dP = (0.5 + std::pow(mm, 2)) * std::pow(K * Te, 2) * std::abs(I[voltageStep]) * E * 1e-9
                          * 1e24; // [pW²/Hz]

        const double NEPs = NUMBER_OF_SINS_IN_CEB * (dI - 2.0 * dPdI + dP); // [pW²/Hz], all terms positive

        //----- NEP SIN integral ---------------------------------
        /*
		mm = NEPInt(DeltaT, V[j] / Vg, tau, tauE, &dI, &dP, &dPdI);

		dI = EL * I0 * dI / std::pow(dIdV * Sv, 2) * 1e9;			// [pW²/Hz]

		dPdI = dPdI / (dIdV * Sv) * EL * std::pow(Vg, 2) / Rsin * 1e12 * 1e9;	// [pW²/Hz]

		dP *= std::pow(Vg, 3) / Rsin * EL * 1e24;			// [pW²/Hz]

		NEPs = NUMBER_OF_SINS_IN_CEB * (dI - 2.0 * dPdI + dP);				// [pW²/Hz], all terms positive
        */
        //--------------------------------------------------------

        // const double NEPph = 1e-6 * (Pbg * totalBolometersNumber);
        // [pW/sqrt(Hz)], at 0 GHz

        const double NEPph = 1e12 * std::sqrt(totalBolometersNumber * 2.0 * (FREQUENCY * 1e9) * (dPbg * 1e-12) * H + std::pow((dPbg * 1e-12) * totalBolometersNumber, 2) / (BANDWIDTH * 1e9));	// [pW/sqrt(Hz)]

        const double NEP = std::sqrt((NEPe_ph2 + NEPs) * totalBolometersNumber + NEPa + std::pow(NEPph, 2));
        //all squares

        file_Noise
                << (2.0 * V[voltageStep] + 1e-9 * I[voltageStep] * Rabs) * bolometersInSeries << SEP
                << 1e9 * std::sqrt(NEPe_ph2 * totalBolometersNumber) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NEPs * totalBolometersNumber) * std::abs(Sv) << SEP
                << 1e9 * std::sqrt(NoiA) << SEP
                << 1e9 * NEP * std::abs(Sv) << SEP
                << 1e9 * NEPph * std::abs(Sv) << SEP
                << 1e9 * std::abs(Sv) * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)) << std::endl;

        file_NEP
                << (2.0 * V[voltageStep] + 1e-9 * I[voltageStep] * Rabs) * bolometersInSeries << SEP
                << 1e-9 * I[voltageStep] * bolometersInParallel << SEP
                << 1e-12 * std::sqrt(NEPe_ph2 * totalBolometersNumber) << SEP
                << 1e-12 * std::sqrt(NEPs * totalBolometersNumber) << SEP
                << 1e-12 * std::sqrt(NEPa) << SEP
                << 1e-12 * NEP << SEP
                << 1e-12 * NEPph << SEP
                << 1e12 * std::abs(Sv) << SEP
                << 1e-12 * std::sqrt(std::pow(NEP, 2) - std::pow(NEPph, 2)) << std::endl;

        file_G
                << (NUMBER_OF_SINS_IN_CEB * V[voltageStep] + 1e-9 * (I[voltageStep] * Rabs)) * bolometersInSeries << SEP
                << G_e << SEP
                << G_NIS << std::endl;

        std::clog
                << std::setw(static_cast<int>(std::ceil(std::log10(voltageStepsCount)))) << voltageStep << '/' <<
                voltageStepsCount - 1 << ':' << ' '
                << "V:" << std::setw(12) << Vnum[voltageStep - 1] << SEP
                << "I:" << std::setw(12) << Inum[voltageStep - 1] << SEP
                << "Sv:" << std::setw(12) << 1e12 * std::abs(Sv) << SEP
                << "Te:" << std::setw(12) << Te << SEP
                << "NEPs:" << std::setw(12) << 1e-12 * std::sqrt(NEPs * totalBolometersNumber) << SEP
                << "NEPt:" << std::setw(12) << 1e-12 * NEP << std::endl;
    }

    file_Noise.close();
    file_Te.close();
    file_NEP.close();
    file_G.close();

    std::clog
            << "Time spent: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            << std::endl;

    return voltageStepsCount - 1;
}
