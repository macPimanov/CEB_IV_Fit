#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "constants.h"
#include "utils.h"
#include "CEBNumericModel.h"
#include "MinimizationAlgorithms.h"

#include "IVParamFitter.h"

struct IterationInput {
    size_t voltageStep;
    double DeltaT;
    double tauSin;
    double Delta;
    double Vg;
    double Rsin;
    double Rabs;
    double Rleak;
    double I0;
    double Sigma;
    double Vol;
    double TephPOW;
    double Tph;
    double Tc;
    double Wt;
    double tm;
    double ii;
    double dT;
    double dV;
    double bolometersInSeries;
    double bolometersInParallel;
    double beta;
    double dPbg;
    double FREQUENCY;
    double BANDWIDTH;
    double totalBolometersNumber;
    std::valarray<double>* V;
};

struct IterationResult {
    size_t voltageStep;
    double Vnum;
    double Inum;
    double I;
    double I_A;
    double V;
    double Te;
    double Tsin;
    double DeltaT;
    double Pe_ph;
    double Pand;
    double Pleak;
    double Pabs;
    double Pcool;
    double NEPe_ph2;
    double NEPs;
    double NoiA;
    double NEPph;
    double NEP;
    double Sv;
    double G_e;
    double G_NIS;
};

IterationResult computeIteration(const IterationInput& input) {
    IterationResult result{};

    constexpr double dT = 0.005;

    double Pe_ph, Pabs, Pleak, Pcool, Ps, Pand;
    double tauELower = 0.0;
    double tauEUpper = 3.0 / BCS_INTEGRAL;
    double tauE;
    double I, I_A;

    for (size_t l = 0; l < 15; ++l) {
        tauE = (tauELower + tauEUpper) / 2.0;

        I = currentIntegral(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE) * input.I0
            + 1e9 * ((*input.V)[input.voltageStep] / input.Rleak);

        I_A = input.ii * AndCurrent(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, tauE, input.Wt, input.tm) * input.I0;

        Pe_ph = input.Sigma * input.Vol
            * (std::pow(input.Tph, input.TephPOW) - std::pow(tauE * input.Delta, input.TephPOW)) * 1e3;

        Pabs = std::pow(I, 2) * input.Rabs * 1e-6;

        Pleak = NUMBER_OF_SINS_IN_CEB * std::pow((*input.V)[input.voltageStep], 2) / input.Rleak * 1e12;

        Pand = std::pow(I_A * 1e-3, 2) * input.Rabs + 2.0 * (I_A * 1e3) * (*input.V)[input.voltageStep];

        std::tie(Pcool, Ps) = PowerCoolInt(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE);

        Pcool *= std::pow(input.Vg, 2) / input.Rsin * 1e12;
        Ps *= std::pow(input.Vg, 2) / input.Rsin * 1e12;

        const double Pheat = Pe_ph + Pabs + Pand + input.dPbg + 2.0 * input.beta * Ps + Pleak;
        if (Pheat < NUMBER_OF_SINS_IN_CEB * Pcool) {
            tauEUpper = tauE;
        } else {
            tauELower = tauE;
        }
    }

    const double Te = tauE * input.Delta;

    result.Inum = 1e-9 * (I + I_A) * input.bolometersInParallel;
    result.Vnum = (NUMBER_OF_SINS_IN_CEB * (*input.V)[input.voltageStep] + 1e-9 * (I + I_A) * input.Rabs) * input.bolometersInSeries;

    result.Vnum = result.Vnum;
    result.Inum = result.Inum;
    result.I = I;
    result.I_A = I_A;
    result.V = (*input.V)[input.voltageStep];
    result.Te = Te;
    result.Tsin = input.Tph;
    result.DeltaT = input.DeltaT;
    result.Pe_ph = Pe_ph;
    result.Pand = Pand;
    result.Pleak = Pleak;
    result.Pabs = Pabs;
    result.Pcool = Pcool;
    result.voltageStep = input.voltageStep;

    const double dPT = std::get<0>(PowerCoolInt(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE + dT / input.Delta))
        - std::get<0>(PowerCoolInt(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE - dT / input.Delta));

    const double dPdT = 1e12 * (std::pow(input.Vg, 2) / input.Rsin) * dPT / (2.0 * dT);

    const double dIdT = input.I0
        * (currentIntegral(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE + dT / input.Delta)
           - currentIntegral(input.DeltaT, (*input.V)[input.voltageStep] / input.Vg, input.tauSin, tauE - dT / input.Delta))
        / (2.0 * dT);

    const double dIdV = input.I0
        * (currentIntegral(input.DeltaT, (*input.V)[input.voltageStep + 1] / input.Vg, input.tauSin, tauE)
           + input.ii * AndCurrent(input.DeltaT, (*input.V)[input.voltageStep + 1] / input.Vg, tauE, input.Wt, input.tm)
           - currentIntegral(input.DeltaT, (*input.V)[input.voltageStep - 1] / input.Vg, input.tauSin, tauE)
           - input.ii * AndCurrent(input.DeltaT, (*input.V)[input.voltageStep - 1] / input.Vg, tauE, input.Wt, input.tm))
        / (2.0 * input.dV);

    const double dPdV = std::pow(input.Vg, 2) / input.Rsin * 1e12
        * (std::get<0>(PowerCoolInt(input.DeltaT, (*input.V)[input.voltageStep + 1] / input.Vg, input.tauSin, tauE))
           - std::get<0>(PowerCoolInt(input.DeltaT, (*input.V)[input.voltageStep - 1] / input.Vg, input.tauSin, tauE)))
        / (2.0 * input.dV);

    const double G_NIS = dPdT;
    const double G_e = 5.0 * input.Sigma * input.Vol * std::pow(Te, 4) * 1e3;
    const double G = G_e + NUMBER_OF_SINS_IN_CEB * (G_NIS - dIdT / dIdV * dPdV);
    const double Sv = -2.0 * dIdT / dIdV / G / input.bolometersInParallel;

    const double NEPe_ph2 = 10.0 * (E * K) * input.Sigma * input.Vol * (std::pow(input.Tph, input.TephPOW) + std::pow(Te, input.TephPOW)) * 1e3 * 1e12;

    const double NoiA = std::pow(VOLTAGE_NOISE_2_AMPS, 2)
        + std::pow(CURRENT_NOISE_2_AMPS * (2.0 * 1e9 / dIdV + input.Rabs) * input.bolometersInSeries / input.bolometersInParallel, 2);

    const double NEPa = NoiA / std::pow(Sv, 2);

    const double dI = 1e9 * (2.0 * E * std::abs(I) / std::pow(dIdV * Sv, 2));

    const double dPdI = 1e9 * (2.0 * 2.0 * E * Pcool / (dIdV * Sv));

    const double mm = std::log(std::sqrt(2.0 * M_PI * K * Te * input.Vg) / (2.0 * std::abs(I) * input.Rsin * 1e-9));

    const double dP = (0.5 + std::pow(mm, 2)) * std::pow(K * Te, 2) * std::abs(I) * E * 1e-9 * 1e24;

    const double NEPs = NUMBER_OF_SINS_IN_CEB * (dI - 2.0 * dPdI + dP);

    const double NEPph = 1e12 * std::sqrt(input.totalBolometersNumber * 2.0 * (input.FREQUENCY * 1e9) * (input.dPbg * 1e-12) * H + std::pow((input.dPbg * 1e-12) * input.totalBolometersNumber, 2) / (input.BANDWIDTH * 1e9));

    const double NEP = std::sqrt((NEPe_ph2 + NEPs) * input.totalBolometersNumber + NEPa + std::pow(NEPph, 2));

    result.NEPe_ph2 = NEPe_ph2;
    result.NEPs = NEPs;
    result.NoiA = NoiA;
    result.NEPph = NEPph;
    result.NEP = NEP;
    result.Sv = Sv;
    result.G_e = G_e;
    result.G_NIS = G_NIS;

    return result;
}

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

    std::vector<IterationResult> results(voltageStepsCount - 1);
    std::vector<std::thread> threads;
    const size_t numThreads = std::thread::hardware_concurrency();

    IterationInput baseInput{
        .DeltaT = DeltaT,
        .tauSin = tauSin,
        .Delta = Delta,
        .Vg = Vg,
        .Rsin = Rsin,
        .Rabs = Rabs,
        .Rleak = Rleak,
        .I0 = I0,
        .Sigma = Sigma,
        .Vol = Vol,
        .TephPOW = TephPOW,
        .Tph = Tph,
        .Tc = Tc,
        .Wt = Wt,
        .tm = tm,
        .ii = ii,
        .dT = 0.005,
        .dV = dV,
        .bolometersInSeries = bolometersInSeries,
        .bolometersInParallel = bolometersInParallel,
        .beta = beta,
        .dPbg = dPbg,
        .FREQUENCY = FREQUENCY,
        .BANDWIDTH = BANDWIDTH,
        .totalBolometersNumber = totalBolometersNumber,
        .V = &V
    };

    for (size_t i = 0; i < voltageStepsCount - 1; ++i) {
        threads.emplace_back([i, &baseInput, &results]() {
            IterationInput input = baseInput;
            input.voltageStep = i + 1;
            results[i] = computeIteration(input);
        });

        if (threads.size() >= numThreads || i == voltageStepsCount - 2) {
            for (auto& t : threads) {
                if (t.joinable()) {
                    t.join();
                }
            }
            threads.clear();
        }
    }

    std::ranges::sort(results, [](const auto& a, const auto& b) {
        return a.voltageStep < b.voltageStep;
    });

    for (const auto& result : results) {
        I[result.voltageStep] = result.I;
        I_A[result.voltageStep] = result.I_A;

        Inum[result.voltageStep - 1] = result.Inum;
        Vnum[result.voltageStep - 1] = result.Vnum;

        file_Te
                << result.Vnum << SEP
                << result.Inum << SEP
                << 1e-9 * result.I * bolometersInParallel << SEP
                << 1e-9 * result.I_A * bolometersInParallel << SEP
                << 1e9 * (result.V / Rleak) * bolometersInParallel << SEP
                << result.Te << SEP
                << result.Tsin << SEP
                << result.DeltaT << SEP
                << result.Pe_ph << SEP
                << result.Pand << SEP
                << result.Pleak << SEP
                << result.Pabs << SEP
                << result.Pcool << std::endl;

        file_Noise
                << (2.0 * result.V + 1e-9 * result.I * Rabs) * bolometersInSeries << SEP
                << 1e9 * std::sqrt(result.NEPe_ph2 * totalBolometersNumber) * std::abs(result.Sv) << SEP
                << 1e9 * std::sqrt(result.NEPs * totalBolometersNumber) * std::abs(result.Sv) << SEP
                << 1e9 * std::sqrt(result.NoiA) << SEP
                << 1e9 * result.NEP * std::abs(result.Sv) << SEP
                << 1e9 * result.NEPph * std::abs(result.Sv) << SEP
                << 1e9 * std::abs(result.Sv) * std::sqrt(std::pow(result.NEP, 2) - std::pow(result.NEPph, 2)) << std::endl;

        file_NEP
                << (2.0 * result.V + 1e-9 * result.I * Rabs) * bolometersInSeries << SEP
                << 1e-9 * result.I * bolometersInParallel << SEP
                << 1e-12 * std::sqrt(result.NEPe_ph2 * totalBolometersNumber) << SEP
                << 1e-12 * std::sqrt(result.NEPs * totalBolometersNumber) << SEP
                << 1e-12 * std::sqrt(result.NoiA) << SEP
                << 1e-12 * result.NEP << SEP
                << 1e-12 * result.NEPph << SEP
                << 1e12 * std::abs(result.Sv) << SEP
                << 1e-12 * std::sqrt(std::pow(result.NEP, 2) - std::pow(result.NEPph, 2)) << std::endl;

        file_G
                << (NUMBER_OF_SINS_IN_CEB * result.V + 1e-9 * (result.I * Rabs)) * bolometersInSeries << SEP
                << result.G_e << SEP
                << result.G_NIS << std::endl;

        std::clog
                << std::setw(static_cast<int>(std::ceil(std::log10(voltageStepsCount)))) << result.voltageStep << '/' <<
                voltageStepsCount - 1 << ':' << ' '
                << "V:" << std::setw(12) << result.Vnum << SEP
                << "I:" << std::setw(12) << result.Inum << SEP
                << "Sv:" << std::setw(12) << 1e12 * std::abs(result.Sv) << SEP
                << "Te:" << std::setw(12) << result.Te << SEP
                << "NEPs:" << std::setw(12) << 1e-12 * std::sqrt(result.NEPs * totalBolometersNumber) << SEP
                << "NEPt:" << std::setw(12) << 1e-12 * result.NEP << std::endl;
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
