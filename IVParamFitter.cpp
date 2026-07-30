#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "constants.h"
#include "utils.h"
#include "CEBNumericModel.h"
#include "MinimizationAlgorithms.h"
#include "ceb_library.h"

#include "IVParamFitter.h"

using CEBLibrary::IterationInput;
using CEBLibrary::IterationResult;
using CEBLibrary::computeIteration;

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

            std::ostringstream oss;
            oss << std::left << std::setw(18) << parname << " = " << parvalue << ",";
            std::clog << oss.str()
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
            std::ostringstream oss;
            oss << "time = " << std::chrono::system_clock::now();
            params << oss.str() << std::endl;
            for (const auto& [parname, parvalue]: par) {
                std::ostringstream oss2;
                oss2 << parname << " = " << parvalue << " (" << (ToFit[parname] ? std::string("fit") : std::string("skip")) << ")";
                params << oss2.str() << std::endl;
            }
            std::ostringstream oss3;
            oss3 << "fmin = " << fmin;
            params << oss3.str() << std::endl;
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

    std::vector<double> V(voltageStepsCount + 1);
    for (size_t i = 0; i < V.size(); ++i) {
        V[i] = Vstr + (static_cast<double>(i) * dV);
    }

    // Call shared library's threaded implementation
    CEBParameters params{
        .M = bolometersInSeries,
        .MP = bolometersInParallel,
        .Pbg = Pbg,
        .beta = beta,
        .TephPOW = TephPOW,
        .Vol = Vol,
        .Z = Sigma,
        .Tc = Tc,
        .Rn = par["Rn"],
        .Rleak = Rleak * bolometersInSeries / bolometersInParallel,
        .Wt = Wt,
        .tm = tm,
        .ii = ii,
        .Ra = Rabs,
        .Tp = Tph,
        .F = FREQUENCY,
        .dF = BANDWIDTH,
        .dVFinVg = dVFinVg,
        .dVStartVg = dVStartVg,
        .dV = dV,
        .voltage_noise = 0.0,
        .current_noise = 0.0
    };

    CEBResult resultStruct{};
    compute_ceb_properties_threaded(&params, &resultStruct);

    if (resultStruct.error_code != 0) {
        throw std::runtime_error(resultStruct.error_message);
    }

    // Copy results from CEBResult to member variables
    for (size_t i = 0; i < resultStruct.Inum.array_size; ++i) {
        Inum[i] = resultStruct.Inum.data[i];
        Vnum[i] = resultStruct.Vnum.data[i];
    }

    std::clog
            << "Time spent: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            << std::endl;

    return voltageStepsCount - 1;
}
