#include "ceb_library.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

#include "constants.h"
#include "CEBNumericModel.h"

CEBLibrary::IterationResult CEBLibrary::computeIteration(const CEBLibrary::IterationInput& input) {
    CEBLibrary::IterationResult result{};

    double Pe_ph, Pabs, Pleak, Pcool, Ps, Pand;
    double tauELower = 0.0;
    double tauEUpper = 3.0 / BCS_INTEGRAL;
    double tauE;
    double I, I_A;

    for (size_t l = 0; l < 15; ++l) {
        tauE = (tauELower + tauEUpper) / 2.0;

        I = currentIntegral(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE) * input.I0
            + 1e9 * (input.V[input.voltageStep] / input.Rleak);

        I_A = input.ii * AndCurrent(input.DeltaT, input.V[input.voltageStep] / input.Vg, tauE, input.Wt, input.tm) * input.I0;

        Pe_ph = input.Sigma * input.Vol
            * (std::pow(input.Tph, input.TephPOW) - std::pow(tauE * input.Delta, input.TephPOW)) * 1e3;

        Pabs = std::pow(I, 2) * input.Rabs * 1e-6;

        Pleak = NUMBER_OF_SINS_IN_CEB * std::pow(input.V[input.voltageStep], 2) / input.Rleak * 1e12;

        Pand = std::pow(I_A * 1e-3, 2) * input.Rabs + 2.0 * (I_A * 1e3) * input.V[input.voltageStep];

        std::tie(Pcool, Ps) = PowerCoolInt(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE);

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
    result.Vnum = (NUMBER_OF_SINS_IN_CEB * input.V[input.voltageStep] + 1e-9 * (I + I_A) * input.Rabs) * input.bolometersInSeries;

    result.Vnum = result.Vnum;
    result.Inum = result.Inum;
    result.I = I;
    result.I_A = I_A;
    result.V = input.V[input.voltageStep];
    result.Te = Te;
    result.Tsin = input.Tph;
    result.DeltaT = input.DeltaT;
    result.Pe_ph = Pe_ph;
    result.Pand = Pand;
    result.Pleak = Pleak;
    result.Pabs = Pabs;
    result.Pcool = Pcool;
    result.voltageStep = input.voltageStep;

    const double dPT = std::get<0>(PowerCoolInt(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE + input.dT / input.Delta))
        - std::get<0>(PowerCoolInt(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE - input.dT / input.Delta));

    const double dPdT = 1e12 * (std::pow(input.Vg, 2) / input.Rsin) * dPT / (2.0 * input.dT);

    const double dIdT = input.I0
        * (currentIntegral(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE + input.dT / input.Delta)
           - currentIntegral(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE - input.dT / input.Delta))
        / (2.0 * input.dT);

    const double dIdV = input.I0
        * (currentIntegral(input.DeltaT, input.V[input.voltageStep + 1] / input.Vg, input.tauSin, tauE)
           + input.ii * AndCurrent(input.DeltaT, input.V[input.voltageStep + 1] / input.Vg, tauE, input.Wt, input.tm)
           - currentIntegral(input.DeltaT, input.V[input.voltageStep - 1] / input.Vg, input.tauSin, tauE)
           - input.ii * AndCurrent(input.DeltaT, input.V[input.voltageStep - 1] / input.Vg, tauE, input.Wt, input.tm))
        / (2.0 * input.dV);

    const double dPdV = std::pow(input.Vg, 2) / input.Rsin * 1e12
        * (std::get<0>(PowerCoolInt(input.DeltaT, input.V[input.voltageStep + 1] / input.Vg, input.tauSin, tauE))
           - std::get<0>(PowerCoolInt(input.DeltaT, input.V[input.voltageStep - 1] / input.Vg, input.tauSin, tauE)))
        / (2.0 * input.dV);

    const double G_NIS = dPdT;
    const double G_e = 5.0 * input.Sigma * input.Vol * std::pow(Te, 4) * 1e3;
    const double G = G_e + NUMBER_OF_SINS_IN_CEB * (G_NIS - dIdT / dIdV * dPdV);
    const double Sv = -2.0 * dIdT / dIdV / G / input.bolometersInParallel;

    const double NEPe_ph2 = 10.0 * (E * K) * input.Sigma * input.Vol * (std::pow(input.Tph, input.TephPOW) + std::pow(Te, input.TephPOW)) * 1e3 * 1e12;

    const double MSQRT2 = std::sqrt(2.0);
    const double MSQRT1_2 = 1.0 / MSQRT2;

    const double NoiA = std::pow(input.voltageNoise * MSQRT2, 2)
        + std::pow(input.currentNoise * MSQRT1_2 * (2.0 * 1e9 / dIdV + input.Rabs) * input.bolometersInSeries / input.bolometersInParallel, 2);

    const double NEPa = NoiA / std::pow(Sv, 2);

    const double MPI = 3.14159265358979323846;
    const double dI = 1e9 * (2.0 * E * std::abs(I) / std::pow(dIdV * Sv, 2));

    const double dPdI = 1e9 * (2.0 * 2.0 * E * Pcool / (dIdV * Sv));

    const double mm = std::log(std::sqrt(2.0 * MPI * K * Te * input.Vg) / (2.0 * std::abs(I) * input.Rsin * 1e-9));

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

void compute_ceb_properties_threaded(const CEBParameters* params, CEBResult* result) {
    std::memset(result, 0, sizeof(CEBResult));
    result->error_code = 0;

    try {
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        const auto bolometersInSeries = params->M;
        const auto bolometersInParallel = params->MP;
        const auto totalBolometersNumber = bolometersInSeries * bolometersInParallel;

        const double Pbg = params->Pbg;
        const double beta = params->beta;
        const double TephPOW = params->TephPOW;
        const double Vol = params->Vol;
        const double Sigma = params->Z;
        const double Tc = params->Tc;
        const double Rn = params->Rn * bolometersInParallel / bolometersInSeries;
        const double Rleak = params->Rleak * bolometersInParallel / bolometersInSeries;
        const double Wt = params->Wt;
        const double tm = params->tm;
        const double ii = params->ii;
        const double Rabs = params->Ra;
        const double Tph = params->Tp;
        const double FREQUENCY = params->F;
        const double BANDWIDTH = params->dF;
        const double dVFinVg = params->dVFinVg;
        const double dVStartVg = params->dVStartVg;
        const double dV = params->dV;
        const double Tsin = Tph;

        const double DeltaT = std::sqrt(1.0 - std::pow(Tsin / Tc, 3.2));
        const double dPbg = Pbg;
        const double Delta = (BCS_INTEGRAL * Tc);

        // Backup existing Te.txt if it exists
        if (std::filesystem::exists("Te.txt")) {
            std::filesystem::rename("Te.txt", "Te_old.txt");
        }

        const double Rsin = (Rn - Rabs) / NUMBER_OF_SINS_IN_CEB;
        const double I0 = 1e9 * (Delta / Rsin * K);
        const double Vg = Delta * K;
        const double tauSin = Tsin / Delta;

        const double Vstr = dVStartVg * Vg;
        const double Vfin = dVFinVg * Vg;

        const auto voltageStepsCount = static_cast<size_t>(std::round((Vfin - Vstr) / dV));
        if (!voltageStepsCount) {
            std::strncpy(result->error_message, "No voltage steps to do", sizeof(result->error_message) - 1);
            result->error_code = 1;
            return;
        }

        // Allocate memory for arrays
        std::vector<double> V(voltageStepsCount + 1);
        for (size_t i = 0; i <= voltageStepsCount; ++i) {
            V[i] = Vstr + (static_cast<double>(i) * dV);
        }

        std::vector<CEBLibrary::IterationResult> results(voltageStepsCount - 1);
        const size_t numThreads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;

        CEBLibrary::IterationInput baseInput{
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
            .voltageNoise = params->voltage_noise,
            .currentNoise = params->current_noise,
            .V = V.data()
        };

        for (size_t i = 0; i < voltageStepsCount - 1; ++i) {
            threads.emplace_back([i, &baseInput, &results]() {
                CEBLibrary::IterationInput input = baseInput;
                input.voltageStep = i + 1;
                results[i] = CEBLibrary::computeIteration(input);
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

        // Write output files
        std::ofstream file_Noise("Noise.txt");
        std::ofstream file_Te("Te.txt");
        std::ofstream file_NEP("NEP.txt");
        std::ofstream file_G("G.txt");

        file_Noise << "Voltage\tNOISEep\tNOISEs\tNOISEa\tNOISE\tNOISEph\tNOISE^2-NOISEph^2\n";
        file_Te << "Voltage\tCurrent\tIqp\tIand\tV/Rleak\tTe\tTs\tDeltaT\tPeph\tPand\tPleak\tPabs\tPcool\n";
        file_NEP << "Voltage\tCurrent\tNEPeph\tNEPs\tNEPa\tNEP\tNEPph\tSv\tNEP^2-NEPph^2\n";
        file_G << "Voltage\tGe\tGnis\n";

        for (const auto& res : results) {
            file_Te << res.Vnum << '\t' << res.Inum << '\t' << 1e-9 * res.I * bolometersInParallel << '\t'
                    << 1e-9 * res.I_A * bolometersInParallel << '\t' << 1e9 * (res.V / Rleak) * bolometersInParallel << '\t'
                    << res.Te << '\t' << res.Tsin << '\t' << res.DeltaT << '\t' << res.Pe_ph << '\t' << res.Pand << '\t'
                    << res.Pleak << '\t' << res.Pabs << '\t' << res.Pcool << '\n';

            file_Noise << (2.0 * res.V + 1e-9 * res.I * Rabs) * bolometersInSeries << '\t'
                      << 1e9 * std::sqrt(res.NEPe_ph2 * totalBolometersNumber) * std::abs(res.Sv) << '\t'
                      << 1e9 * std::sqrt(res.NEPs * totalBolometersNumber) * std::abs(res.Sv) << '\t'
                      << 1e9 * std::sqrt(res.NoiA) << '\t' << 1e9 * res.NEP * std::abs(res.Sv) << '\t'
                      << 1e9 * res.NEPph * std::abs(res.Sv) << '\t'
                      << 1e9 * std::abs(res.Sv) * std::sqrt(std::pow(res.NEP, 2) - std::pow(res.NEPph, 2)) << '\n';

            file_NEP << (2.0 * res.V + 1e-9 * res.I * Rabs) * bolometersInSeries << '\t'
                     << 1e-9 * res.I * bolometersInParallel << '\t' << 1e-12 * std::sqrt(res.NEPe_ph2 * totalBolometersNumber) << '\t'
                     << 1e-12 * std::sqrt(res.NEPs * totalBolometersNumber) << '\t' << 1e-12 * std::sqrt(res.NoiA) << '\t'
                     << 1e-12 * res.NEP << '\t' << 1e-12 * res.NEPph << '\t' << 1e12 * std::abs(res.Sv) << '\t'
                     << 1e-12 * std::sqrt(std::pow(res.NEP, 2) - std::pow(res.NEPph, 2)) << '\n';

            file_G << (NUMBER_OF_SINS_IN_CEB * res.V + 1e-9 * (res.I * Rabs)) * bolometersInSeries << '\t'
                   << res.G_e << '\t' << res.G_NIS << '\n';

            std::clog << std::setw(3) << res.voltageStep << '/' << voltageStepsCount - 1 << ": "
                      << "V:" << std::setw(12) << res.Vnum << "\t"
                      << "I:" << std::setw(12) << res.Inum << "\t"
                      << "Sv:" << std::setw(12) << 1e12 * std::abs(res.Sv) << "\t"
                      << "Te:" << std::setw(12) << res.Te << "\t"
                      << "NEPs:" << std::setw(12) << 1e-12 * std::sqrt(res.NEPs * totalBolometersNumber) << "\t"
                      << "NEPt:" << std::setw(12) << 1e-12 * res.NEP << '\n';
        }

        file_Noise.close();
        file_Te.close();
        file_NEP.close();
        file_G.close();

        // Allocate memory for results
        const size_t resultSize = voltageStepsCount - 1;
        result->Inum.data = new double[resultSize];
        result->Inum.array_size = resultSize;
        result->Vnum.data = new double[resultSize];
        result->Vnum.array_size = resultSize;

        for (size_t i = 0; i < resultSize; ++i) {
            result->Inum.data[i] = results[i].Inum;
            result->Vnum.data[i] = results[i].Vnum;
        }

        result->time_spent = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    } catch (const std::exception& e) {
        std::strncpy(result->error_message, e.what(), sizeof(result->error_message) - 1);
        result->error_code = 1;
    }
}

void free_ceb_result(CEBResult* result) {
    if (result) {
        delete[] result->Inum.data;
        delete[] result->Vnum.data;
        result->Inum.data = nullptr;
        result->Inum.array_size = 0;
        result->Vnum.data = nullptr;
        result->Vnum.array_size = 0;
    }
}