#include "ceb_library.h"
#include <algorithm>
#include <chrono>
#include <cstring>
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

        I_A = AndCurrent(input.ii, input.DeltaT, input.V[input.voltageStep] / input.Vg, tauE, input.Wt, input.tm) * input.I0;

        Pe_ph = input.Sigma * input.Vol
            * (std::pow(input.Tph, input.TephPOW) - std::pow(tauE * input.Delta, input.TephPOW)) * 1e3;

        Pabs = std::pow(I, 2) * input.Rabs * 1e-6;

        Pleak = NUMBER_OF_SINS_IN_CEB * std::pow(input.V[input.voltageStep], 2) / input.Rleak * 1e12;

        Pand = std::pow(I_A * 1e-3, 2) * input.Rabs + 2.0 * (I_A * 1e3) * input.V[input.voltageStep];

        std::tie(Pcool, Ps) = PowerCoolInt(input.DeltaT, input.V[input.voltageStep] / input.Vg, input.tauSin, tauE);

        Pcool *= std::pow(input.Vg, 2) / input.Rsin * 1e12;
        Ps *= std::pow(input.Vg, 2) / input.Rsin * 1e12;

        // std::cout << "VSHAMPOR: dPbg is " << input.dPbg << std::endl;
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
           + AndCurrent(input.ii, input.DeltaT, input.V[input.voltageStep + 1] / input.Vg, tauE, input.Wt, input.tm)
           - currentIntegral(input.DeltaT, input.V[input.voltageStep - 1] / input.Vg, input.tauSin, tauE)
           - AndCurrent(input.ii, input.DeltaT, input.V[input.voltageStep - 1] / input.Vg, tauE, input.Wt, input.tm))
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
        const double dPbg = Pbg; // (... / totalBolometersNumber? )
        const double Delta = (BCS_INTEGRAL * Tc);

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

        // Allocate memory for all arrays
        const size_t resultSize = voltageStepsCount - 1;
        
        // Basic numerical results
        result->Inum.data = new double[resultSize];
        result->Inum.array_size = resultSize;
        result->Vnum.data = new double[resultSize];
        result->Vnum.array_size = resultSize;
        
        // Detailed per-iteration results
        result->I.data = new double[resultSize];
        result->I.array_size = resultSize;
        result->I_A.data = new double[resultSize];
        result->I_A.array_size = resultSize;
        result->Te.data = new double[resultSize];
        result->Te.array_size = resultSize;
        result->Tsin.data = new double[resultSize];
        result->Tsin.array_size = resultSize;
        result->DeltaT.data = new double[resultSize];
        result->DeltaT.array_size = resultSize;
        result->Pe_ph.data = new double[resultSize];
        result->Pe_ph.array_size = resultSize;
        result->Pand.data = new double[resultSize];
        result->Pand.array_size = resultSize;
        result->Pleak.data = new double[resultSize];
        result->Pleak.array_size = resultSize;
        result->Pabs.data = new double[resultSize];
        result->Pabs.array_size = resultSize;
        result->Pcool.data = new double[resultSize];
        result->Pcool.array_size = resultSize;
        result->NEPe_ph2.data = new double[resultSize];
        result->NEPe_ph2.array_size = resultSize;
        result->NEPs.data = new double[resultSize];
        result->NEPs.array_size = resultSize;
        result->NoiA.data = new double[resultSize];
        result->NoiA.array_size = resultSize;
        result->NEPph.data = new double[resultSize];
        result->NEPph.array_size = resultSize;
        result->NEP.data = new double[resultSize];
        result->NEP.array_size = resultSize;
        result->Sv.data = new double[resultSize];
        result->Sv.array_size = resultSize;
        result->G_e.data = new double[resultSize];
        result->G_e.array_size = resultSize;
        result->G_NIS.data = new double[resultSize];
        result->G_NIS.array_size = resultSize;

        std::vector<double> V(voltageStepsCount + 1);
        for (size_t i = 0; i <= voltageStepsCount; ++i) {
            V[i] = Vstr + (static_cast<double>(i) * dV);
        }

        std::vector<CEBLibrary::IterationResult> results(resultSize);
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

        for (size_t i = 0; i < resultSize; ++i) {
            threads.emplace_back([i, &baseInput, &results]() {
                CEBLibrary::IterationInput input = baseInput;
                input.voltageStep = i + 1;
                results[i] = CEBLibrary::computeIteration(input);
            });

            if (threads.size() >= numThreads || i == resultSize - 1) {
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

        // Write console output and populate result arrays
        for (size_t i = 0; i < resultSize; ++i) {
            const auto& res = results[i];
            
            result->Inum.data[i] = res.Inum;
            result->Vnum.data[i] = res.Vnum;
            result->I.data[i] = res.I;
            result->I_A.data[i] = res.I_A;
            result->Te.data[i] = res.Te;
            result->Tsin.data[i] = res.Tsin;
            result->DeltaT.data[i] = res.DeltaT;
            result->Pe_ph.data[i] = res.Pe_ph;
            result->Pand.data[i] = res.Pand;
            result->Pleak.data[i] = res.Pleak;
            result->Pabs.data[i] = res.Pabs;
            result->Pcool.data[i] = res.Pcool;
            result->NEPe_ph2.data[i] = res.NEPe_ph2;
            result->NEPs.data[i] = res.NEPs;
            result->NoiA.data[i] = res.NoiA;
            result->NEPph.data[i] = res.NEPph;
            result->NEP.data[i] = res.NEP;
            result->Sv.data[i] = res.Sv;
            result->G_e.data[i] = res.G_e;
            result->G_NIS.data[i] = res.G_NIS;
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
        delete[] result->I.data;
        delete[] result->I_A.data;
        delete[] result->Te.data;
        delete[] result->Tsin.data;
        delete[] result->DeltaT.data;
        delete[] result->Pe_ph.data;
        delete[] result->Pand.data;
        delete[] result->Pleak.data;
        delete[] result->Pabs.data;
        delete[] result->Pcool.data;
        delete[] result->NEPe_ph2.data;
        delete[] result->NEPs.data;
        delete[] result->NoiA.data;
        delete[] result->NEPph.data;
        delete[] result->NEP.data;
        delete[] result->Sv.data;
        delete[] result->G_e.data;
        delete[] result->G_NIS.data;
        
        result->Inum.data = nullptr;
        result->Inum.array_size = 0;
        result->Vnum.data = nullptr;
        result->Vnum.array_size = 0;
        result->I.data = nullptr;
        result->I.array_size = 0;
        result->I_A.data = nullptr;
        result->I_A.array_size = 0;
        result->Te.data = nullptr;
        result->Te.array_size = 0;
        result->Tsin.data = nullptr;
        result->Tsin.array_size = 0;
        result->DeltaT.data = nullptr;
        result->DeltaT.array_size = 0;
        result->Pe_ph.data = nullptr;
        result->Pe_ph.array_size = 0;
        result->Pand.data = nullptr;
        result->Pand.array_size = 0;
        result->Pleak.data = nullptr;
        result->Pleak.array_size = 0;
        result->Pabs.data = nullptr;
        result->Pabs.array_size = 0;
        result->Pcool.data = nullptr;
        result->Pcool.array_size = 0;
        result->NEPe_ph2.data = nullptr;
        result->NEPe_ph2.array_size = 0;
        result->NEPs.data = nullptr;
        result->NEPs.array_size = 0;
        result->NoiA.data = nullptr;
        result->NoiA.array_size = 0;
        result->NEPph.data = nullptr;
        result->NEPph.array_size = 0;
        result->NEP.data = nullptr;
        result->NEP.array_size = 0;
        result->Sv.data = nullptr;
        result->Sv.array_size = 0;
        result->G_e.data = nullptr;
        result->G_e.array_size = 0;
        result->G_NIS.data = nullptr;
        result->G_NIS.array_size = 0;
    }
}