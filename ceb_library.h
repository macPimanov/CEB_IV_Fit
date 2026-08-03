#ifndef CEB_LIBRARY_H
#define CEB_LIBRARY_H

#include <cstddef>

// C++ internal structures and functions (not part of C interface)
#ifdef __cplusplus
namespace CEBLibrary {

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
    double voltageNoise;
    double currentNoise;
    const double* V;
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

IterationResult computeIteration(const IterationInput& input);

} // namespace CEBLibrary
#endif // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif

// C-compatible interface for CEB computation library

typedef struct {
    double* data;
    size_t array_size;
} DoubleArray;

typedef struct {
    // Physical parameters
    double M;          // number of bolometers in series
    double MP;         // number of bolometers in parallel
    double Pbg;        // incoming power [pW]
    double beta;       // returning power ratio
    double TephPOW;    // exponent for Te-ph, 7, 6, or 5
    double Vol;        // volume of absorber [um³]
    double Z;         // heat exchange in normal metal
    double Tc;         // critical temperature [K]
    double Rn;        // normal resistance for 1 bolometer [Ohm]
    double Rleak;     // leakage resistance per 1 bolometer [Ohm]
    double Wt;        // transparency of the barrier
    double tm;        // depairing energy
    double ii;        // coefficient for Andreev current
    double Ra;        // normal resistance of 1 absorber [Ohm]
    double Tp;        // phonon temperature [K]
    double F;          // main frequency [GHz]
    double dF;        // bandwidth [GHz]
    double dVFinVg;   // voltage range end [Vg units]
    double dVStartVg; // voltage range start [Vg units]
    double dV;        // voltage step [V]
    
    // Amplifier noise parameters
    double voltage_noise;  // [V/sqrt(Hz)]
    double current_noise;  // [A/sqrt(Hz)]
} CEBParameters;

typedef struct {
    // Basic numerical results
    DoubleArray Inum;  // Numerical current values
    DoubleArray Vnum;  // Numerical voltage values
    
    // Detailed per-iteration results (optional, may be NULL)
    DoubleArray I;           // Current values per step
    DoubleArray I_A;         // Andreev current values per step
    DoubleArray Te;           // Electron temperature per step
    DoubleArray Tsin;         // Superconductor temperature per step
    DoubleArray DeltaT;       // DeltaT per step
    DoubleArray Pe_ph;        // Electron-phonon power per step
    DoubleArray Pand;         // Andreev power per step
    DoubleArray Pleak;        // Leakage power per step
    DoubleArray Pabs;         // Absorbed power per step
    DoubleArray Pcool;        // Cooling power per step
    DoubleArray NEPe_ph2;     // Electron-phonon NEP squared per step
    DoubleArray NEPs;         // SIN NEP squared per step
    DoubleArray NoiA;         // Amplifier noise squared per step
    DoubleArray NEPph;        // Photon NEP per step
    DoubleArray NEP;          // Total NEP per step
    DoubleArray Sv;           // Responsivity per step
    DoubleArray G_e;          // Electron thermal conductance per step
    DoubleArray G_NIS;        // NIS thermal conductance per step
    
    double time_spent; // Computation time in seconds
    int error_code;    // 0 for success, non-zero for error
    char error_message[256];
} CEBResult;

// Main computation function
void compute_ceb_properties_threaded(const CEBParameters* params, CEBResult* result);

// Memory management utilities
void free_ceb_result(CEBResult* result);

#ifdef __cplusplus
}
#endif

#endif // CEB_LIBRARY_H