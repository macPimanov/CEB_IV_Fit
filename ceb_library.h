#ifndef CEB_LIBRARY_H
#define CEB_LIBRARY_H

#include <cstddef>

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
    double Vol;        // volume of the absorber [um³]
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
    DoubleArray Inum;  // Numerical current values
    DoubleArray Vnum;  // Numerical voltage values
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