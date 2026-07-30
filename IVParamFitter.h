#ifndef CFOO_H
#define CFOO_H

#include "ceb_library.h"
#include <string>
#include <unordered_map>
#include <valarray>

class IVParamFitter {
    std::valarray<double> Iexp, Vexp;
    std::valarray<double> Inum, Vnum;

    std::unordered_map<std::string, double> par;
    std::unordered_map<std::string, bool> ToFit;
    std::string parameterName;
    
    // Store detailed computation results temporarily for file writing
    std::vector<double> detailed_I, detailed_I_A, detailed_Te, detailed_Tsin, detailed_DeltaT;
    std::vector<double> detailed_Pe_ph, detailed_Pand, detailed_Pleak, detailed_Pabs, detailed_Pcool;
    std::vector<double> detailed_NEPe_ph2, detailed_NEPs, detailed_NoiseA, detailed_NEPph, detailed_NEP;
    std::vector<double> detailed_Sv, detailed_G_e, detailed_G_NIS;
    std::string output_dir;

public:
    double operator()(double dParam);

    size_t computeCEBProperties();

    void SeqFit(size_t runCount, const std::valarray<double>& Irex);

    size_t loadExperimentData(const std::string& filename, bool removeOffset = false);

    [[nodiscard]] std::tuple<std::valarray<double>, std::valarray<double>> resample() const;

    void set_output_dir(const std::string& output_dir);
    
    explicit IVParamFitter();
private:
    void writeOutputFiles(const CEBResult& result, const double bolometersInSeries, const double bolometersInParallel);
};

#endif
