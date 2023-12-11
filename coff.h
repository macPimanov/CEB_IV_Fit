#ifndef COOF_H
#define COOF_H

#include <valarray>

class COff {
    std::valarray<double> Iofx, Vofx, Iref, Vref;

public:
    COff(const std::valarray<double>& Iexp, const std::valarray<double>& Vexp);

    double operator()(double dOffset);
};

#endif
