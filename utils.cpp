#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "coff.h"
#include "mins.h"
#include "utils.h"

size_t GetExp(const std::string& fname,
              std::valarray<double>& Iexp,
              std::valarray<double>& Vexp,
              const bool bRemOffset) {
    /*
     *  Get experimental values `Iexp` and `Vexp` from file `fname`
     *
     *  'bRemOffset' indicates to remove offsets
     */

    std::vector<double> vVexp;
    std::vector<double> vIexp;
    try {
        std::ifstream in;
        in.open(fname, std::ifstream::in);

        double tmpV, tmpI;
        while (!in.eof()) {
            in >> std::skipws >> tmpV >> tmpI;
            vVexp.push_back(tmpV);
            vIexp.push_back(tmpI);
        }
    } catch (std::exception& e) {
        std::cerr << "Error while reading experimental data: " << e.what() << std::endl;
        getchar();
        exit(0);
    }

    std::clog << "Fitting for " << fname << " started!" << std::endl;

    std::ofstream params;
    params.open("fitparameters_new.txt", std::ofstream::app);
    params << "Filename: " << fname << std::endl;
    params.close();

    size_t count = std::min(vVexp.size(), vIexp.size());
    Vexp = std::valarray<double>(vVexp.data(), count);
    Iexp = std::valarray<double>(vIexp.data(), count);

    if (bRemOffset) {
        double dOffset = 0;
        double dLBound = -0.0005;
        double dRBound = 0.0005;
        Golden broff;
        broff.ax = dLBound;
        broff.bx = dOffset;
        broff.cx = dRBound;
        COff Off(Iexp, Vexp);
        dOffset = broff.minimize(Off);
        Vexp -= dOffset;

        writeIV("IVoffset.txt", Iexp, Vexp);
    }

    return count;
}

void Resample(std::valarray<double>& Irex,
              std::valarray<double>& Vrex,
              const std::valarray<double>& Iexp,
              const std::valarray<double>& Vexp,
              const std::valarray<double>& Inum,
              const std::valarray<double>& Vnum) {
    /*
     *  Resample the IV-curve to the specified range
     *
     *  Adapt of experimental data for the fitting algorithm.
     *  Make the points differ by the same voltage.
     */

    if (Iexp.size() != Vexp.size()) {
        throw std::length_error("Experimental I and V must be of the same size");
    }
    if (Inum.size() != Vnum.size()) {
        throw std::length_error("Numeric I and V must be of the same size");
    }

    const size_t countnum = Vnum.size();
    const size_t countexp = Vexp.size();
    Vrex = Vnum;
    Irex.resize(countnum);

    size_t lLeft, lRight;
    for (size_t i = 0; i < countnum; ++i) {
        const double dCurVal = Vnum[i];
        size_t lCurPos = 0;
        double tmp = std::abs(dCurVal - Vexp[0]);
        for (size_t j = 1; j < countexp; j++) {
            if (double tmp1 = std::abs(dCurVal - Vexp[j]); tmp1 < tmp) {
                tmp = tmp1;
                lCurPos = j;
            }
        }
        if ((lCurPos < countexp - 1)
            && ((Vexp[lCurPos] <= dCurVal) == (dCurVal <= Vexp[lCurPos + 1]))) {
            lLeft = lCurPos;
            lRight = lCurPos + 1;
        } else {
            lLeft = lCurPos - 1;
            lRight = lCurPos;
        }

        if (lRight == 0)
            Irex[i] = Iexp[0];
        else if (lRight > countexp - 1)
            Irex[i] = Iexp[countexp - 1];
        else
            Irex[i] = Iexp[lLeft]
                      + (dCurVal - Vexp[lLeft]) * (Iexp[lRight] - Iexp[lLeft])
                      / (Vexp[lRight] - Vexp[lLeft]);
    }

    writeIV("IVresampled.txt", Irex, Vrex);
}

double ChiSq(const std::valarray<double>& Inum, const std::valarray<double>& Irex) {
    /*
     *  Calculate Chi² figure of merit for current difference
     */

    if (Inum.size() != Irex.size()) {
        throw std::length_error("Numeric I and Recalculated I must be of the same size");
    }

    return std::pow((Inum - Irex) / Irex, 2).sum() / static_cast<double>(Inum.size());
}

double ChiSqHi(const std::valarray<double>& Inum, const std::valarray<double>& Irex) {
    /*
     *  Calculate Chi² figure of merit
     */

    if (Inum.size() != Irex.size()) {
        throw std::length_error("Numeric I and Recalculated I must be of the same size");
    }

    const std::valarray<double> dI2 = std::pow(Inum - Irex, 2);
    const size_t countnum = Inum.size();
    double sum = 0;
    for (size_t i = 0; i < countnum; i++) {
        sum += dI2[i] * static_cast<double>(i);
    }
    return 1e8 * sum / static_cast<double>(countnum);
}

double ChiSqDer(const std::valarray<double>& Vnum,
                const std::valarray<double>& Inum,
                const std::valarray<double>& Irex) {
    /*
     *  Calculate Chi² figure of merit for differential resistance
     */

    if (Inum.size() != Vnum.size()) {
        throw std::length_error("Numeric I and V must be of the same size");
    }
    if (Inum.size() != Irex.size()) {
        throw std::length_error("Numeric I and Recalculated I must be of the same size");
    }

    const size_t countnum = Vnum.size();
    const std::valarray dVnum = Vnum[std::slice(1, countnum - 1, 1)]
                                - Vnum[std::slice(0, countnum - 1, 1)];
    const std::valarray dInum = Inum[std::slice(1, countnum - 1, 1)]
                                - Inum[std::slice(0, countnum - 1, 1)];
    const std::valarray dIrex = Irex[std::slice(1, countnum - 1, 1)]
                                - Irex[std::slice(0, countnum - 1, 1)];
    return ChiSq(Inum, Irex)
           + (std::pow((dInum - dIrex) / dVnum, 2) / std::pow(dIrex / dVnum, 2)).sum()
           / static_cast<double>(countnum - 1);
}

/*
 * Not used for fitting now
 *
void GetDBStartPoint(double *par, const char *fname)
{
    CFoo foo(0);
    FILE *DB = fopen(fname, "r");
    int num, numentries = 0, flag = 0, indmin = 0;
    double tmp1, tmp2, tmp3, chimin = 0;
    fscanf(DB, "%d", &num);
    while (flag != EOF) {
        flag = fscanf(DB, "%lf %lf %lf", &tmp1, &tmp2, &tmp3);
        for (int i = 0; i < num; i++)
            flag = fscanf(DB, "%lf %lf", &tmp1, &tmp2);
        if (flag != EOF)
            numentries++;
    }
    freopen(fname, "r", DB);
    fscanf(DB, "%d", &num);
    double *arrI, *arrV, *arrpar1, *arrpar2, *arrpar3;
    arrI = new double[num * numentries];
    arrV = new double[num * numentries];
    arrpar1 = new double[numentries];
    arrpar2 = new double[numentries];
    arrpar3 = new double[numentries];
    for (int i = 0; i < numentries; i++) {
        fscanf(DB, "%lf %lf %lf", &arrpar1[i], &arrpar2[i], &arrpar3[i]);
        for (int j = 0; j < num; j++)
            fscanf(DB, "%lf %lf", &arrV[i * num + j], &arrI[i * num + j]);
    }
    printf("DB loaded successfully: %d entries, %d ppe\n", numentries, num);
    fclose(DB);
    double *tmpI, *tmpV;
    tmpI = new double[num];
    tmpV = new double[num];

    Resample(num, foo.countexp, tmpI, tmpV, foo.Iexp, foo.Vexp, arrI, arrV);
    chimin = ChiSq(num, arrV, arrI, tmpI);
    for (int i = 1; i < numentries; i++) {
        Resample(num, foo.countexp, tmpI, tmpV, foo.Iexp, foo.Vexp, arrI + num * i, arrV + num * i);
        tmp1 = ChiSq(num, arrV + num * i, arrI + num * i, tmpI);
        if (tmp1 < chimin) {
            chimin = tmp1;
            indmin = i;
        }
    }
    par[0] = arrpar1[indmin];
    par[1] = arrpar2[indmin];
    par[2] = arrpar3[indmin];
    printf("Found starting point from DB: %lf %lf %lf, chisq: %e \n", par[0], par[1], par[2], chimin);
    delete[] arrI;
    delete[] arrV;
    delete[] arrpar1;
    delete[] arrpar2;
    delete[] arrpar3;
    delete[] tmpI;
    delete[] tmpV;
    return;
}
*/

// Write two columns of data into a file
void writeIV(const std::string& filename,
             const std::valarray<double>& I,
             const std::valarray<double>& V) {
    /*
     *  Write two columns of data into a file
     */
    if (I.size() != V.size()) {
        throw std::length_error("I and V must be of the same size");
    }

    if (std::ofstream f(filename.c_str(), std::ios::out); f) {
        for (auto pv = begin(V), pi = begin(I); pv != end(V) && pi != end(I); ++pv, ++pi) {
            f << *pv << ' ' << *pi << std::endl;
        }
        f.close();
    }
}

void writeConverg(const double x, const double f, const time_t start) {
    std::clog << '\n' << "CURRENT XMIN = " << x << '\t' << "CHISQMIN = " << f << std::endl;

    if (std::ofstream conv("converg.txt", std::ios::app); conv) {
        conv << x << f << (double) std::difftime(time(nullptr), start) << std::endl;
        conv.close();
    }
}
