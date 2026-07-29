#include <cmath>
#include <chrono>
#include <iostream>

#include "constants.h"
#include "IVParamFitter.h"

int main() {
    const std::chrono::time_point<std::chrono::steady_clock> beginning = std::chrono::steady_clock::now();
    try {
        IVParamFitter foo;
        foo.loadExperimentData(DATA_FILE_NAME, false);
        foo.computeCEBProperties();

        auto [Irex, Vrex] = foo.resample();

        foo.SeqFit(3, Irex);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    std::clog << '\n' << "Finished." << std::endl;
    std::clog << "Spent " << std::chrono::duration<double>(std::chrono::steady_clock::now() - beginning) << std::endl;

    // std::getchar();
    return 0;
}
