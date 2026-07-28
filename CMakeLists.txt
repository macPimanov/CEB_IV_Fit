cmake_minimum_required(VERSION 3.5)

project(CEB_IV_Fit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(CEB_IV_Fit main.cpp
        constants.h
        CEBNumericModel.h
        utils.h utils.cpp
        IVParamFitter.h IVParamFitter.cpp
        OffsetEliminator.h OffsetEliminator.cpp
        MinimizationAlgorithms.h MinimizationAlgorithms.cpp
        README.md
        SINS1_53_240_303_mK.txt
        startparams.txt
)

include(GNUInstallDirs)
install(TARGETS CEB_IV_Fit
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
