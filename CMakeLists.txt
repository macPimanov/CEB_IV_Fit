cmake_minimum_required(VERSION 3.5)

project(CEB_IV_Fit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Threads REQUIRED)

add_executable(CEB_IV_Fit main.cpp
        constants.h
        CEBNumericModel.h
        utils.h utils.cpp
        IVParamFitter.h IVParamFitter.cpp
        OffsetEliminator.h OffsetEliminator.cpp
        MinimizationAlgorithms.h MinimizationAlgorithms.cpp
        ceb_library.h ceb_library.cpp
)

target_link_libraries(CEB_IV_Fit Threads::Threads)

# Create shared library for Python bindings
add_library(ceb_library SHARED
        ceb_library.cpp
        ceb_library.h
        constants.h
        CEBNumericModel.h
        MinimizationAlgorithms.cpp
        MinimizationAlgorithms.h
        OffsetEliminator.cpp
        OffsetEliminator.h
        utils.cpp
        utils.h
)

target_link_libraries(ceb_library Threads::Threads)

# Set output directory for the shared library
set_target_properties(ceb_library PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python"
        POSITION_INDEPENDENT_CODE ON
)

include(GNUInstallDirs)
install(TARGETS CEB_IV_Fit
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
