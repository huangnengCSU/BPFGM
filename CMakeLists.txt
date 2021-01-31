cmake_minimum_required(VERSION 3.10)
project(bmp)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -maes -mavx -mfma -mavx2 -D_FORTIFY_SOURCE=2 -ffunction-sections -fdata-sections -O3 -Wl,--gc-sections -Wl,--as-needed -Wno-undef")

include_directories(bamtools/src/ src kalign/src)

link_directories(./)

add_executable(block main.cpp src/consensus.cpp)

target_link_libraries(block bamtools z kalign pthread)

add_executable(trivial_training ProvideDatasetTrivial.cpp src/consensus.cpp src/dataset.cpp src/dataset.h)

target_link_libraries(trivial_training bamtools z pthread)

add_executable(complex_training ProvideDatasetComplex.cpp src/consensus.cpp src/dataset.cpp src/dataset.h)

target_link_libraries(complex_training bamtools z pthread)
