cmake_minimum_required(VERSION 3.18)
project(LbmLdcCuda LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Build for the GPU in this machine. Override with e.g. -DCMAKE_CUDA_ARCHITECTURES=86
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(lbm_cuda main.cu inputHandler.cpp)
# cudaUtils.h lives next to main.cu and is included by the matrixMethods2 headers
target_include_directories(lbm_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
