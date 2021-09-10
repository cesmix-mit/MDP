################################################################################                   
#                    Molecular Dynamics Potentials (MDP)
#                           CESMIX-MIT Project  
# 
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
################################################################################
# mkdir exec
# cd exec
# cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation 
# cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON  ../Installation  
# cmake --build .
################################################################################

cmake_minimum_required(VERSION 3.12)

project(mdp LANGUAGES CXX)

get_filename_component(MDP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

set(MDP_SOURCE_DIR     ${MDP_DIR}/C++)
set(MDP_MAIN_DIR       ${MDP_DIR}/C++/Main)
set(MDP_CORE_DIR       ${MDP_DIR}/C++/Core)
set(MDP_POTENTIALS_DIR ${MDP_DIR}/C++/Potentials)

file(GLOB MAIN_SOURCES ${MDP_MAIN_DIR}/main.cpp)
file(GLOB CPUCORE_SOURCES ${MDP_CORE_DIR}/cpuCore.cpp)
file(GLOB CPUPOTENTIAL_SOURCES ${MDP_POTENTIALS_DIR}/cpuEmpiricalPotentials.cpp)
file(GLOB GPUCORE_SOURCES ${MDP_CORE_DIR}/gpuCore.cu)
file(GLOB GPUPOTENTIAL_SOURCES ${MDP_POTENTIALS_DIR}/gpuEmpiricalPotentials.cu)    

if(MDP_POTENTIALS)
    add_library(cpuEmpiricalPotentials SHARED ${CPUPOTENTIAL_SOURCES})
    target_compile_options(cpuEmpiricalPotentials PRIVATE -ffast-math -O3 -DNDEBUG) 

    if(MDP_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)        

        cuda_add_library(gpuEmpiricalPotentials SHARED ${GPUPOTENTIAL_SOURCES} OPTIONS -arch sm_70 -std c++11 -O3 -w -D_FORCE_INLINES -DNDEBUG --use_fast_math)    
        target_compile_options(gpuEmpiricalPotentials PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++11 --use_fast_math -arch=sm_70 -D_FORCE_INLINES -O3 -w -DNDEBUG>) 
    endif()
endif()

if(MDP_CORES)
    add_library(cpuCore STATIC ${CPUCORE_SOURCES})
    target_compile_options(cpuCore PRIVATE -ffast-math -O3 -DNDEBUG) 

    if(MDP_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)        

        cuda_add_library(gpuCore STATIC ${GPUCORE_SOURCES} OPTIONS -arch sm_70 -std c++11 -O3 -w -D_FORCE_INLINES -DNDEBUG --use_fast_math)
        target_compile_options(gpuCore PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++11 --use_fast_math -arch=sm_70 -D_FORCE_INLINES -O3 -w -DNDEBUG>) 
    endif()
endif()

if(MDP_EXECUTABLES)

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND AND BLAS_FOUND)
       set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()

    if(MDP_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)

        add_executable(gpuMDP ${MAIN_SOURCES})
        target_compile_definitions(gpuMDP PRIVATE -DUSE_CUDA)    
        target_link_libraries(gpuMDP PRIVATE cpuCore cpuEmpiricalPotentials gpuCore gpuEmpiricalPotentials ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
        target_compile_options(gpuMDP PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)     
    else()
        add_executable(cpuMDP ${MAIN_SOURCES})
        target_link_libraries(cpuMDP PRIVATE cpuCore cpuEmpiricalPotentials ${lapackblas_libraries})
        target_compile_options(cpuMDP PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG) 
    endif()
endif()

