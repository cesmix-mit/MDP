# MDP.jl


# Dependency
MDP requires Clang++/LLVM compiler, Enzyme, CUDA Toolkit, and Blas/Lapack libaries. 


# Installation
 
On a Linux machine, open the terminal and go to the directory MDP.jl/src/enzyme

  sudo apt update

  sudo apt install clang

  cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -DLLVM_DIR=/path/to/llvm/lib/cmake/clang
  
  make
 



