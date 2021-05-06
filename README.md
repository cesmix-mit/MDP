# MDP.jl


# Dependency
MDP requires Clang/LLVM compiler, Enzyme, Blas/Lapack libaries, and CUDA Toolkit (optional). As MDP is optimized for Nvidia GPUs, it is highly recommended to use MDP on GPUs. 

# Installation
 
On a Linux machine, open the terminal and go to the directory MDP.jl/src/enzyme, then install Clang/LLVM comipler and Enzyme:

    sudo apt update

    sudo apt install clang

    cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -DLLVM_DIR=/path/to/llvm/lib/cmake/llvm
  
    make
 
You can skip installing Clang/LLVM compiler, if it is already installed on your machine. Note that you need to know where LLVM is installed on your computer in order to link LLVM libarries located in /path/to/llvm/lib/cmake/llvm to Enzyme.  
 
On a MacOS computer, both LLVM and Enzyme can be installed using 

    brew install enzyme

Installing Enzyme will create two dynamic libaries LLVMEnzyme and ClangEnzyme. MDP will need the path to the Enzyme library ClangEnzyme in order to compile the code.  

