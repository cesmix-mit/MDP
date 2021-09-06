# Molecular Dynamics Potentials (MDP)

# Installation

mkdir exec 

cd exec 

cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON ../Installation 

cmake --build .

For CUDA support, please use 

cmake -D MDP_POTENTIALS=ON -D MDP_CORES=ON -D MDP_EXECUTABLES=ON -D MDP_CUDA=ON  ../Installation  

# Applications

Run snapTA06A script in the Applications/snap folder 
