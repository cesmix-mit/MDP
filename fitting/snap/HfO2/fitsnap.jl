cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setpath.jl");

using MDP

filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MDP/fitting/snap/HfO2/EXYZ/HfO2_cpmd_1.xyz";
species = ["Hf", "O"];
config = MDP.Preprocessing.readEXTXYZ(filename, species);




