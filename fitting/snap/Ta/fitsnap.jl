ii = findlast("MDP", pwd()); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setpath.jl");

using MDP

species = ["Ta"];
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MDP/fitting/snap/Ta/JSON/Displaced_A15";
config = MDP.Preprocessing.readJSONfromFolder(datapath, species);



