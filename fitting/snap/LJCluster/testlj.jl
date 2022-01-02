cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setpath.jl");

using Preprocessing, Potential 

# Read xyz file to obtain configurations 
filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster/curated_lj_cluster.xyz";
species = ["Ar"];
config = Preprocessing.readEXTXYZ(filename, species);

# shift atom positions to put them in the box [0,10]x[0,10]x[0,10]
config.x = config.x .+ 5.0

# cumalative sum of numbers of atoms 
natom = [0; cumsum(config.natom[:])];

ci = 1; # configuration index in the data 
x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
t = config.t[:,(natom[ci]+1):natom[ci+1]]                        # atom types of configuration ci
q = [];                                   # atom charges of configuration ci

# 3 vectors define the simulation box
a = [10.0; 0.0; 0.0]
b = [0.0; 10.0; 0.0]
c = [0.0; 0.0; 10.0]

# periodic boundary conditions
pbc = [0; 0; 0]

# LJ parameters
epsilon = 1.0; 
sigma = 1.0;
A = 2*epsilon*sigma^12;
B = 2*epsilon*sigma^6;
rcut = 4.0;

# define a list of potentials
potential = Array{Any}(undef, 1)
potential[1] = Potential.addpotential("pair", "nonbonded", [0], "ljpotential", [A; B], rcut);
for i = 1:length(potential)
    include(potential[i].potentialfunction * ".jl");  # include the potential file
end

# compute total energy, peratom energies, and forces
rcutmax = rcut; eta = [0]; kappa = [0];
etot, e, f, v, vtot = Potential.empiricalpotential(config.x[:,1:13], q, t, a, b, c, pbc, rcutmax, eta, kappa, potential);

# using LAMMPS.jl to compute energy and forces
unitstyle = "lj"
pair_style = "lj/cut $rcut"  
pair_coeff = ["* * $epsilon $sigma"]
pe, forces, id, pos, type = Potential.lammpspotential(x, t, a, b, c, pbc, unitstyle, pair_style, pair_coeff)

# compare with LAMMPS.jl
print("Error in potential energy: "); display(abs(etot - pe));
print("Error in force: "); display(maximum(abs.(f[:,id] .- forces)));
print("Error in stress: "); display(maximum(abs.(vtot .- config.stress[:,ci])));

