cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setup.jl");

# Read xyz file to obtain configurations 
filename = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster/curated_lj_cluster.xyz";
species = ["Ar"];
config = readEXTXYZ(filename, species);

# shift atom positions to put them in the box [0,10]x[0,10]x[0,10]
config.x = config.x .+ 5.0

# cumalative sum of numbers of atoms 
natom = [0; cumsum(config.natom[:])];

ci = 1; # choose a configuration to start MD simulation
x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci
q = [];                                   # atom charges of configuration ci
v = 0.0*x # zero velocity 

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
rcut = 3.2934420858948856; 

# define LJ potential
potentials[1] = addpotential("pair", "nonbonded", [0], "ljpotential", [A; B], rcut);
include(potentials[1].potentialfunction * ".jl");  # include the potentials file

app.unitstyle = "lj"; # unit system
app.ensemblemode = "nve";
app.time = 0.0;             # initial time
app.dt = 0.5e-3;            # time step size
app.ntimesteps = 1000;       # time steps
app.globalfreq = 10;        # time frequency to print out global outputs 
app.neighskin = 0.5;        # add a skin distance for neighborlist  
app.atommasses = [1.0 1.0];

# run MD simulation using velocity Verlet algorithm 
VelocityVerlet(x, v, q, t, a, b, c, pbc, app, potentials)



