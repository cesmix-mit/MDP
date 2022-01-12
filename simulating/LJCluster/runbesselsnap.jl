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

# Bessel parameters
mu = [9080.776024081284, 22740.90481568728, 22826.411376542143, 12272.607095393349, 3542.9517448941156, 437.16635703260755]
rcut = 3.294988099200801
alpha = 4.902241343306524

# define Bessel potential
potentials[1] = addpotential("pair", "nonbonded", [0], "scaledbessel", mu, rcut);
include(potentials[1].potentialfunction * ".jl");  # include the potentials file

# SNAP parameters
twojmax = 4
ntypes = 1
rcutfac = rcut;  
rfac0 = 1.0;
rmin0 = 0;
bzeroflag = 0;
quadraticflag = 0;
switchflag = 1;
chemflag = 0;
bnormflag = 0;
wselfallflag = 0;

snaparam = [ntypes, twojmax, rcutfac, rfac0, rmin0, bzeroflag, switchflag, quadraticflag, chemflag, bnormflag, wselfallflag]
elemradius = 0.5*ones(Float64, ntypes)
elemweight = ones(Float64, ntypes)

snapcoeff = [0.0
9.51697531984316e-5
 -3.503209826706362e-6
 -0.00015477632687816527
  8.864246681893688e-5
 -6.228723743014278e-5
  0.0002572563742710612
  1.0535143289316965e-6
  9.529421498586589e-5
  9.123169496032785e-6
 -0.0001693491133185178
  0.00013696751121984198
  3.328433092023331e-5
 -6.500144488870558e-5
 -3.968313413722244e-5]

# define SNAP potential
potentials[2] = initsna(snaparam, elemradius, elemweight, snapcoeff)

app.unitstyle = "lj"; # unit system
app.ensemblemode = "nve";
app.time = 0.0;             # initial time
app.dt = 0.5e-3;            # time step size
app.ntimesteps = 1000;       # time steps
app.globalfreq = 10;        # time frequency to print out global outputs 
app.neighskin = 0.5;        # add a skin distance for neighborlist  
app.eta = [rcut alpha]
app.atommasses = [1.0 1.0];

# run MD simulation using velocity Verlet algorithm 
VelocityVerlet(x, v, q, t, a, b, c, pbc, app, potentials)



