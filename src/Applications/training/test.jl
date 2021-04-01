version = "0";

# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Applications", cdir);
sourcepath = cdir[1:(ii[1]-1)];
include(sourcepath * "/Installation/setpath.jl");

# MDP packages
using Preprocessing, Gencode, Postprocessing

app = Preprocessing.initializeapp(sourcepath,version);
app.appname = "test";
app.currentdir = cdir; # current directory
app.cpucompiler = "/usr/local/Cellar/llvm/11.1.0/bin/clang++";
app.cpuflags = "-lblas -llapack -Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib";
app.gpucompiler = "";  # GPU is not supported yet
app.configfile = "configs/CONFIGS1"; # configuration file
app.configmode = 4;   # mode for reading configuration file
app.traininglist = reshape(0:79,1,80);  # a list of configurations for training the potential
app.validatelist = reshape(60:100,1,41);# a list of configurations for validating the potential
app.dim = 3;          # physical dimension
app.maxnumneighbors = 100; # maximum number of neighbors allowed
app.runMD = 0;        # 0 no MD simulation, 1 -> run MD simulation
app.training = 1;     # 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.dftdata = 3;      # 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
app.potentialform = 1;# 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
app.neighpair = 0;    # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.neighcell = 0;    # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.decomposition = 0;# 0 -> force decomposition, 1 -> atom decomposition
app.energycal = 1;    # turns energy calculation on or off
app.forcecal = 1;     # turns force calculation on or off
app.pbc = [1 1 1];    # periodic boundary conditions

# atom types
app.natomtype = 1;  # number of atom types
app.atomnumbers = reshape([18],1,1);
app.atommasses = reshape([39.948],1,1);
app.atomcharges = reshape([0],1,1);

# Machine learning descriptors
app.descriptor = 0;   # descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 0;     # spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.chemtype = 0;     # 0 -> single atom-type basis functions, 1 -> pair atom-type basis functions  
app.K = 6;            # number roots of radial basis functions
app.L = 4;            # degree of radial basis functions and spherical harmonics     
app.rcutml = 8.5;     # cut-off radius for machine learning potential

# Empirical potential descriptors
app.potentialfile = "LJpotential.jl";   # empirical potential file
include(app.potentialfile);  # include the potential file
epsilon = 0.01029849;
sigma = 3.4;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
app.pot2a = reshape([1],1,1);     # list of nonbonded pair potentials 
app.mu2a = [A B];    # parameters for  LJ potential
app.rcut2a = reshape([8.5],1,1);  # cut-off radius 

# train and validate the potential
app, config = Postprocessing.mdp(app);




# app, config = Preprocessing.preprocessing(app); # preprocess input data
# Gencode.gencode(app); # Generate C++ code
# Gencode.compile(app);
# Gencode.runcode(app);











# app.appname = "jl";
# app.configfile = "configs/CONFIGS1"; # configuration file
# app.configmode = 4;   # mode for reading configuration file
# app.dim = 3;          # physical dimension
# app.training = 1;     # 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
# app.maxnumneighbors = 100; # maximum number of neighbors allowed
# app.runMD = 0;        # 0 no MD simulation, 1 -> run MD simulation
# app.potentialform = 0;# 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
# app.neighpair = 0;    # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
# app.neighcell = 0;    # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
# app.decomposition = 0;# 0 -> force decomposition, 1 -> atom decomposition
# app.descriptor = 0;   # descriptor flag: 0 -> Spherical Harmonics Bessel
# app.spectrum = 0;     # spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
# app.chemtype = 0;     # 0 -> single atom-type basis functions, 1 -> double atom-type basis functions  
# app.K = 2;            # number roots of radial basis functions
# app.L = 6;            # degree of radial basis functions and spherical harmonics     
# app.rcutml = 8.5;     # cut-off radius for machine learning potential
# app.energycal = 1;    # turns energy calculation on or off
# app.forcecal = 1;     # turns force calculation on or off
# app.bcs = [0 0 0 0 0 0]; # boundary conditions
# app.pbc = [1 1 1];       # periodic boundary conditions

# # LJ potential
# app.potentialfile = "LJpotential";
# app.natomtype = 1;  # number of atom types
# app.atomnumbers = reshape([18],1,1);
# app.atommasses = reshape([39.948],1,1);
# app.atomcharges = reshape([0],1,1);

# epsilon = 0.01029849;
# sigma = 3.4;
# A = 4*epsilon*sigma^12;
# B = 4*epsilon*sigma^6;
# app.pot2a = reshape([1],1,1);     # list of nonbonded pair potentials 
# app.mu2a = [A B];    # parameters for  LJ potential
# app.rcut2a = reshape([8.5],1,1);  # cut-off radius 

# app,config = preprocessing(app); # preprocess input data
# # gencode(app);                      # generate code

# # cd(cppsource);
# # disp("Compiling C++ source code");
# # compile("/usr/local/Cellar/llvm/11.1.0/bin/clang++");
# # disp("Running C++ source code");
# # #!./cpuMDP test out
# # cd(cdir);


