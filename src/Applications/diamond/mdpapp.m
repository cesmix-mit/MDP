version = "0";

cdir = pwd(); 
ii = strfind(cdir, "Applications");
sourcepath = cdir(1:(ii-1));
run(sourcepath + "Installation/setpath.m");

app = initializeapp(sourcepath,version);
app.cpucompiler = "/usr/local/Cellar/llvm/11.1.0/bin/clang++";
app.cpuflags = "-Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib -lblas -llapack";
%app.cpucompiler = "g++";
%app.cpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
%app.cpuflags = "-lblas -llapack";
app.gpucompiler = [];  % GPU is not supported yet
app.appname = "test";
app.currentdir = cdir; % current directory
app.configfile = "diamond.configs"; % configuration file
app.configmode = 4;   % LAMMPS mode for reading configuration file
app.traininglist = 0:101;  % a list of configurations for training the potential
app.validatelist = 0:101;% a list of configurations for validating the potential
app.dim = 3;          % physical dimension
app.maxnumneighbors = 160; % maximum number of neighbors allowed
app.runMD = 0;        % 0 no MD simulation, 1 -> run MD simulation
app.training = 1;     % 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.dftdata = 3;      % 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
app.potentialform = 1;% 0 -> empirical potential, 1 -> ML potential, 2 -> hybrid empirical+ML potential    
app.neighpair = 0;    % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.neighcell = 0;    % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.decomposition = 0;% 0 -> force decomposition, 1 -> atom decomposition
app.energycal = 1;    % turns energy calculation on or off
app.forcecal = 1;     % turns force calculation on or off
app.pbc = [1 1 1];    % periodic boundary conditions

% atom types
app.natomtype = 1;  % number of atom types
app.atomnumbers = [18];
app.atommasses = [39.948];
app.atomcharges = [0];

% Machine learning potential descriptors
app.descriptor = 0;   % descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 0;     % spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.chemtype = 0;     % 0 -> single atom-type basis functions, 1 -> pair atom-type basis functions  
app.K = 4;            % number roots of radial basis functions
app.L = 3;            % degree of radial basis functions and spherical harmonics     
app.rcutml = 5.0;     % cut-off radius for machine learning potential

% Empirical potential descriptors
app.potentialfile = "LJpotential";   % empirical potential file
epsilon = 1.0; %0.01029849;
sigma = 2.0;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
app.pot2a = [1];     % list of nonbonded pair potentials 
app.mu2a = [A B];    % parameters for  LJ potential
app.rcut2a = [5.0];  % cut-off radius 

% train and validate the potential
app = mdp(app);

% get validation data
valid = validate(app);




