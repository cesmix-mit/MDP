function app = initializeapp(version)

app.codename = "MDP";
app.version = version;
app.appname = "app";
app.platform = "cpu";
app.cpucompiler = "g++";
app.mpicompiler = "mpicxx";
app.gpucompiler = "nvcc";
app.mpirun = "mpirun";
app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
app.gpuflags = "-lcudart -lcublas";
app.potentialfile = "";
app.configfile = "";
app.preprocessmode = 1;
app.mpiprocs = 1;

app.dim = 3;         % physical dimension
app.nconfigs = 1;    % number of configurations
app.ntimesteps = 0;  % number of time steps
app.nab = 4096;      % number of atoms per block
app.backend = 0;     % backend computing platform
app.training = 0;     % 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.runMD = 0;        % 0 no MD simulation, 1 -> run MD simulation
app.potentialform = 0;% 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
app.neightype = 0;    % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.energycal = 1;    % turns energy calculation on or off
app.forcecal = 1;     % turns force calculation on or off
app.stresscal = 0;    % turns stress calculation on or off
app.time = 0.0;       % initial time
app.dt = 0.0;         % time step

% boundary conditions
app.bcs = [0 0 0 0 0 0]';
app.pbc = [1 1 1]'; % periodic boundary conditions

% atom types
app.natomtype = 2;               % number of atom types
app.atomtypes = [1 1 1; 2 2 2]'; % information about each atom type

% molecue types
app.nmoletype = 1;      % number of molecule types
app.moletypes = [1 1]'; % information about each molecule type

% simulation and solver parameters
app.ndims = [1 2]';     % a list of integers indicating dimensions 
app.flag = [0 0]';      % a list of flags 
app.simparam = [0.0 0.0]'; % simulation parameters
app.solparam = [0.0 0.0]'; % solver parameters

% potential parameters
app.mu = [0.0 0.0]';    % regular parameters
app.eta = [0.0 0.0]';   % hyperparameters
app.kappa = [0 0]';     % interger parameters 

% machine learning potentials
app.mu0 = [0 0]';   % coefficients of ML potential
app.K = 3;          % degree of radial basis functions
app.L = 3;          % degree of spherical harmonics     
app.rcut0 = 0.1;    % cut-off radius for machine learning potential
app.descriptor = 0; % descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 2;   % spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.ncmu0 = 2;      % length of mu0

% nonbonded single potentials
app.np1a = 2;           % number of nonbonded single potentials 
app.pot1a = [1 2]';     % list of nonbonded single potentials 
app.mu1a = [0.0 0.0]';  % parameters for nonbonded single potentials 
app.ncmu1a = 2;         % length of mu1a

% bonded single potentials
app.np1b = 2;           % number of bonded single potentials 
app.pot1b = [1 2]';      % list of bonded single potentials 
app.mu1b = [0.0 0.0]';   % parameters for bonded single potentials 
app.ncmu1b = 2;         % length of mu1a
app.atom1b = [1 2]';     % list of bonded atom types

% nonbonded pair potentials
app.np2a = 2;           % number of nonbonded pair potentials 
app.pot2a = [1 2]';     % list of nonbonded pair potentials 
app.mu2a = [0.0 0.0]';  % parameters for all nonbonded pair potentials 
app.ncmu2a = 2;         % length of mu2a
app.rcut2a = [0.1 0.1]'; % cut-off radius for each nonbonded pair potential

% bonded pair potentials
app.np2b = 3;               % number of bonded pair potentials 
app.pot2b = [1 2 3]';       % list of bonded pair potentials 
app.mu2b = [0.0 0.0 0.0]';  % parameters for all bonded pair potentials 
app.ncmu2b = 3;             % length of mu2b
app.rcut2b = [0.1 0.1 0.1]';   % cut-off radius for each nonbonded pair potential
app.atom2b = [1 1; 2 2; 1 2]'; % list of bonded atom pair types

% two-body bond order potentials
app.np2c = 2;               % number of two-body bond order potentials
app.pot2c = [1 2]';         % list of two-body bond order potentials
app.mu2c = [0.0 0.0 0.0]';  % parameters for all two-body bond order potentials
app.ncmu2c = 3;             % length of mu2c
app.rcut2c = [0.1 0.1]';    % cut-off radius for each two-body bond order potential
app.atom2c = [1 2];        % list of bond porder atom types

% nonbonded triplet potentials
app.np3a = 2;           % number of nonbonded triplet potentials
app.pot3a = [1 2]';     % list of nonbonded triplet potentials
app.mu3a = [0.0 0.0]';  % parameters for all nonbonded triplet potentials
app.ncmu3a = 2;         % length of mu3a
app.rcut3a = [0.1 0.1]'; % cut-off radius for each nonbonded triplet potential

% bonded triplet potentials
app.np3b = 2;               % number of bonded triplet potentials 
app.pot3b = [1 2]';         % list of bonded triplet potentials 
app.mu3b = [0.0 0.0 0.0]';  % parameters for all bonded triplet potentials 
app.ncmu3b = 3;             % length of mu3b
app.rcut3b = [0.1 0.1]';    % cut-off radius for each nonbonded triplet potential
app.atom3b = [1 2 1; 2 1 2]'; % list of bonded atom triplet types

% three-body bond order potentials
app.np3c = 2;               % number of three-body bond order potentials
app.pot3c = [1 2]';         % list of three-body bond order potentials
app.mu3c = [0.0 0.0 0.0]';  % parameters for all three-body bond order potentials
app.ncmu3c = 3;             % length of mu3c
app.rcut3c = [0.1 0.1]';    % cut-off radius for each three-body bond order potential
app.atom3c = [1 1; 2 2]';   % list of three-body bond order atom types

% nonbonded quadruplet potentials
app.np4a = 2;           % number of nonbonded quadruplet potentials
app.pot4a = [1 2]';     % list of nonbonded quadruplet potentials
app.mu4a = [0.0 0.0]';  % parameters for all nonbonded quadruplet potentials
app.ncmu4a = 2;          % length of mu4a
app.rcut4a = [0.1 0.1]'; % cut-off radius for each nonbonded quadruplet potential

% bonded quadruplet potentials
app.np4b = 2;               % number of bonded quadruplet potentials 
app.pot4b = [1 2]';         % list of bonded quadruplet potentials 
app.mu4b = [0.0 0.0 0.0]';  % parameters for all bonded quadruplet potentials 
app.ncmu4b = 3;             % length of mu3b
app.rcut4b = [0.1 0.1]';    % cut-off radius for each nonbonded quadruplet potential
app.atom4b = [1 2 1 1; 2 1 2 1]'; % list of bonded atom quadruplet types

app.ncx = 3; % number of compoments of x
app.ncv = 3; % number of compoments of v
app.ncf = 0; % number of compoments of f
app.ncq = 0; % number of compoments of q
app.ncmu = 1; % number of compoments of mu
app.nceta = 1; % number of compoments of eta
app.nckappa = 1; % number of compoments of kappa

end

