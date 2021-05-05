function app = initializeapp(sourcepath,version)

app.sourcepath = sourcepath;
app.codename = "MDP";
app.version = version;
app.appname = "app";
app.platform = "cpu";
app.cpucompiler = "g++";
app.mpicompiler = "mpicxx";
app.gpucompiler = "";
app.mpirun = "mpirun";
app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
app.gpuflags = "-lcudart -lcublas";
app.cpumacros = "";
app.gpumacros = "";
app.potentialfile = "";
app.configfile = "";
app.weightfile = "";
app.configmode = 0;
app.weightmode = 0;
app.preprocessmode = 1;
app.mpiprocs = 1;
app.traininglist = []; % a list of configurations for training the potential
app.validatelist = []; % a list of configurations for validating the potential

app.dim = 3;         % physical dimension
app.nconfigs = 1;    % number of configurations
app.maxnumneighbors = 10; % maximum number of neighbors allowed
app.ntimesteps = 0;  % number of time steps
app.nab = 4096;      % number of atoms per block
app.backend = 0;     % backend computing platform
app.training = 0;     % 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.dftdata = 0;      % 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
app.runMD = 0;        % 0 no MD simulation, 1 -> run MD simulation
app.potentialform = 0;% 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
app.neighpair = 0;    % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.neighcell = 0;    % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.decomposition = 0;% 0 -> force decomposition, 1 -> atom decomposition
app.energycal = 0;    % turns energy calculation on or off
app.forcecal = 0;     % turns force calculation on or off
app.stresscal = 0;    % turns stress calculation on or off
app.time = 0.0;       % initial time
app.dt = 0.0;         % time step
app.rcutsqmax = 0.0;  % square of maximum cutoff radius
app.boxoffset = [0.0 0.0 0.0]; % offset for simulation box to account for periodic boundary conditions

app.bcs = [0 0 0 0 0 0]; % boundary conditions
app.pbc = [1 1 1];       % periodic boundary conditions

app.we = 1; % energy weight
app.wf = 1; % force weight

% atom types
app.natomtype = 0;  % number of atom types
app.atommasses = [];
app.atomcharges = [];
app.atomnumbers = [];

% molecue types
app.nmoletype = 0;  % number of molecule types
app.moleinfo = []; % information about each molecule type

% simulation and solver parameters
app.ndims = [];     % a list of integers indicating dimensions 
app.flag = [];      % a list of flags 
app.simparam = []; % simulation parameters
app.solparam = []; % solver parameters

% machine learning potentials
app.muml = [];   % coefficients of ML potential
app.K = 0;          % degree of radial basis functions
app.L = 0;          % degree of spherical harmonics     
app.rcutml = 0.0;    % cut-off radius for machine learning potential
app.descriptor = 0; % descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 2;   % spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.chemtype = 0;   % 0 -> single atom-type basis functions, 1 -> pair atom-type basis functions 
app.ncmuml = 0;     % length of mu0

% empirical potential parameters
app.eta = [];   % hyperparameters for emprical potentials
app.kappa = []; % flag parameters for emprical potentials

% nonbonded single potentials
app.npot1a = 0;           % number of nonbonded single potentials 
app.pot1a = [];     % list of nonbonded single potentials 
app.mu1a = [];  % parameters for nonbonded single potentials 
app.ncmu1a = 0;         % length of mu1a

% bonded single potentials
app.npot1b = 0;           % number of bonded single potentials 
app.pot1b = [];      % list of bonded single potentials 
app.mu1b = [];   % parameters for bonded single potentials 
app.ncmu1b = 0;         % length of mu1a
app.atom1b = [];     % list of bonded atom types

% nonbonded pair potentials
app.npot2a = 0;           % number of nonbonded pair potentials 
app.pot2a = [];     % list of nonbonded pair potentials 
app.mu2a = [];  % parameters for all nonbonded pair potentials 
app.ncmu2a = 0;         % length of mu2a
app.rcut2a = []; % cut-off radius for each nonbonded pair potential

% bonded pair potentials
app.npot2b = 0;               % number of bonded pair potentials 
app.pot2b = [];       % list of bonded pair potentials 
app.mu2b = [];  % parameters for all bonded pair potentials 
app.ncmu2b = 0;             % length of mu2b
app.rcut2b = [];   % cut-off radius for each nonbonded pair potential
app.atom2b = []; % list of bonded atom pair types

% two-body bond order potentials
app.npot2c = 0;               % number of two-body bond order potentials
app.pot2c = [];         % list of two-body bond order potentials
app.mu2c = [];  % parameters for all two-body bond order potentials
app.ncmu2c = 0;             % length of mu2c
app.rcut2c = [];    % cut-off radius for each two-body bond order potential
app.atom2c = [];        % list of bond porder atom types

% nonbonded triplet potentials
app.npot3a = 0;           % number of nonbonded triplet potentials
app.pot3a = [];     % list of nonbonded triplet potentials
app.mu3a = [];  % parameters for all nonbonded triplet potentials
app.ncmu3a = 0;         % length of mu3a
app.rcut3a = []; % cut-off radius for each nonbonded triplet potential

% bonded triplet potentials
app.npot3b = 0;               % number of bonded triplet potentials 
app.pot3b = [];         % list of bonded triplet potentials 
app.mu3b = [];  % parameters for all bonded triplet potentials 
app.ncmu3b = 0;             % length of mu3b
app.rcut3b = [];    % cut-off radius for each nonbonded triplet potential
app.atom3b = []; % list of bonded atom triplet types

% three-body bond order potentials
app.npot3c = 0;               % number of three-body bond order potentials
app.pot3c = [];         % list of three-body bond order potentials
app.mu3c = [];  % parameters for all three-body bond order potentials
app.ncmu3c = 0;             % length of mu3c
app.rcut3c = []';    % cut-off radius for each three-body bond order potential
app.atom3c = []';   % list of three-body bond order atom types

% nonbonded quadruplet potentials
app.npot4a = 0;           % number of nonbonded quadruplet potentials
app.pot4a = [];     % list of nonbonded quadruplet potentials
app.mu4a = [];  % parameters for all nonbonded quadruplet potentials
app.ncmu4a = 0;          % length of mu4a
app.rcut4a = []; % cut-off radius for each nonbonded quadruplet potential

% bonded quadruplet potentials
app.npot4b = 0;               % number of bonded quadruplet potentials 
app.pot4b = [];         % list of bonded quadruplet potentials 
app.mu4b = [];  % parameters for all bonded quadruplet potentials 
app.ncmu4b = 0;             % length of mu3b
app.rcut4b = [];    % cut-off radius for each nonbonded quadruplet potential
app.atom4b = []; % list of bonded atom quadruplet types

app.ncx = 0; % number of compoments of x
app.ncv = 0; % number of compoments of v
app.nce = 0; % number of compoments of e
app.ncf = 0; % number of compoments of f
app.ncq = 0; % number of compoments of q
app.ncmu = 0; % number of compoments of mu
app.nceta = 0; % number of compoments of eta
app.nckappa = 0; % number of compoments of kappa

end

