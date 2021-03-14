version = "0";
cdir = pwd(); ii = strfind(cdir, "Applications");
run(cdir(1:(ii-1)) + "Installation/setpath.m");
cppsource = cdir(1:(ii-1)) + "C++/Main";

%[app,config] = initializemdp(version);
app = initializeapp(version);

app.dim = 3;          % physical dimension
app.training = 1;     % 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.maxnumneighbors = 100; % maximum number of neighbors allowed
app.runMD = 0;        % 0 no MD simulation, 1 -> run MD simulation
app.potentialform = 0;% 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
app.neighpair = 0;    % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.neighcell = 0;    % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.decomposition = 0;% 0 -> force decomposition, 1 -> atom decomposition
app.descriptor = 0;   % descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 2;     % spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.K = 3;            % degree of radial basis functions
app.L = 3;            % degree of spherical harmonics     
app.rcutml = 0.1;     % cut-off radius for machine learning potential
app.energycal = 1;    % turns energy calculation on or off
app.forcecal = 1;     % turns force calculation on or off
app.bcs = [0 0 0 0 0 0]; % boundary conditions
app.pbc = [1 1 1];       % periodic boundary conditions

app.natomtype = 1;  % number of atom types
app.atomnumbers = [18];
app.atommasses = [39.948];
app.atomcharges = [0];
app.appname = "test";
app.configfile = "configs/CONFIGS1";

% LJ potential
epsilon = 0.01029849;
sigma = 3.4;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
app.pot2a = [1];     % list of nonbonded pair potentials 
app.mu2a = [A B];    % parameters for  LJ potential
app.rcut2a = [8.5];  % cut-off radius 

disp("Read configuration file");
config = readconfig(app, 4);              % read configurations from data file
app.nconfigs = config.nconfigs;           % number of configurations
[app,config] = preprocessing(app,config); 
movefile('testapp.bin',char(cppsource));
movefile('testconfig.bin',char(cppsource));

cd(cppsource);
disp("Compiling C++ source code");
compile("g++");
disp("Running C++ source code");
!./cpuMDP test out
cd(cdir);


