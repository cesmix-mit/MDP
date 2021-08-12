cdir = pwd(); ii = strfind(cdir, "Applications");
sourcepath = cdir(1:(ii-1));
run(sourcepath + "Installation/setpath.m");

app = initializeapp(sourcepath,version);
app.appname = "snapTA06A";
app.currentdir = cdir; % current directory
app.cpucompiler = "/usr/local/Cellar/llvm/11.1.0/bin/clang++";
app.cpumacros = "-std=c++11 -D _ENZYME";
app.cpuflags = "-Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib";
% app.cpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
% app.cpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
% app.cpumacros = "-std=c++11 -D _ENZYME";
% app.gpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
% app.gpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
% app.gpumacros = "--cuda-gpu-arch=sm_60 -std=c++11 -D _ENZYME";
app.buildcorelib = 1;

% unit, boundary condition, and force calculation
app.unitstyle = "metal";
app.pbc = [1 1 1];    % periodic boundary conditions
app.decomposition = 0;% 0 -> force decomposition, 1 -> atom decomposition
app.potentialform = 0;% 0 -> empirical potential, 1 -> ML potential, 2 -> hybrid empirical+ML potential    

%  neighbor list
app.neighskin = 1.0;
app.neighcell = 0;       % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.neighmax = 200;      % maximum number of neighbors allowed
app.neighevery = 1;      % perform neighbor rebuild list for "every" iterations
app.neighdelay = 0;      % delay neighbor rebuild
app.neighcheck = 1;      % 0 -> neighbor rebuild by delay and every settings, 
                         % 1 -> neighbor rebuild occurs if some atom has moved more than half the skin distance
app.neighpair = 0;       % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials

% atom 
app.atommasses = [180.88];
app.atomcharges = [0];

% lattice, region, domain
app.lattice = setlattice("bcc", 3.316);
app.region = setregion([0 0 0], [4 4 4]);

% nonparametric potential descriptors
app.descriptor = "snap";    % (snap, shp)
app.snaprcutfac = 4.67637;  
app.snaptwojmax = 6;
app.snaprfac0 = 0.99363;
app.snaprmin0 = 0;
app.snapbzeroflag = 0;
app.snapquadraticflag = 0;
app.snapswitchflag = 1;
app.snapchemflag = 0;
app.snapbnormflag = 0;
app.snapwselfallflag = 0;
app.snapnelem = 1;    % # of elements
app.snapncoeff = 31;  % # of coeffients per element
app.snapelemradius = 0.5; % radius per element
app.snapelemweight = 1;   % weight per element
app.snapcoefffile = "snapTA06A.coeff";
app.rcutml = app.snaprcutfac;       % cut-off radius for nonparametric potential

% parametric potential descriptors
app.potentialfile = "ZBLpotential";   % empirical potential file
zblcutinner = 4;
zblcutouter = 4.8;
zblz = 73;
app.pot2a = [1];     % list of nonbonded pair potentials 
app.mu2a = [zblcutinner zblcutouter zblz];    % parameters for ZBL potential
app.rcut2a = [zblcutouter];  % cut-off radius 

% time integration and outputs 
app.ensemblemode = "nve";
app.time = 0.0;             % initial time
app.dt = 0.5e-3;            % time step size
app.ntimesteps = 100;       % # time steps
app.globalfreq = 10;        % frequency to print and save default global outputs (pe, ke, ce, te, temp, press)
app.peratomfreq = 20;       % frequency to save default peratom outputs (id, t, x, f)
% app.globaloutputs = ["stresses"];    % additional global outputs (stresses, com, vcm)
% app.peratomoutputs = ["velocity", "mass", "image"];   % additional peratom outputs (velocity, mass, virial, image)

temp = 300;  seed0 = 4928459; distflag = 0; sumflag = 0; loopflag = 2; momflag = 1; rotflag = 0;
app.createvelocity = [temp, seed0, distflag, sumflag, loopflag, momflag, rotflag];

% generate input files and run simulation 
app = mdp(app);

