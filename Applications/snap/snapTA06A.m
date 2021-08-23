app = initializeapp();
app.unitstyle = "metal"; % unit system
app.pbc = [1 1 1];       % periodic boundary conditions

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
nlattice = 4;  % note that # of atoms = 2 * nlattice^3 
app.lattice = setlattice("bcc", 3.316);
app.region = setregion([0 0 0], [nlattice nlattice nlattice]);

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

temp = 300;  seed0 = 4928459; distflag = 0; sumflag = 0; loopflag = 2; momflag = 1; rotflag = 0;
app.createvelocity = [temp, seed0, distflag, sumflag, loopflag, momflag, rotflag];

% generate input files and run simulation 
app = mdp(app);


