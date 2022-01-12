cdir = pwd(); ii = strfind(cdir, "MDP");
if isempty(ii) == 0
    MDPpath = cdir(1:(ii+2)) + "/";
    run(MDPpath + "Installation/setpath.m");
else
    % MDPpath = /path/to/MDP;
    disp("MDP's path is not found. Please uncomment the above line and set the path to MDP.");
end

app = initializeapp(cdir,MDPpath);
app.buildexec = 1;
app.platform = "cpu";    % app.platform = "gpu" => MDP runs on Nvidia GPU 
app.unitstyle = "metal"; % unit system
app.pbc = [1 1 1];       % periodic boundary conditions
nlattice = 4;            % note that # of atoms = 2 * nlattice^3 
app.appname = "snapWBe" + num2str(nlattice);

%  neighbor list
app.neighskin = 1.0;
app.neighcell = 0;       % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.neighmax = 100;      % maximum number of neighbors allowed
app.neighevery = 1;      % perform neighbor rebuild list for "every" iterations
app.neighdelay = 0;      % delay neighbor rebuild
app.neighcheck = 1;      % 0 -> neighbor rebuild by delay and every settings, 
                         % 1 -> neighbor rebuild occurs if some atom has moved more than half the skin distance
app.neighpair = 0;       % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials

% atom 
app.atommasses = [183.84 9.012182];
app.atomcharges = [0 0];

% lattice, region, domain
latticetype = [1 1]; % 1 => W, 2 => Be
app.lattice = setlattice("bcc", 3.1803,[],[],[],[],latticetype);
app.region = setregion([nlattice nlattice nlattice]);
app.setatomtypefraction = [2 0.05 3590153]; % change 5% of W to Be

% nonparametric potential descriptors
app.descriptor = "snap";    % (snap, shp)
app.snaprcutfac = 4.8123;  
app.snaptwojmax = 8;
app.snaprfac0 = 0.99363;
app.snaprmin0 = 0;
app.snapbzeroflag = 1;
app.snapquadraticflag = 0;
app.snapswitchflag = 1;
app.snapchemflag = 0;
app.snapbnormflag = 0;
app.snapwselfallflag = 0;
app.snapnelem = 2;    % # of elements
app.snapncoeff = 56;  % # of coeffients per element
app.snapelemradius = [0.5 0.417932]; % radius per element
app.snapelemweight = [1 0.959049];   % weight per element
app.snapcoefffile = "snapWBe.coeff";
app.rcutml = app.snaprcutfac;       % cut-off radius for nonparametric potential

% parametric potential descriptors
app.potentialfile = "ZBLpotentialWBe";   % empirical potential file
zblcutinner = 4;
zblcutouter = 4.8;
zblz1 = 74;
zblz2 = 4;
app.pot2a = [1];     % list of nonbonded pair potentials 
app.mu2a = [zblcutinner zblcutouter zblz1 zblz2];    % parameters for ZBL potential
app.rcut2a = [zblcutouter];  % cut-off radius 

% time integration and outputs 
app.ensemblemode = "nve";
app.time = 0.0;             % initial time
app.dt = 0.5e-3;            % time step size
app.ntimesteps = 10000;       % # time steps
app.globalfreq = 10;

temp = 300;  seed0 = 4928459; distflag = 0; sumflag = 0; loopflag = 2; momflag = 1; rotflag = 0;
app.createvelocity = [temp, seed0, distflag, sumflag, loopflag, momflag, rotflag];

% generate input files and run simulation 
app = mdp(app);


