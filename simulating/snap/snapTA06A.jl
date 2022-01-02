# Find and add MDP's path to Julia search path
cdir = pwd(); ii = findlast("MDP", cdir);
if ii===nothing
    # MDPpath = /path/to/MDP;
    display("MDP's path is not found. Please uncomment the above line and set the path to MDP.");
else
    MDPpath = cdir[1:ii[end]] * "/";    
    include(MDPpath * "Installation/setpath.jl");
end

# External packages
using Revise, DelimitedFiles, SymPy

# MDP packages
using Preprocessing, Gencode, Postprocessing

app = Preprocessing.initializeapp(cdir,MDPpath);
app.platform = "cpu";    # app.platform = "gpu" => MDP runs on Nvidia GPU 
app.unitstyle = "metal"; # unit system
app.pbc = [1 1 1];       # periodic boundary conditions
nlattice = 4;            # note that # of atoms = 2 * nlattice^3 
app.appname = "snapTA06A" * string(nlattice);

#  neighbor list
app.neighskin = 1.0;
app.neighcell = 0;       # 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.neighmax = 100;      # maximum number of neighbors allowed
app.neighevery = 1;      # perform neighbor rebuild list for "every" iterations
app.neighdelay = 0;      # delay neighbor rebuild
app.neighcheck = 1;      # 0 -> neighbor rebuild by delay and every settings, 
                         # 1 -> neighbor rebuild occurs if some atom has moved more than half the skin distance
app.neighpair = 0;       # 0 -> full neighbor list, 1 -> half neighbor list for pair potentials

# atom 
app.atommasses = reshape([180.88],1,1);
app.atomcharges = reshape([0.0],1,1);

# lattice, region, domain
app.lattice = Preprocessing.setlattice("bcc", 3.316);
app.region = Preprocessing.setregion([nlattice nlattice nlattice]);

# nonparametric potential descriptors
app.descriptor = "snap";    # (snap, shp)
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
app.snapnelem = 1;    # # of elements
app.snapncoeff = 31;  # # of coeffients per element
app.snapelemradius = reshape([0.5],1,1); # radius per element
app.snapelemweight = reshape([1.0],1,1);   # weight per element
app.snapcoefffile = "snapTA06A.coeff";
app.rcutml = app.snaprcutfac;       # cut-off radius for nonparametric potential

# parametric potential descriptors
app.potentialfile = "ZBLpotential.jl";   # empirical potential file
include(app.potentialfile);  # include the potential file
zblcutinner = 4;
zblcutouter = 4.8;
zblz = 73;
app.pot2a = reshape([1],1,1);     # list of nonbonded pair potentials 
app.mu2a = [zblcutinner zblcutouter zblz];    # parameters for ZBL potential
app.rcut2a = reshape([zblcutouter],1,1);  # cut-off radius 

# time integration and outputs 
app.ensemblemode = "nve";
app.time = 0.0;             # initial time
app.dt = 0.5e-3;            # time step size
app.ntimesteps = 100;       # # time steps
app.globalfreq = 10;

temp = 300;  seed0 = 4928459; distflag = 0; sumflag = 0; loopflag = 2; momflag = 1; rotflag = 0;
app.createvelocity = [temp seed0 distflag sumflag loopflag momflag rotflag];

# generate input files and run simulation 
app, config = mdp(app);



