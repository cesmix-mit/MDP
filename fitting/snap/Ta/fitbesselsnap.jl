cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setup.jl");

# periodic boundary conditions
pbc = [1; 1; 1]

# path to the database 
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/Research/Software/MDP/fitting/snap/Ta/JSON/";
dataformat = "json"
fileextension = "json"
atomspecies = ["Ta"];

# folders in the datapath 
folders = ["Displaced_A15", "Displaced_BCC", "Displaced_FCC", "Elastic_BCC",
            "Elastic_FCC", "GSF_110", "GSF_112", "Liquid", "Surface",
            "Volume_A15", "Volume_BCC", "Volume_FCC"];

# weights for energies, forces, and stresses in the linear fit
weightinner=[100            1               1.00E-08
            100             1               1.00E-08
            100             1               1.00E-08
            1.00E-08        1.00E-08        0.0001
            1.00E-09        1.00E-09        1.00E-09
            100             1               1.00E-08
            100             1               1.00E-08
            4.67E+02        1               1.00E-08
            100             1               1.00E-08
            1.00E+00        1.00E-09        1.00E-09
            1.00E+00        1.00E-09        1.00E-09
            1.00E+00        1.00E-09        1.00E-09]

# randomly selecting the configurations in the database
randomize = false;

# use all the data 
percentage = 100.0;

# translate atom positions 
translationvector = nothing

# rotate atom positions 
rotationmatrix = nothing

# transpose lattice vectors because the lattice ordering is python-based
transposelattice = true 

# training data 
for i = 1:length(folders)
    traindata[i] = adddata(datapath * folders[i], dataformat, fileextension, 
                percentage, randomize, atomspecies, weightinner[i,:], translationvector, 
                rotationmatrix, transposelattice)
end

rcut = 4.67637   # cut-off radiuss 

# single descriptors to catpure energy of isolated pure elements
descriptors[1] = PotentialStruct("single", "nonbonded", [1], "single", [0.0], rcut);
include(descriptors[1].potentialfunction * ".jl");  # include the descriptors file
      
besselscale = 2.5;
eta = [rcut; besselscale]; 
kappa = [6];      # number of Bessel descriptors

# Bessel descriptors
descriptors[2] = PotentialStruct("pair", "nonbonded", [1, 1], "scaledbessel", [0.0], rcut);
include(descriptors[2].potentialfunction * ".jl");  # include the descriptors file

# SNAP parameters
twojmax = 6
ntypes = 1
rcutfac = rcut;  
rfac0 = 0.99363;
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

# SNAP descriptors
descriptors[3] = initsna(snaparam, elemradius, elemweight)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforcestress"
# use least-square method
method = "lsq" 
# divide energies by the number of atoms
normalizeenergy = true
# divide virials by the box volume
normalizestress = true
# optimization parameters
optim = setoptim(lossfunc, method, normalizeenergy, normalizestress, eta, kappa)

# linear fit to compute SNAP coefficients
coeff,~ = linearfit2(traindata, descriptors, potentials, optim, pbc)

print("SNAP Coeffficients: "), show(stdout, "text/plain", coeff)

# compute unweighted MAE, RMSE, RSQ errors 
energyerrors, forceerrors, stresserrors = validate2(traindata, descriptors, potentials, optim, coeff, pbc)    

printerrors(folders, energyerrors, "Energy Errors")
printerrors(folders, forceerrors, "Force Errors")
printerrors(folders, stresserrors, "Stress Errors")

# ----------------------------------------------------------------------------------------
# Energy Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.07752835884217    |    0.22705475009736    |    0.99596448956663   |
# Displaced_A15 |   0.00049192909655    |    0.00058897195332    |    0.95326317157895   |
# Displaced_BCC |   0.00179938127019    |    0.00240364187430    |    0.99758372333633   |
# Displaced_FCC |   0.00229495671435    |    0.00240573323954    |    0.63958365406662   |
# Elastic_BCC   |   0.01227624775096    |    0.01229967922620    |    -735.09749881215   |
# Elastic_FCC   |   0.00281118185418    |    0.00369493891024    |    -2.3502044635248   |
# GSF_110       |   0.00473898165805    |    0.00555933602248    |    0.77894211169715   |
# GSF_112       |   0.00556624505684    |    0.00617549243995    |    0.89842350462463   |
# Liquid        |   0.00160658042785    |    0.00180561530815    |    0.99960899813060   |
# Surface       |   0.01298185518923    |    0.01588038702008    |    0.93534466141737   |
# Volume_A15    |   0.20707055289531    |    0.25553937978906    |    0.99375825428829   |
# Volume_BCC    |   0.29043994595705    |    0.39373648391637    |    0.99855177887171   |
# Volume_FCC    |   0.45029070090874    |    0.65941196368284    |    0.98171731349904   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Force Errors  |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.07435126021971    |    0.12848376278885    |    0.96636534845803   |
# Displaced_A15 |   0.09159588698128    |    0.11543137635400    |    0.96709623259060   |
# Displaced_BCC |   0.11631098454386    |    0.14961150381635    |    0.98750172427778   |
# Displaced_FCC |   0.08406970708774    |    0.10543003421090    |    0.95619944089030   |
# Elastic_BCC   |   0.05860645433728    |    0.07186259158171    |    -47819.762694298   |
# Elastic_FCC   |   0.04317009096991    |    0.05538270877199    |    -18726.563076726   |
# GSF_110       |   0.02727025920704    |    0.04553237860642    |    0.87898335977552   |
# GSF_112       |   0.04987378103932    |    0.07922494107237    |    0.82408980275735   |
# Liquid        |   0.27366116024599    |    0.34985442714443    |    0.95494518581465   |
# Surface       |   0.04259395186436    |    0.09563703931568    |    0.82445969671280   |
# Volume_A15    |   3.37240113741997    |    6.15659213169514    |    -Inf000000000000   |
# Volume_BCC    |   3.32863251469136    |    6.98673598115050    |    -Inf000000000000   |
# Volume_FCC    |   1.97839064823913    |    3.85646088762313    |    -Inf000000000000   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Stress Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   46622.9349273975    |    235182.209214730    |    0.98144926194258   |
# Displaced_A15 |   1221.04850726689    |    1553.55605717264    |    0.98563925123290   |
# Displaced_BCC |   7341.61299136141    |    9909.83900242742    |    0.86467137841218   |
# Displaced_FCC |   29864.5485879079    |    41957.0091551672    |    -2.6414119710704   |
# Elastic_BCC   |   425.430825216202    |    504.480804809002    |    0.99619587954770   |
# Elastic_FCC   |   39292.2213685134    |    54711.6947589831    |    -11.003934696752   |
# GSF_110       |   2084.38759763698    |    3183.70255860428    |    0.50643293174301   |
# GSF_112       |   1707.54624727980    |    2859.07610929099    |    0.71666105292148   |
# Liquid        |   42390.0156468890    |    58778.5888928048    |    -6.0401903250365   |
# Surface       |   1662.20520087705    |    2917.99248883842    |    0.90799433870165   |
# Volume_A15    |   70617.2614162872    |    140043.228853478    |    0.99284697772326   |
# Volume_BCC    |   239272.502211254    |    672446.741169593    |    0.98672888320201   |
# Volume_FCC    |   169066.027472050    |    558413.504176687    |    0.94963712438035   |
# ----------------------------------------------------------------------------------------
