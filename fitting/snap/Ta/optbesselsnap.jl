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

rcut = 5.0   # cut-off radiuss 

# single descriptors to catpure energy of isolated pure elements
descriptors[1] = PotentialStruct("single", "nonbonded", [1], "single", [0.0], rcut);
include(descriptors[1].potentialfunction * ".jl");  # include the descriptors file
      
# Bessel descriptors
besselscale = 2.5;
eta = [rcut; besselscale]; 
kappa = [6];      # number of Bessel descriptors
descriptors[2] = Potential.PotentialStruct("pair", "nonbonded", [1, 1], "scaledbessel", [0.0], rcut);
include(descriptors[2].potentialfunction * ".jl");  # include the descriptors file

# SNAP parameters
twojmax = 6
ntypes = 1
rcutfac = 4.67637;  
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

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.25, 0.75, 2e-7]

# bounds for cut-off radius to be optimized
rcutrange = [4.6, 5.4]

# bounds for scaling parameleter to be optimized
scalerange = [2.0, 4.0]

# define range for nonlinear parameters to be optimized
etarange =  vcat(reshape(rcutrange,(1,2)), reshape(scalerange,(1,2)))

# number of interpolation points for each nonlinear parameters 
N = [8, 8]
etaspace = [etarange N]

# optimization parameters
optim = setoptim(lossfunc, method, normalizeenergy, normalizestress, eta, kappa, weightouter, etaspace)

# optimize the Bessel/SNAP potential 
eta, coeff, fmin, iter, polycoeff = optimize2(traindata, traindata, descriptors, potentials, optim, pbc)

print("SNAP Coeffficients: "), show(stdout, "text/plain", coeff)

# compute unweighted MAE, RMSE, RSQ errors 
energyerrors, forceerrors, stresserrors = validate2(traindata, descriptors, potentials, optim, coeff, pbc)    

printerrors(folders, energyerrors, "Energy Errors")
printerrors(folders, forceerrors, "Force Errors")
printerrors(folders, stresserrors, "Stress Errors")

using DelimitedFiles
etaspace = [etarange N [32, 32]]
etaout, fout = mpolyeval(etaspace, polycoeff)
out = [etaout fout]
writedlm("optbesselsnaplossfunction.txt", out)


# ----------------------------------------------------------------------------------------
# Energy Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.04407201624074    |    0.13895234006755    |    0.99848863618353   |
# Displaced_A15 |   0.00064295683681    |    0.00084732580228    |    0.90326784300319   |
# Displaced_BCC |   0.00266815814537    |    0.00279961784024    |    0.99672203251388   |
# Displaced_FCC |   0.00159030140295    |    0.00172387396161    |    0.81493658950080   |
# Elastic_BCC   |   0.01091363883685    |    0.01100632694440    |    -588.43048075524   |
# Elastic_FCC   |   0.00302729308650    |    0.00306579881163    |    -1.3064502487690   |
# GSF_110       |   0.00366406390478    |    0.00400799472712    |    0.88510147584172   |
# GSF_112       |   0.00325432404242    |    0.00354119150458    |    0.96659981570708   |
# Liquid        |   0.00084617992731    |    0.00103532130901    |    0.99987144804575   |
# Surface       |   0.01079599316896    |    0.01500070756280    |    0.94230931245908   |
# Volume_A15    |   0.07848650524447    |    0.10962449687720    |    0.99885130263424   |
# Volume_BCC    |   0.24761485607900    |    0.34764026600377    |    0.99887102660169   |
# Volume_FCC    |   0.21855205786430    |    0.36344650749071    |    0.99444598091703   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Force Errors  |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.07363846374494    |    0.12575317788972    |    0.96777978926350   |
# Displaced_A15 |   0.09446578972838    |    0.11852562750111    |    0.96530855378974   |
# Displaced_BCC |   0.10660973760926    |    0.13590313156901    |    0.98968714163891   |
# Displaced_FCC |   0.09629424008693    |    0.12054105144872    |    0.94274401239914   |
# Elastic_BCC   |   0.05808574102001    |    0.07122021832773    |    -46968.652761763   |
# Elastic_FCC   |   0.03751218866409    |    0.04863178583668    |    -14439.203632284   |
# GSF_110       |   0.03187074058294    |    0.05363994853004    |    0.83204957019260   |
# GSF_112       |   0.04774908683784    |    0.07436887168531    |    0.84499362823215   |
# Liquid        |   0.26273384173361    |    0.33833743022054    |    0.95786271517765   |
# Surface       |   0.04000342054590    |    0.09161969982319    |    0.83889748233340   |
# Volume_A15    |   2.89993087276951    |    4.85293370653430    |    -Inf000000000000   |
# Volume_BCC    |   2.76832107868402    |    5.25683612049748    |    -Inf000000000000   |
# Volume_FCC    |   2.09738287182340    |    4.06096724650454    |    -Inf000000000000   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Stress Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   16575.0798602557    |    87940.0330737213    |    0.99740626554616   |
# Displaced_A15 |   15490.4558881970    |    21516.7061763517    |    -1.7547031528251   |
# Displaced_BCC |   11778.4200009532    |    16358.5489271636    |    0.63123812456294   |
# Displaced_FCC |   915.502519376297    |    1232.52113746799    |    0.99685768641770   |
# Elastic_BCC   |   680.783042315670    |    906.291087830642    |    0.98772277016923   |
# Elastic_FCC   |   8359.52276223574    |    11754.3177725659    |    0.44593828295390   |
# GSF_110       |   1828.42673919399    |    2901.92378103316    |    0.58993452570032   |
# GSF_112       |   1754.88638599766    |    2747.82369773695    |    0.73828261524334   |
# Liquid        |   43409.1051199719    |    60496.6513110316    |    -6.4577661295236   |
# Surface       |   2344.03298702141    |    4224.36580573058    |    0.80717237641878   |
# Volume_A15    |   25323.0104544580    |    56682.3988563785    |    0.99882817761392   |
# Volume_BCC    |   93126.9641553286    |    232208.518467598    |    0.99841748447201   |
# Volume_FCC    |   61878.6776523236    |    223359.641907389    |    0.99194235547032   |
# ----------------------------------------------------------------------------------------
