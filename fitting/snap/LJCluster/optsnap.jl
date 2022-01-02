cdir = pwd(); ii = findlast("MDP", cdir); MDPpath = cdir[1:ii[end]] * "/";    
include(MDPpath * "Installation/setpath.jl");

# path to the database 
datapath = "/Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster";
dataformat = "extxyz"
fileextension = "xyz"
atomspecies = ["Ar"];

# randomly selecting the configurations in the database
randomize = true;

# weights for energies, forces, and stresses in the inner linear fit
weightinner = [1.0; 1.0; 0.25]

# weights for energies, forces, and stresses in the outer nonlinear optimization
weightouter = [0.2, 0.8, 0.005]

# shift atom positions to put them in the box [0,10]x[0,10]x[0,10]
translation = [5.0; 5.0; 5.0]

# create data sets for training, testing and validation
# data set to compute the loss function in the inner linear optimization
percentage = 10.0;
traindata = Array{Any}(undef, 1)
traindata[1] = Preprocessing.adddata(datapath, dataformat, fileextension, percentage, 
                        randomize, atomspecies, weightinner, translation)

# data set to compute the loss function in the outer nonlinear optimization
percentage = 10.0;           
testdata = Array{Any}(undef, 1)             
testdata[1] = Preprocessing.adddata(datapath, dataformat, fileextension, percentage, 
                         randomize, atomspecies, weightinner, translation)

# data set to validate the fitted bessel potential after training and testing
percentage = 20.0;
validdata = Array{Any}(undef, 1)       
validdata[1] = Preprocessing.adddata(datapath, dataformat, fileextension, percentage, 
                        randomize, atomspecies, weightinner, translation)

# 3 vectors define the simulation box
a = [10.0; 0.0; 0.0]
b = [0.0; 10.0; 0.0]
c = [0.0; 0.0; 10.0]

# periodic boundary conditions
pbc = [0; 0; 0]

# bounds for cut-off radius to be optimized
rcutrange = [2.5, 3.5]

# define range for nonlinear parameters to be optimized
etarange =  reshape(rcutrange,(1,2))

# number of interpolation points for each nonlinear parameters 
N = [10]
etaspace = [etarange N]

# number of empirical descriptors
kappa = [0];      

# SNAP parameters
twojmax = 8
ntypes = 1
rcutfac = 3.0;  
rfac0 = 1.0;
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
descriptors = Array{Any}(undef, 1)
descriptors[1] = Potential.initsna(snaparam, elemradius, elemweight)

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforcestress"

# optimize the SNAP potential
eta, coeff, fmin, iter, polycoeff = Optimization.optimize(traindata, testdata, 
                    descriptors, etaspace, kappa, weightouter, lossfunc, pbc, a, b, c)    
                                      
# verify the optimal solution against the minimum point on a fine grid 
etagrid, fgrid = Optimization.mpolyeval([etaspace [100]], polycoeff)
fout, ind = findmin(fgrid)
etaout = etagrid[ind,:]

print("                Gradient-descent optimization        Fine-grid search  \n");
for i = 1:length(eta)
    etai = eta[i]
    etaouti = etaout[i]
    print("parameter $i:       $etai                  $etaouti  \n");
end
print("objective value:   $fmin              $fout  \n");
                    
# validate the potential 
emae, fmae, smae = Optimization.validate(validdata, descriptors, coeff, eta, kappa, pbc, a, b, c)
print("     Mean absolute errors on the validation data set  \n");
print("Energy MAE:  $emae  \n");
print("Force  MAE:  $fmae \n");
print("Stress MAE:  $smae \n");


#                 Gradient-descent optimization        Fine-grid search  
# parameter 1:       3.4874049149170774                  3.4865147728462356  
# objective value:   0.585635934018706              0.5856682939000912  
# Read configuration from /Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster
#      Mean absolute errors on the validation data set  
# Energy MAE:  0.21826706448821329  
# Force  MAE:  0.6892192811999076 
# Stress MAE:  1.5131521670494448 
