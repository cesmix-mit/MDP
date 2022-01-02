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

# bounds for scaling parameleter to be optimized
scalerange = [3.0, 7.0]

# define range for nonlinear parameters to be optimized
etarange =  vcat(reshape(rcutrange,(1,2)), reshape(scalerange,(1,2)))

# number of interpolation points for each nonlinear parameters 
N = [6, 10]
etaspace = [etarange N]

# number of descriptors
kappa = [5];      

# define a list of empirical and ML descriptors
descriptors = Array{Any}(undef, 1)
descriptors[1] = Potential.PotentialStruct("pair", "nonbonded", [1, 1], "scaledbessel", [0.0], rcutrange[1]);
for i = 1:length(descriptors)
    include(descriptors[i].potentialfunction * ".jl");  # include the descriptors file
end

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforcestress"

# optimize the potential
eta, coeff, fmin, iter, polycoeff = Optimization.optimize(traindata, testdata, 
                    descriptors, etaspace, kappa, weightouter, lossfunc, pbc, a, b, c)    

# verify the optimal solution against the minimum point on a fine grid 
etagrid, fgrid = Optimization.mpolyeval([etaspace [36, 54]], polycoeff)
fout, ind = findmin(fgrid)
etaout = etagrid[ind,:]

print("                Gradient-descent optimization        Fine-grid search  \n");
for i = 1:length(eta)
    etai = eta[i]
    etaouti = etaout[i]
    print("parameter $i:       $etai                  $etaouti  \n");
end
print("objective value:   $fmin              $fout  \n");
                    
# validate the Bessel potential 
emae, fmae, smae = Optimization.validate(validdata, descriptors, coeff, eta, kappa, pbc, a, b, c)
print("     Mean absolute errors on the validation data set  \n");
print("Energy MAE:  $emae  \n");
print("Force  MAE:  $fmae \n");
print("Stress MAE:  $smae \n");

# Read configuration from /Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster
#                 Gradient-descent optimization        Fine-grid search  
# parameter 1:       3.460806619112657                  3.4643629203327673  
# parameter 2:       5.855284463760658                  5.882802215561238  
# objective value:   0.004451193410442596              0.00445912531843542  
# Read configuration from /Users/ngoccuongnguyen/Dropbox (MIT)/MDP/data/LJCluster
#      Mean absolute errors on the validation data set  
# Energy MAE:  0.002377523806746342  
# Force  MAE:  0.005204441463953111 
# Stress MAE:  0.0122515068463677 

