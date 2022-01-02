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
descriptors[2] = initsna(snaparam, elemradius, elemweight)

# ZBL reference potential
potentials[1] = addpotential("pair", "nonbonded", [1, 1], "zbl", [4.0, 4.8, 73], 4.8);
include(potentials[1].potentialfunction * ".jl");  # include the potential file

# loss function style must be: "energy", "force", "energyforce", "energyforcestress"
lossfunc = "energyforcestress"
# use least-square method
method = "lsq" 
# divide energies by the number of atoms
normalizeenergy = true
# divide virials by the box volume
normalizestress = true
# optimization parameters
optim = setoptim(lossfunc, method, normalizeenergy, normalizestress)

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
# ALL           |   0.11265174691267    |    0.37977565433503    |    0.99419038677079   |
# Displaced_A15 |   0.00235792938939    |    0.00246855977428    |    -1.2156363735745   |
# Displaced_BCC |   0.00188371436433    |    0.00273064614363    |    0.99587522319726   |
# Displaced_FCC |   0.00061741229597    |    0.00075009168170    |    0.97160166434659   |
# Elastic_BCC   |   0.00599579098947    |    0.00604563796454    |    0.99230919870267   |
# Elastic_FCC   |   0.00354421105172    |    0.00445843637691    |    0.99333159772842   |
# GSF_110       |   0.00652425814607    |    0.00883023587747    |    0.77951696510583   |
# GSF_112       |   0.00868251914489    |    0.00961453139000    |    0.87318809800944   |
# Liquid        |   0.00566811969259    |    0.00655584081006    |    0.99601358415720   |
# Surface       |   0.01352360826034    |    0.01571242406414    |    0.98816795126889   |
# Volume_A15    |   0.24578877837149    |    0.34490495333075    |    0.99781261498874   |
# Volume_BCC    |   0.32857920069640    |    0.51271284190945    |    0.99853838108547   |
# Volume_FCC    |   0.81209055282782    |    1.18123504852115    |    0.98443646179057   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Force Errors  |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   0.08359239990607    |    0.16311584659554    |    0.98512903562289   |
# Displaced_A15 |   0.09806643270309    |    0.12820169392009    |    0.94957955563905   |
# Displaced_BCC |   0.14339589323270    |    0.18486360986945    |    0.98517811287921   |
# Displaced_FCC |   0.05740637316446    |    0.07166197131771    |    0.98061307246834   |
# Elastic_BCC   |   0.06439647995478    |    0.07899854689009    |    0.41298668776769   |
# Elastic_FCC   |   0.05065258662499    |    0.06490610225971    |    0.63719969979068   |
# GSF_110       |   0.02651057125534    |    0.05418215537551    |    0.99889872481223   |
# GSF_112       |   0.05698006876321    |    0.10025780716990    |    0.99737013744671   |
# Liquid        |   0.35379720579374    |    0.48364455518089    |    0.94779367837188   |
# Surface       |   0.04760235280733    |    0.10895743153042    |    0.99671811972583   |
# Volume_A15    |   2.79473115136733    |    6.25088426150716    |    0.69215396831169   |
# Volume_BCC    |   2.88326079591402    |    6.06857960431047    |    -4.7705594018581   |
# Volume_FCC    |   1.39818002339307    |    2.88854746833441    |    0.32360003124122   |
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Stress Errors |          MAE          |          RMSE          |          RSQ          |
# ----------------------------------------------------------------------------------------
# ALL           |   68331.2692226880    |    381756.218018856    |    0.93961426214254   |
# Displaced_A15 |   19553.8208834055    |    26783.3681374323    |    0.99755940235506   |
# Displaced_BCC |   14933.0136738189    |    20558.4972892687    |    0.99870935725821   |
# Displaced_FCC |   35323.5274494020    |    49684.5017070339    |    0.99113396279263   |
# Elastic_BCC   |   170.405533015481    |    268.553751189739    |    0.99999977908789   |
# Elastic_FCC   |   48132.1295201740    |    66784.2466913797    |    0.98419798844541   |
# GSF_110       |   1623.00572832491    |    2322.56048399292    |    0.99993839179953   |
# GSF_112       |   2068.90488220449    |    3116.39773656099    |    0.99987644651027   |
# Liquid        |   42920.1434033482    |    58562.0314191833    |    0.98628360736784   |
# Surface       |   2309.28255361132    |    4277.50368379242    |    0.99977822976790   |
# Volume_A15    |   132939.655632498    |    357432.429798387    |    0.95261560122466   |
# Volume_BCC    |   271466.858168366    |    770712.743079258    |    0.97424571846356   |
# Volume_FCC    |   304211.704271920    |    1.07923030094437    |    0.79235217127926   |
# ----------------------------------------------------------------------------------------
