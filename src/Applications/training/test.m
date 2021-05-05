version = "0";

cdir = pwd(); 
ii = strfind(cdir, "Applications");
sourcepath = cdir(1:(ii-1));
run(sourcepath + "Installation/setpath.m");

app = initializeapp(sourcepath,version);
%app.cpucompiler = "/usr/local/Cellar/llvm/11.1.0/bin/clang++";
%app.cpuflags = "-Xclang -load -Xclang /usr/local/lib/ClangEnzyme-11.dylib";
app.cpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
app.cpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
app.cpumacros = "-std=c++11 -D _ENZYME";
app.gpucompiler = "/home/linuxbrew/.linuxbrew/bin/clang++";
app.gpuflags = "-Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so";
app.gpumacros = "--cuda-gpu-arch=sm_60 -std=c++11 -D _ENZYME";
app.appname = "test";
app.currentdir = cdir; % current directory
app.configfile = "configs/CONFIGS1"; % configuration file
app.configmode = 4;   % LAMMPS mode for reading configuration file
app.traininglist = 0:79;  % a list of configurations for training the potential
app.validatelist = 0:100;% a list of configurations for validating the potential
app.dim = 3;          % physical dimension
app.maxnumneighbors = 100; % maximum number of neighbors allowed
app.runMD = 0;        % 0 no MD simulation, 1 -> run MD simulation
app.training = 3;     % 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
app.dftdata = 3;      % 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
app.potentialform = 1;% 0 -> empirical potential, 1 -> ML potential, 2 -> hybrid empirical+ML potential    
app.neighpair = 0;    % 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
app.neighcell = 0;    % 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
app.decomposition = 0;% 0 -> force decomposition, 1 -> atom decomposition
app.energycal = 1;    % turns energy calculation on or off
app.forcecal = 1;     % turns force calculation on or off
app.pbc = [1 1 1];    % periodic boundary conditions
app.we = 1;           % energy weight
app.wf = 10;          % force weight

% atom types
app.natomtype = 1;  % number of atom types
app.atomnumbers = [18];
app.atommasses = [39.948];
app.atomcharges = [0];

% Machine learning potential descriptors
app.descriptor = 0;   % descriptor flag: 0 -> Spherical Harmonics Bessel
app.spectrum = 0;     % spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
app.chemtype = 0;     % 0 -> single atom-type basis functions, 1 -> pair atom-type basis functions  
app.K = 6;            % number roots of radial basis functions
app.L = 4;      % degree of radial basis functions and spherical harmonics     
app.rcutml = 8.5;     % cut-off radius for machine learning potential

% Empirical potential descriptors
app.potentialfile = "LJpotential";   % empirical potential file
epsilon = 1.0; %0.01029849;
sigma = 3.4;
A = 4*epsilon*sigma^12;
B = 4*epsilon*sigma^6;
app.pot2a = [1];     % list of nonbonded pair potentials 
app.mu2a = [A B];    % parameters for  LJ potential
app.rcut2a = [8.5];  % cut-off radius 

% train and validate the potential
app = mdp(app);

% get validation data
valid = validate(app);

%setenv('LD_LIBRARY_PATH','/usr/bin');

% /home/linuxbrew/.linuxbrew/bin/clang++ -std=c++11 -D _ENZYME main.cpp -o cpuMDP -ffast-math -O3 -Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so -lblas -llapack
% /home/linuxbrew/.linuxbrew/bin/clang++  --cuda-gpu-arch=sm_60 -std=c++11 -D_FORCE_INLINES -O3 -c -fPIC -w ../Core/gpuCore.cu -o ../Core/gpuCore.o
% /home/linuxbrew/.linuxbrew/bin/clang++  --cuda-gpu-arch=sm_60 -std=c++11 -D_FORCE_INLINES -O3 -c -fPIC -w -Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so ../Potentials/gpuEmpiricalPotentials.cu -o ../Core/gpuEmpiricalPotentials.o
% /home/linuxbrew/.linuxbrew/bin/clang++ --cuda-gpu-arch=sm_60 -std=c++11 -D _DEBUG -D _CUDA main.cpp -o gpuMDP ../Core/gpuCore.o ../Core/gpuEmpiricalPotentials.o -Xclang -load -Xclang /home/cuongng/enzyme/Enzyme/ClangEnzyme-11.so -ffast-math -O2 -fno-vectorize -fno-unroll-loops -fPIC -lcudart -lcublas -lblas -llapack
% cuda-memcheck ./gpuMDP test out
% cuda-gdb --args ./gpuMDP test out

% mkdir build
% cd build
% cmake -D CMAKE_C_COMPILER=/home/linuxbrew/.linuxbrew/bin/clang -D CMAKE_CXX_COMPILER=/home/linuxbrew/.linuxbrew/bin/clang++ .. -DLLVM_DIR=/home/linuxbrew/.linuxbrew/opt/llvm@11/lib/cmake/clang
% cd ..
% make
