%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function config = initializeconfig(app)

nconfigs = app.nconfigs;            % number of configurations
dim = app.dim;
ncx = app.ncx;
ncv = app.ncv;
ncf = app.ncf;
ncq = app.ncq;
nce = app.nce;

config.dim = dim;
config.ncx = ncx;
config.ncv = ncv;
config.ncq = ncq;
config.nce = nce;
config.ncf = ncf;
config.nconfigs = nconfigs;
config.natom = ones(1, nconfigs);   % number of atoms per each configuration
config.natomall = sum(config.natom);% number of atoms for all configurations

% simulation box for each configuration
config.a = zeros(dim, nconfigs); % the 1st principal vector of the simulation box
config.b = zeros(dim, nconfigs); % the 2nd principal vector of the simulation box
config.c = zeros(dim, nconfigs); % the 3rd principal vector of the simulation box
config.e = ones(nce, nconfigs);  % energies for all configurations

config.t = zeros(1, config.natomall);   % atom types for all configurations
config.x = zeros(ncx, config.natomall); % atom positions for all configurations
config.q = zeros(ncq, config.natomall); % atom charges for all configurations
config.v = zeros(ncv, config.natomall); % atom velocities for all configurations
config.f = zeros(ncf, config.natomall); % atom forces for all configurations


