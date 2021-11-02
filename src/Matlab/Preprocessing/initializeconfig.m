%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function config = initializeconfig(app)

nconfigs = app.nconfigs;            % number of configurations
dim = app.dim;
% ncx = dim;
% ncv = dim;
% ncf = dim;
% ncq = 1;
% nce = 1;
nci = app.nci;
nct = app.nct;
ncg = app.ncg;
ncx = app.dim;
ncv = app.ncv;
ncf = app.ncf;
ncq = app.ncq;
nce = app.nce;
ncs = app.ncs;
ncz = app.ncz;
ncm = app.ncm;
nco = app.nco;
ncl = app.ncl;
ncp = app.ncp;

config.nconfigs = nconfigs;
config.dim = dim;
config.ncx = ncx;
config.ncv = ncv;
config.ncq = ncq;
config.nce = nce;
config.ncf = ncf;
config.ncs = ncs;
config.nci = nci;
config.nct = nct;
config.ncg = ncg;
config.ncz = ncz;
config.ncm = ncm;
config.nco = nco;
config.ncl = ncl;
config.ncp = ncp;

config.natom = ones(1, nconfigs);   % number of atoms per each configuration
config.lattice = zeros(ncl, nconfigs); % stresses for all configurations
config.pbc = zeros(ncp, nconfigs); % periodic boundary conditions
config.e = ones(nce, nconfigs);  % potential energies for all configurations
config.stress = zeros(ncs, nconfigs); % stresses for all configurations

config.natomall = sum(config.natom);% number of atoms for all configurations
config.we = ones(1, nconfigs);      % energy weight per each configuration
config.wf = ones(1, nconfigs);      % force weight per each configuration
config.ws = ones(1, nconfigs);      % stress weight per each configuration

% simulation box for each configuration
config.a = zeros(3, nconfigs); % the 1st principal vector of the simulation box
config.b = zeros(3, nconfigs); % the 2nd principal vector of the simulation box
config.c = zeros(3, nconfigs); % the 3rd principal vector of the simulation box

config.Z = zeros(ncz, config.natomall);    % atom numbers for all configurations
config.mass = zeros(ncm, config.natomall); % atom masses for all configurations
config.move = zeros(nco, config.natomall);   % atom move masks for all configurations
% config.eatom = zeros(1, config.natomall); % atom energies for all configurations
% config.vatom = zeros(6, config.natomall); % atom virial for all configurations
config.tags = zeros(nci, config.natomall); % atom tags for all configurations
config.t = zeros(nct, config.natomall);   % atom types for all configurations
config.group = zeros(ncg, config.natomall);   % atom groups for all configurations
config.x = zeros(ncx, config.natomall); % atom positions for all configurations
config.q = zeros(ncq, config.natomall); % atom charges for all configurations
config.v = zeros(ncv, config.natomall); % atom velocities for all configurations
config.f = zeros(ncf, config.natomall); % atom forces for all configurations


