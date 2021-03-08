function config = initializeconfig(app)

nconfigs = app.nconfigs;            % number of configurations
dim = app.dim;
ncx = app.ncx;
ncv = app.ncv;
ncf = app.ncf;
ncq = app.ncq;

config.natom = ones(1, nconfigs);   % number of atoms per each configuration
config.natomall = sum(config.natom);% number of atoms for all configurations

% simulation box for each configuration
config.a = zeros(dim, nconfigs); % the 1st principal vector of the simulation box
config.b = zeros(dim, nconfigs); % the 2nd principal vector of the simulation box
config.c = zeros(dim, nconfigs); % the 3rd principal vector of the simulation box

config.t = zeros(1, config.natomall);
config.x = zeros(ncx, config.natomall);
config.q = zeros(ncq, config.natomall);
config.v = zeros(ncv, config.natomall);
config.f = zeros(ncf, config.natomall);


