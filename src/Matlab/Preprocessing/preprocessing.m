function [app,config] = preprocessing(app)

% read configurations from data file
config = readconfig(app, app.configmode);     

 % number of configurations
app.nconfigs = config.nconfigs;          

% nonbonded single potentials
app.npot1a = length(app.pot1a);  % number of nonbonded single potentials 
app.ncmu1a = length(app.mu1a); % length of mu1a

% bonded single potentials
app.npot1b = length(app.pot1b);  % number of nonbonded single potentials 
app.ncmu1b = length(app.mu1b); % length of mu1a

% nonbonded pair potentials
app.npot2a = length(app.pot2a);  % number of nonbonded single potentials 
app.ncmu2a = length(app.mu2a); % length of mu1a

% bonded pair potentials
app.npot2b = length(app.pot2b);  % number of nonbonded single potentials 
app.ncmu2b = length(app.mu2b); % length of mu1a

% two-body bond order potentials
app.npot2c = length(app.pot2c);  % number of nonbonded single potentials 
app.ncmu2c = length(app.mu2c); % length of mu1a

% nonbonded triplet potentials
app.npot3a = length(app.pot3a);  % number of nonbonded single potentials 
app.ncmu3a = length(app.mu3a); % length of mu1a

% bonded triplet potentials
app.npot3b = length(app.pot3b);  % number of nonbonded single potentials 
app.ncmu3b = length(app.mu3b); % length of mu1a

% three-body bond order potentials
app.npot3c = length(app.pot3c);  % number of nonbonded single potentials 
app.ncmu3c = length(app.mu3c); % length of mu1a

% nonbonded quadruplet potentials
app.npot4a = length(app.pot4a);  % number of nonbonded single potentials 
app.ncmu4a = length(app.mu4a); % length of mu1a

% bonded quadruplet potentials
app.npot4b = length(app.pot4b);  % number of nonbonded single potentials 
app.ncmu4b = length(app.mu4b); % length of mu1a

% app.ncx = 0; % number of compoments of x
% app.ncv = 0; % number of compoments of v
% app.nce = 0; % number of compoments of e
% app.ncf = 0; % number of compoments of f
% app.ncq = 0; % number of compoments of q
% app.ncmu = 0; % number of compoments of mu
% app.nceta = 0; % number of compoments of eta
% app.nckappa = 0; % number of compoments of kappa

app.flag = [app.descriptor app.spectrum app.training app.runMD app.potentialform app.neighpair...
            app.energycal app.forcecal app.stresscal app.neighcell app.decomposition app.chemtype ....
            app.dftdata];

rcut = [app.rcutml; app.rcut2a(:); app.rcut2b(:); app.rcut2c(:); ...
              app.rcut3a(:); app.rcut3b(:); app.rcut3c(:); app.rcut4a(:); app.rcut4b(:)];
rcutmax = max(rcut);          
app.rcutsqmax = max(rcutmax.^2);        
app.boxoffset = [rcutmax rcutmax rcutmax];
app.simparam = [app.time app.dt];        

app.ndims = zeros(20,1);
app.ndims(1) = app.dim;
app.ndims(2) = app.L;
app.ndims(3) = app.K;
app.ndims(4) = app.ntimesteps;
app.ndims(5) = app.nab;
app.ndims(6) = app.mpiprocs;  
app.ndims(7) = app.backend;
app.ndims(8) = app.nconfigs;
app.ndims(9) = app.natomtype;
app.ndims(10) = app.nmoletype; 
app.ndims(11) = app.maxnumneighbors;

% check configurations
m = 1;
for i = 1:config.nconfigs
    [B2C, C2B] = cubemapping(config.a(:,i), config.b(:,i), config.c(:,i));
    ximages = boxperiodicimages(app.pbc, config.a(:,i), config.b(:,i), config.c(:,i));    
    n = config.natom(i);
    config.x(:,m:(m+n-1)) = checkconfig(config.x(:,m:(m+n-1)), ximages, B2C, C2B);    
    m = m + n;
end

if app.training == 0
    config.we = [];
    config.wf = [];
else
    config = readweight(app, config);    
end

app = writeapp(app, app.appname + "app.bin");
writeconfig(config, app.appname + "config.bin");
   
movefile(char(app.appname + "app.bin"), char(app.sourcepath + "C++/Main"));
movefile(char(app.appname + "config.bin"), char(app.sourcepath + "C++/Main"));

end

