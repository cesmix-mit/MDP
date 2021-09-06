%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [app,config] = preprocessing(app)

% read configurations from data file
if strcmpi(app.configfile,"")==0 && app.configmode>0
    config = readconfig(app, app.configmode);     
    % number of configurations
    app.nconfigs = config.nconfigs;          
else    
    app.nconfigs = 0;          
    config = 0;
end

if strcmpi(app.snapcoefffile,"")==0 
    fileID = fopen(app.snapcoefffile);
    coeff = textscan(fileID,'%f','Delimiter',',');
    app.snapcoeff = coeff{1};
    fclose(fileID);   
end

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

 %lj or real or metal or si or cgs or electron or micro or nano
style = app.unitstyle;
if (strcmpi(style,"lj")) 
    app.unitstyle = 0;
elseif (strcmpi(style,"real")) 
    app.unitstyle = 1;    
elseif (strcmpi(style,"metal")) 
    app.unitstyle = 2;
elseif (strcmpi(style,"si")) 
    app.unitstyle = 3;
elseif (strcmpi(style,"cgs")) 
    app.unitstyle = 4;
elseif (strcmpi(style,"electron")) 
    app.unitstyle = 5;
elseif (strcmpi(style,"micro")) 
    app.unitstyle = 6;
elseif (strcmpi(style,"nano")) 
    app.unitstyle = 7;
else
    error("Invalid unit style");
end

% descriptor flag: 0 -> Spherical Harmonics Bessel, 1-> snap
descriptor = app.descriptor;
if (strcmpi(descriptor,"shp")) 
    app.descriptor = 0;
elseif (strcmpi(descriptor,"snap")) 
    app.descriptor = 1;    
else
    app.descriptor = -1; 
end

ensemblemode = app.ensemblemode;
if (strcmpi(ensemblemode,"nve")) 
    app.ensemblemode = 0;
    app.runMD = 1;  
elseif (strcmpi(ensemblemode,"nvelimit")) 
    app.ensemblemode = 1;    
    app.runMD = 1;  
elseif (strcmpi(ensemblemode,"nvt")) 
    app.ensemblemode = 2;        
    app.runMD = 1;  
else
    app.ensemblemode = -1;
    app.runMD = 0;  
end

app.natomtype = length(app.atommasses);
app.atomnumbers = [0; app.atomnumbers(:)];
app.atommasses = [0; app.atommasses(:)];
app.atomcharges = [0; app.atomcharges(:)];

app.flag = [app.descriptor app.spectrum app.training app.runMD app.potentialform app.neighpair...
            app.energycal app.forcecal app.stresscal app.neighcell app.decomposition app.chemtype ....
            app.dftdata app.unitstyle app.ensemblemode app.neighcheck];

rcut = [app.rcutml; app.rcut2a(:); app.rcut2b(:); app.rcut2c(:); ...
              app.rcut3a(:); app.rcut3b(:); app.rcut3c(:); app.rcut4a(:); app.rcut4b(:)];
rcutmax = max(rcut)+app.neighskin;          
app.rcutsqmax = max(rcutmax.^2);        
app.boxoffset = [rcutmax rcutmax rcutmax];
app.simparam = [app.time app.dt];        
app.solparam = [rcutmax app.neighskin];
app.snaparam = [app.snapnelem app.snapncoeff app.snaptwojmax app.snaprcutfac app.snaprfac0 app.snaprmin0 ...
                app.snapbzeroflag app.snapswitchflag app.snapquadraticflag app.snapchemflag ...
                app.snapbnormflag app.snapwselfallflag];
            
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
app.ndims(11) = app.neighmax;
app.ndims(12) = app.neighevery;
app.ndims(13) = app.neighdelay;
app.ndims(14) = app.globalfreq;
app.ndims(15) = app.peratomfreq;
 
% check configurations
m = 1;
for i = 1:app.nconfigs
    [B2C, C2B] = cubemapping(config.a(:,i), config.b(:,i), config.c(:,i));
    ximages = boxperiodicimages(app.pbc, config.a(:,i), config.b(:,i), config.c(:,i));    
    n = config.natom(i);
    config.x(:,m:(m+n-1)) = checkconfig(config.x(:,m:(m+n-1)), ximages, B2C, C2B);    
    m = m + n;
end

%bindir = "C++/Main";
bindir = "exec";
if ~isempty(app.lattice)
    writelattice(app.lattice, app.appname + "lattice.bin");
    movefile(char(app.appname + "lattice.bin"), char(app.sourcepath + bindir));
end
if ~isempty(app.region)
    writeregion(app.region, app.appname + "region.bin");
    movefile(char(app.appname + "region.bin"), char(app.sourcepath + bindir));
end
if ~isempty(app.domain)
    writedomain(app.domain, app.appname + "domain.bin");
    movefile(char(app.appname + "domain.bin"), char(app.sourcepath + bindir));
end

if app.nconfigs>0
    if app.training == 0
        config.we = [];
        config.wf = [];
    else
        config = readweight(app, config);    
    end    
    writeconfig(config, app.appname + "config.bin");
    movefile(char(app.appname + "config.bin"), char(app.sourcepath + bindir));
end   

app = writeapp(app, app.appname + "app.bin");
movefile(char(app.appname + "app.bin"), char(app.sourcepath + bindir));

end

