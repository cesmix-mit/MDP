function [app,config] = preprocessing(app,config)

app.flag = [app.descriptor app.spectrum app.training app.runMD app.potentialform app.neighpair...
            app.energycal app.forcecal app.stresscal app.neighcell app.decomposition app.chemtype];

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

m = 1;
for i = 1:config.nconfigs
    [B2C, C2B] = cubemapping(config.a(:,i), config.b(:,i), config.c(:,i));
    ximages = boxperiodicimages(app.pbc, config.a(:,i), config.b(:,i), config.c(:,i));    
    n = config.natom(i);
    config.x(:,m:(m+n-1)) = checkconfig(config.x(:,m:(m+n-1)), ximages, B2C, C2B);    
    m = m + n;
end

writeapp(app, app.appname + "app.bin");
writeconfig(config, app.appname + "config.bin", 2);
               
end

