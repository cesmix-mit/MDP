function config = readMDP(filename, mode)

fileID = fopen(filename,'r');
if mode==0 % binary    
    frame = fread(fileID,'double');
else % text    
    frame = fscanf(fileID, '%f');
end
fclose(fileID);

natomall = 0;
nframe = 0;
n2 = 0;
n = length(frame);
while (n2 < n)
    nframe = nframe + 1;
    
    % id, type, group, move, atomic number, atomic mass, charge, position, velocity, force
    natom = int64(frame(n2+1));  % number of atoms
    ncl = int64(frame(n2+2));    % lattice  
    ncp = int64(frame(n2+3));    % periodic boundary conditions     
    nce = int64(frame(n2+4));    % energy
    ncs = int64(frame(n2+5));    % stress 
    nci = int64(frame(n2+6));    % atom id
    nct = int64(frame(n2+7));    % atom type
    ncg = int64(frame(n2+8));    % atom groups
    nco = int64(frame(n2+9));    % atom moved/fixed flags
    ncz = int64(frame(n2+10));   % atomic numbers
    ncm = int64(frame(n2+11));   % atomic masses     
    ncq = int64(frame(n2+12));   % atom charges
    ncx = int64(frame(n2+13));   % positions
    ncv = int64(frame(n2+14));   % atom velocities 
    ncf = int64(frame(n2+15));   % atom forces
    n2 = n2+15;
    
    if nframe== 1
        app.ncx = ncx;
        app.ncv = ncv;
        app.ncq = ncq;
        app.nce = nce;
        app.ncf = ncf;
        app.ncs = ncs;
        app.nci = nci;
        app.nct = nct;
        app.ncg = ncg;
        app.ncz = ncz;
        app.ncm = ncm;
        app.nco = nco;
        app.ncl = ncl;
        app.ncp = ncp;
        app.nconfigs = 1; 
        app.dim = ncx;
        config = initializeconfig(app);    
    end
    
    config.natom(nframe) = natom;
    
    % lattice, pbc, energy, stress
    n1 = n2 + 1; n2 = n2 + ncl;
    if ncl>0    
        config.lattice(:,nframe) = frame(n1:n2);        
    end

    n1 = n2 + 1; n2 = n2 + ncp;
    if ncp>0
        config.pbc(:,nframe) = frame(n1:n2);        
    end

    n1 = n2 + 1; n2 = n2 + nce;
    if nce>0
        config.e(:,nframe) = frame(n1:n2);        
    end

    n1 = n2 + 1; n2 = n2 + ncs;
    if ncs>0
        config.stress(:,nframe) = frame(n1:n2);        
    end
    
    k = nci + nct + ncg + nco + ncz + ncm + ncq + ncx + ncv + ncf;
    n1 = n2 + 1; n2 = n2 + natom*k;    
    tmp = reshape(frame(n1:n2), [k natom]);            
    if nci>0    
        config.tags((natomall+1):(natomall+natom)) = tmp(1:nci,:);        
    end
    if nct>0            
        config.t((natomall+1):(natomall+natom)) = tmp((1+nci):(nci+nct),:);               
    end    
    if ncg>0
        config.group((natomall+1):(natomall+natom)) = tmp((1+nci+nct):(nci+nct+ncg),:);        
    end
    if nco>0
        config.move((natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg):(nci+nct+ncg+nco),:);        
    end
    if ncz>0
        config.Z((natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg+nco):(nci+nct+ncg+nco+ncz),:);        
    end
    if ncm>0
        config.mass((natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg+nco+ncz):(nci+nct+ncg+nco+ncz+ncm),:);        
    end
    if ncq>0
        config.q(:,(natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg+nco+ncz+ncm):(nci+nct+ncg+nco+ncz+ncm+ncq),:);        
    end
    if ncx>0
        config.x(:,(natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg+nco+ncz+ncm+ncq):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx),:);        
    end
    if ncv>0
        config.v(:,(natomall+1):(natomall+natom)) = tmp((1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv),:);        
    end
    if ncf>0
        config.f(:,(natomall+1):(natomall+natom)) = tmp(:,(1+nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv):(nci+nct+ncg+nco+ncz+ncm+ncq+ncx+ncv+ncf))';        
    end
    
    natomall = natomall + natom;
end
config.natomall = natomall;




