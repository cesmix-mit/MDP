function writeapp(app,filename)

    nsize = zeros(50,1);
    nsize[1] = length(app.ndims[:]);
    nsize[2] = length(app.flag[:]);  
    nsize[3] = length(app.bcs[:]); 
    nsize[4] = length(app.pbc[:]); 
    nsize[5] = length(app.boxoffset[:]); 
    nsize[6] = length(app.atomnumbers[:]);
    nsize[7] = length(app.atommasses[:]);
    nsize[8] = length(app.atomcharges[:]);
    nsize[9] = length(app.simparam[:]); 
    nsize[10] = length(app.solparam[:]); 
    nsize[11] = length(app.eta[:]);
    nsize[12] = length(app.kappa[:]);
    nsize[13] = length(app.muml[:]); 
    nsize[14] = length(app.mu1a[:]);
    nsize[15] = length(app.mu1b[:]);
    nsize[16] = length(app.mu2a[:]);
    nsize[17] = length(app.mu2b[:]);
    nsize[18] = length(app.mu2c[:]);
    nsize[19] = length(app.mu3a[:]);
    nsize[20] = length(app.mu3b[:]);
    nsize[21] = length(app.mu3c[:]);
    nsize[22] = length(app.mu4a[:]);
    nsize[23] = length(app.mu4b[:]);
    nsize[24] = length(app.pot1a[:]);
    nsize[25] = length(app.pot1b[:]);
    nsize[26] = length(app.pot2a[:]);
    nsize[27] = length(app.pot2b[:]);
    nsize[28] = length(app.pot2c[:]);
    nsize[29] = length(app.pot3a[:]);
    nsize[30] = length(app.pot3b[:]);
    nsize[31] = length(app.pot3c[:]);
    nsize[32] = length(app.pot4a[:]);
    nsize[33] = length(app.pot4b[:]);
    nsize[34] = length(app.rcutml);
    nsize[35] = length(app.rcut2a[:]);
    nsize[36] = length(app.rcut2b[:]);
    nsize[37] = length(app.rcut2c[:]);
    nsize[38] = length(app.rcut3a[:]);
    nsize[39] = length(app.rcut3b[:]);
    nsize[40] = length(app.rcut3c[:]);
    nsize[41] = length(app.rcut4a[:]);
    nsize[42] = length(app.rcut4b[:]);
    nsize[43] = length(app.rcutsqmax);
    nsize[44] = length(app.atom1b[:]);
    nsize[45] = length(app.atom2b[:]);
    nsize[46] = length(app.atom2c[:]);
    nsize[47] = length(app.atom3b[:]);
    nsize[48] = length(app.atom3c[:]);
    nsize[49] = length(app.atom4b[:]);


    #app.nsize = nsize;
    fileID = open(filename,"w");
    write(fileID,Float64(length(nsize[:])));
    write(fileID,Float64.(nsize[:]));
    write(fileID,Float64.(app.ndims[:]));
    write(fileID,Float64.(app.flag[:]));
    write(fileID,Float64.(app.bcs[:]));
    write(fileID,Float64.(app.pbc[:]));
    write(fileID,Float64.(app.boxoffset[:]));
    write(fileID,Float64.(app.atomnumbers[:]));
    write(fileID,Float64.(app.atommasses[:]));
    write(fileID,Float64.(app.atomcharges[:]));
    write(fileID,Float64.(app.simparam[:]));
    write(fileID,Float64.(app.solparam[:]));
    write(fileID,Float64.(app.eta[:]));
    write(fileID,Float64.(app.kappa[:]));
    write(fileID,Float64.(app.muml[:]));
    write(fileID,Float64.(app.mu1a[:]));
    write(fileID,Float64.(app.mu1b[:]));
    write(fileID,Float64.(app.mu2a[:]));
    write(fileID,Float64.(app.mu2b[:]));
    write(fileID,Float64.(app.mu2c[:]));
    write(fileID,Float64.(app.mu3a[:]));
    write(fileID,Float64.(app.mu3b[:]));
    write(fileID,Float64.(app.mu3c[:]));
    write(fileID,Float64.(app.mu4a[:]));
    write(fileID,Float64.(app.mu4b[:]));
    write(fileID,Float64.(app.pot1a[:]));
    write(fileID,Float64.(app.pot1b[:]));
    write(fileID,Float64.(app.pot2a[:]));
    write(fileID,Float64.(app.pot2b[:]));
    write(fileID,Float64.(app.pot2c[:]));
    write(fileID,Float64.(app.pot3a[:]));
    write(fileID,Float64.(app.pot3b[:]));
    write(fileID,Float64.(app.pot3c[:]));
    write(fileID,Float64.(app.pot4a[:]));
    write(fileID,Float64.(app.pot4b[:]));
    write(fileID,Float64.(app.rcutml^2));
    write(fileID,Float64.(app.rcut2a[:].^2));
    write(fileID,Float64.(app.rcut2b[:].^2));
    write(fileID,Float64.(app.rcut2c[:].^2));
    write(fileID,Float64.(app.rcut3a[:].^2));
    write(fileID,Float64.(app.rcut3b[:].^2));
    write(fileID,Float64.(app.rcut3c[:].^2));
    write(fileID,Float64.(app.rcut4a[:].^2));
    write(fileID,Float64.(app.rcut4b[:].^2));
    write(fileID,Float64.(app.rcutsqmax));
    write(fileID,Float64.(app.atom1b[:]));
    write(fileID,Float64.(app.atom2b[:]));
    write(fileID,Float64.(app.atom2c[:]));
    write(fileID,Float64.(app.atom3b[:]));
    write(fileID,Float64.(app.atom3c[:]));
    write(fileID,Float64.(app.atom4b[:]));

    close(fileID);

end