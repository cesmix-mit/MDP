
for L = 1:8
    K = L+1;
    
    n = 100;
%     x=rand(n,1);
%     y=rand(n,1);
%     z=rand(n,1);
%     r = sqrt(x.^2 + y.^2 + z.^2);
    the = pi*rand(n,1);
    phi = -pi + 2*pi*rand(n,1);
    r = 0.1 + 0.95*rand(n,1);
    [x,y,z] = sphere2cart(the,phi,r);
    
    euler = [2*pi*rand pi*rand 2*pi*rand];
    R = euler2rotm(euler);    
    R = squeeze(R);
    [X,Y,Z] = rotc(R, x, y, z);
        
    [~, indl] = shspectrum(x,y,z,L);
    [cg,indm,rowm] = cgcoefficients(indl);
    [ar, ai] = sssum(x,y,z,K,L);
    p = sspower(ar, ai);
    b = ssbispectrum(ar, ai, cg, indl, indm, rowm);
            
    [~, indl] = shspectrum(X,Y,Z,L);
    [cg,indm,rowm] = cgcoefficients(indl);
    [ar, ai] = sssum(X,Y,Z,K,L);
    p2 = sspower(ar, ai);
    b2 = ssbispectrum(ar, ai, cg, indl, indm, rowm);

    disp(" ");
    fprintf("Euler's angles: [%g %g %g]\n", euler);
    fprintf("Minium radius: %g\n", min(abs(r(:))));
    fprintf("Maximum power spectrum error: %g\n", max(abs((p(:)-p2(:))./p(:))));
    fprintf("Maximum bispectrum error: %g\n", max(abs((b(:)-b2(:))./b(:))));        
    fprintf("Power spectrum range: [%g %g]\n", [min(abs(p(:))) max(abs(p(:)))]);      
    fprintf("Bipectrum range: [%g %g]\n", [min(abs(b(:))) max(abs(b(:)))]);             
end

K = 8; L = 8;
n = 100;
dx = 1e-8;
the = pi*rand(n,1);
phi = -pi + 2*pi*rand(n,1);
r = 0.05 + 0.99*rand(n,1);
[x,y,z] = sphere2cart(the,phi,r);

[Sr, Si, Srx, Sry, Srz, Six, Siy, Siz] = ssharmonicsderiv(x,y,z,K,L);

[Sr2, Si2] = sphericalshellbasis(x+dx,y,z,K,L);
[Sr1, Si1] = sphericalshellbasis(x-dx,y,z,K,L);
for l = 1:(L+1)
    e = (Sr2{l}-Sr1{l})/(2*dx) - Srx{l};
    fprintf("Maximum derivative error Srx: %g\n", max(abs(e(:)))/max(abs(Srx{l}(:))));    
    e = (Si2{l}-Si1{l})/(2*dx) - Six{l};
    fprintf("Maximum derivative error Six: %g\n", max(abs(e(:)))/max(abs(Six{l}(:))));    
end

[Sr2, Si2] = sphericalshellbasis(x,y+dx,z,K,L);
[Sr1, Si1] = sphericalshellbasis(x,y-dx,z,K,L);
for l = 1:(L+1)
    e = (Sr2{l}-Sr1{l})/(2*dx) - Sry{l};
    fprintf("Maximum derivative error Sry: %g\n", max(abs(e(:)))/max(abs(Sry{l}(:))));    
    e = (Si2{l}-Si1{l})/(2*dx) - Siy{l};
    fprintf("Maximum derivative error Siy: %g\n", max(abs(e(:)))/max(abs(Siy{l}(:))));    
end

[Sr2, Si2] = sphericalshellbasis(x,y,z+dx,K,L);
[Sr1, Si1] = sphericalshellbasis(x,y,z-dx,K,L);
for l = 1:(L+1)
    e = (Sr2{l}-Sr1{l})/(2*dx) - Srz{l};
    fprintf("Maximum derivative error Srz: %g\n", max(abs(e(:)))/max(abs(Srz{l}(:))));    
    e = (Si2{l}-Si1{l})/(2*dx) - Siz{l};
    fprintf("Maximum derivative error Siz: %g\n", max(abs(e(:)))/max(abs(Siz{l}(:))));    
end

[ar, ai, arx, ary, arz, aix, aiy, aiz]= sssumderiv(x,y,z,K,L);
darx = arx;
dary = ary;
darz = arz;
daix = aix;
daiy = aiy;
daiz = aiz;
for i = 1:n
    x2 = x; x2(i) = x2(i) + dx;
    [ar2, ai2] = sssum(x2,y,z,K,L);
    x1 = x; x1(i) = x1(i) - dx;
    [ar1, ai1] = sssum(x1,y,z,K,L);    
    for l = 1:(L+1)
        darx{l}(i,:,:) = (ar2{l}-ar1{l})/(2*dx);
        daix{l}(i,:,:) = (ai2{l}-ai1{l})/(2*dx);
    end
    
    y2 = y; y2(i) = y2(i) + dx;
    [ar2, ai2] = sssum(x,y2,z,K,L);
    y1 = y; y1(i) = y1(i) - dx;
    [ar1, ai1] = sssum(x,y1,z,K,L);    
    for l = 1:(L+1)
        dary{l}(i,:,:) = (ar2{l}-ar1{l})/(2*dx);
        daiy{l}(i,:,:) = (ai2{l}-ai1{l})/(2*dx);
    end
    
    z2 = z; z2(i) = z2(i) + dx;
    [ar2, ai2] = sssum(x,y,z2,K,L);
    z1 = z; z1(i) = z1(i) - dx;
    [ar1, ai1] = sssum(x,y,z1,K,L);    
    for l = 1:(L+1)
        darz{l}(i,:,:) = (ar2{l}-ar1{l})/(2*dx);
        daiz{l}(i,:,:) = (ai2{l}-ai1{l})/(2*dx);
    end
end

for l = 1:(L+1)
    e = darx{l} - arx{l};
    fprintf("Maximum derivative error arx: %g\n", max(abs(e(:)))/max(abs(arx{l}(:))));    
    e = daix{l} - aix{l};
    fprintf("Maximum derivative error aix: %g\n", max(abs(e(:)))/max(abs(aix{l}(:))));        
end

for l = 1:(L+1)
    e = dary{l} - ary{l};
    fprintf("Maximum derivative error ary: %g\n", max(abs(e(:)))/max(abs(ary{l}(:))));    
    e = daiy{l} - aiy{l};
    fprintf("Maximum derivative error aiy: %g\n", max(abs(e(:)))/max(abs(aiy{l}(:))));        
end

for l = 1:(L+1)
    e = darz{l} - arz{l};
    fprintf("Maximum derivative error arz: %g\n", max(abs(e(:)))/max(abs(arz{l}(:))));    
    e = daiz{l} - aiz{l};
    fprintf("Maximum derivative error aiz: %g\n", max(abs(e(:)))/max(abs(aiz{l}(:))));        
end


[p,px,py,pz] = sspowerderiv(ar, ai, arx, ary, arz, aix, aiy, aiz);
dpx = 0*px;
dpy = 0*py;
dpz = 0*pz;
for i = 1:n
    x2 = x; x2(i) = x2(i) + dx;    
    [ar2, ai2] = sssum(x2,y,z,K,L);
    p2 = sspower(ar2, ai2);
    x1 = x; x1(i) = x1(i) - dx;
    [ar1, ai1] = sssum(x1,y,z,K,L);    
    p1 = sspower(ar1, ai1);
    dpx(i,:) = (p2-p1)/(2*dx);          
    
    y2 = y; y2(i) = y2(i) + dx;
    [ar2, ai2] = sssum(x,y2,z,K,L);
    p2 = sspower(ar2, ai2);
    y1 = y; y1(i) = y1(i) - dx;
    [ar1, ai1] = sssum(x,y1,z,K,L);    
    p1 = sspower(ar1, ai1);        
    dpy(i,:) = (p2-p1)/(2*dx);                
    
    z2 = z; z2(i) = z2(i) + dx;
    [ar2, ai2] = sssum(x,y,z2,K,L);
    p2 = sspower(ar2, ai2);
    z1 = z; z1(i) = z1(i) - dx;
    [ar1, ai1] = sssum(x,y,z1,K,L);    
    p1 = sspower(ar1, ai1);        
    dpz(i,:) = (p2-p1)/(2*dx);     
end

e = dpx-px;
fprintf("Maximum derivative error px: %g\n", max(abs(e(:)))/max(abs(px(:))));    
e = dpy-py;
fprintf("Maximum derivative error py: %g\n", max(abs(e(:)))/max(abs(py(:))));    
e = dpz-pz;
fprintf("Maximum derivative error pz: %g\n", max(abs(e(:)))/max(abs(pz(:))));    


[~, indl] = shspectrum(x,y,z,L);
[cg,indm,rowm] = cgcoefficients(indl);
[b,bx,by,bz] = ssbispectrumderiv(ar, ai, arx, ary, arz, aix, aiy, aiz, cg, indl, indm, rowm);
dbx = 0*bx;
dby = 0*by;
dbz = 0*bz;
for i = 1:n
    x2 = x; x2(i) = x2(i) + dx;    
    [ar2, ai2] = sssum(x2,y,z,K,L);
    b2 = ssbispectrum(ar2, ai2, cg, indl, indm, rowm); 
    x1 = x; x1(i) = x1(i) - dx;
    [ar1, ai1] = sssum(x1,y,z,K,L);    
    b1 = ssbispectrum(ar1, ai1, cg, indl, indm, rowm); 
    dbx(i,:) = (b2-b1)/(2*dx);          
    
    y2 = y; y2(i) = y2(i) + dx;
    [ar2, ai2] = sssum(x,y2,z,K,L);
    b2 = ssbispectrum(ar2, ai2, cg, indl, indm, rowm); 
    y1 = y; y1(i) = y1(i) - dx;
    [ar1, ai1] = sssum(x,y1,z,K,L);    
    b1 = ssbispectrum(ar1, ai1, cg, indl, indm, rowm); 
    dby(i,:) = (b2-b1)/(2*dx);                
    
    z2 = z; z2(i) = z2(i) + dx;
    [ar2, ai2] = sssum(x,y,z2,K,L);
    b2 = ssbispectrum(ar2, ai2, cg, indl, indm, rowm); 
    z1 = z; z1(i) = z1(i) - dx;
    [ar1, ai1] = sssum(x,y,z1,K,L);    
    b1 = ssbispectrum(ar1, ai1, cg, indl, indm, rowm); 
    dbz(i,:) = (b2-b1)/(2*dx);     
end

e = dbx-bx;
fprintf("Maximum derivative error bx: %g\n", max(abs(e(:)))/max(abs(bx(:))));    
e = dby-by;
fprintf("Maximum derivative error by: %g\n", max(abs(e(:)))/max(abs(by(:))));    
e = dbz-bz;
fprintf("Maximum derivative error bz: %g\n", max(abs(e(:)))/max(abs(bz(:))));    



