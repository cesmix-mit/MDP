
% n = 100;
% the = pi*rand(n,1);
% phi = -pi + 2*pi*rand(n,1);
% r = 0.1 + 0.95*rand(n,1);
% data = [the phi r];
% 
% fileID = fopen("datain.bin",'w');
% fwrite(fileID,data,'double','native');
% fclose(fileID);

% x0 = besselzero((1:25)'-1/2, 20, 2);
% fileID = fopen("besselzeros.bin",'w');
% fwrite(fileID,x0,'double','native');
% fclose(fileID);

n = 100;
dim = 3;
L = 3;
K = 5;

fileID = fopen("datain.bin",'r');
data = fread(fileID,'double');
fclose(fileID);
data = reshape(data,[n dim]);
the = data(:,1);
phi = data(:,2);
r = data(:,3);

[Ylmr, Ylmi] = sphericalharmonics(the,phi,L);
[x,y,z] = sphere2cart(the,phi,r);
[ar, ai]= shsum(x,y,z,L);
[ar, ai, arx, ary, arz, aix, aiy, aiz] = sssumderiv(x,y,z,K,L);
b = shbispectrum(ar, ai);
indl = uniquebispectrum(b);
[cg,indm,rowm] = cgcoefficients(indl);
[ar, ai]= sssum(x,y,z,K,L);
[ar{1} ar{2} ar{3} ar{4}]
[ai{1} ai{2} ai{3} ai{4}]
[arx{1} reshape(arx{2},n,[]) reshape(arx{3},n,[]) reshape(arx{4},n,[])]

p = sspower(ar, ai);
reshape(p,L+1,[])
[p,px,py,pz] = sspowerderiv(ar, ai, arx, ary, arz, aix, aiy, aiz);

b = ssbispectrum(ar, ai, cg, indl, indm, rowm);
reshape(p,[],K*(K+1)/2)
[b,bx,by,bz] = ssbispectrumderiv(ar, ai, arx, ary, arz, aix, aiy, aiz, cg, indl, indm, rowm);


%g = sphericalbessel2(r,L,K);

