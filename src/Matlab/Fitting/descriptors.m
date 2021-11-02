%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [Ae, Af, Av] = descriptors(app, config)

bindir = "exec";
execpath = app.sourcepath + bindir;

filename = execpath + "/out_p0_ai.bin";
fileID = fopen(filename,'r'); Ae = fread(fileID,'double'); fclose(fileID);

filename = execpath + "/out_p0_ad.bin";
fileID = fopen(filename,'r'); tm = fread(fileID,'double'); fclose(fileID);

filename = execpath + "/out_p0_av.bin";
fileID = fopen(filename,'r'); Av = fread(fileID,'double'); fclose(fileID);

nconfigs =  config.nconfigs;
Ae = reshape(Ae, [], nconfigs);
Ae = permute(Ae, [2 1]);

N = size(Ae,2);
Av = reshape(Av, 6, N, nconfigs);
Av = reshape(permute(Av, [1 3 2]),[], N);

m1 = 0;
n1 = 0;
Af = zeros(3*sum(config.natom),N);
for i = 1:nconfigs
    natom = config.natom(i);
    m = 3*natom;
    n = 3*natom*N;    
    m2 = m1 + m;
    n2 = n1 + n;
    Af((m1+1):m2,:) = reshape(tm((n1+1):n2),3*natom,N);
    n1 = n1 + n;
    m1 = m1 + m;
end

filename = execpath + "/out_p0_bi.bin";
fileID = fopen(filename,'r'); Be = fread(fileID,'double'); fclose(fileID);

filename = execpath + "/out_p0_bd.bin";
fileID = fopen(filename,'r'); tm = fread(fileID,'double'); fclose(fileID);

filename = execpath + "/out_p0_bv.bin";
fileID = fopen(filename,'r'); Bv = fread(fileID,'double'); fclose(fileID);

Be = reshape(Be, [], nconfigs);
Be = permute(Be, [2 1]);

M = size(Be,2);
Bv = reshape(Bv, 6, M, nconfigs);
Bv = reshape(permute(Bv, [1 3 2]),[], M);

m1 = 0;
n1 = 0;
Bf = zeros(3*sum(config.natom),M);
for i = 1:nconfigs
    natom = config.natom(i);
    m = 3*natom;
    n = 3*natom*M;    
    m2 = m1 + m;
    n2 = n1 + n;
    Bf((m1+1):m2,:) = reshape(tm((n1+1):n2),3*natom,M);
    n1 = n1 + n;
    m1 = m1 + m;
end

[size(Ae) size(Be)]

Ae = [Ae Be];
Af = [Af Bf];
Av = [Av Bv];



