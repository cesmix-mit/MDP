%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function writelattice(lat,filename)

tmp = [lat.style size(lat.basis,1) lat.spaceflag lat.scale];
basis = lat.basis';

nsize = zeros(11,1);
nsize(1) = length(tmp);
nsize(2) = length(lat.origin(:));
nsize(3) = length(lat.orientx(:));
nsize(4) = length(lat.orienty(:));
nsize(5) = length(lat.orientz(:));
nsize(6) = length(lat.spacing(:));
nsize(7) = length(lat.a1(:));
nsize(8) = length(lat.a2(:));
nsize(9) = length(lat.a3(:));
nsize(10) = length(basis(:));
nsize(11) = length(lat.type(:));

endian = 'native';
fileID = fopen(filename,'w');
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,tmp(:),'double',endian);
fwrite(fileID,lat.origin(:),'double',endian);
fwrite(fileID,lat.orientx(:),'double',endian);
fwrite(fileID,lat.orienty(:),'double',endian);
fwrite(fileID,lat.orientz(:),'double',endian);
fwrite(fileID,lat.spacing(:),'double',endian);
fwrite(fileID,lat.a1(:),'double',endian);
fwrite(fileID,lat.a2(:),'double',endian);
fwrite(fileID,lat.a3(:),'double',endian);
fwrite(fileID,basis(:),'double',endian);
fwrite(fileID,lat.type(:),'double',endian);

fclose(fileID);


