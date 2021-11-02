%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function writeregion(reg,filename)

tmp = [reg.triclinic];

nsize = zeros(4,1);
nsize(1) = length(tmp);
nsize(2) = length(reg.boxlo(:));
nsize(3) = length(reg.boxhi(:));
nsize(4) = length(reg.boxtilt(:));

endian = 'native';
fileID = fopen(filename,'w');
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,tmp(:),'double',endian);
fwrite(fileID,reg.boxlo(:),'double',endian);
fwrite(fileID,reg.boxhi(:),'double',endian);
fwrite(fileID,reg.boxtilt(:),'double',endian);

fclose(fileID);


