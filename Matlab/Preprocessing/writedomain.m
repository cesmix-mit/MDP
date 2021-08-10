%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function writedomain(dom,filename)

tmp = [dom.triclinic];

nsize = zeros(6,1);
nsize(1) = length(tmp);
nsize(2) = length(dom.boxlo(:));
nsize(3) = length(dom.boxhi(:));
nsize(4) = length(dom.boxtilt(:));
nsize(5) = length(dom.pbc(:));
nsize(6) = length(dom.bcs(:));

endian = 'native';
fileID = fopen(filename,'w');
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize(:),'double',endian);
fwrite(fileID,tmp(:),'double',endian);
fwrite(fileID,dom.boxlo(:),'double',endian);
fwrite(fileID,dom.boxhi(:),'double',endian);
fwrite(fileID,dom.boxtilt(:),'double',endian);
fwrite(fileID,dom.pbc(:),'double',endian);
fwrite(fileID,dom.bcs(:),'double',endian);

fclose(fileID);


