%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function app = writeapp(app,filename)

%disp('Writing app into a binary file...'); 
% filename = "testapp.bin";
% fileID = fopen(filename,'r');
% tmp = fread(fileID,'double');
% fclose(fileID);

nsize = zeros(60,1);
nsize(1) = length(app.ndims(:));
nsize(2) = length(app.flag(:));  
nsize(3) = length(app.bcs(:)); 
nsize(4) = length(app.pbc(:)); 
nsize(5) = length(app.boxoffset(:)); 
nsize(6) = length(app.atomnumbers(:));
nsize(7) = length(app.atommasses(:));
nsize(8) = length(app.atomcharges(:));
nsize(9) = length(app.simparam(:)); 
nsize(10) = length(app.solparam(:)); 
nsize(11) = length(app.eta(:));
nsize(12) = length(app.kappa(:));
nsize(13) = length(app.muml(:)); 
nsize(14) = length(app.mu1a(:));
nsize(15) = length(app.mu1b(:));
nsize(16) = length(app.mu2a(:));
nsize(17) = length(app.mu2b(:));
nsize(18) = length(app.mu2c(:));
nsize(19) = length(app.mu3a(:));
nsize(20) = length(app.mu3b(:));
nsize(21) = length(app.mu3c(:));
nsize(22) = length(app.mu4a(:));
nsize(23) = length(app.mu4b(:));
nsize(24) = length(app.pot1a(:));
nsize(25) = length(app.pot1b(:));
nsize(26) = length(app.pot2a(:));
nsize(27) = length(app.pot2b(:));
nsize(28) = length(app.pot2c(:));
nsize(29) = length(app.pot3a(:));
nsize(30) = length(app.pot3b(:));
nsize(31) = length(app.pot3c(:));
nsize(32) = length(app.pot4a(:));
nsize(33) = length(app.pot4b(:));
nsize(34) = length(app.rcutml(:));
nsize(35) = length(app.rcut2a(:));
nsize(36) = length(app.rcut2b(:));
nsize(37) = length(app.rcut2c(:));
nsize(38) = length(app.rcut3a(:));
nsize(39) = length(app.rcut3b(:));
nsize(40) = length(app.rcut3c(:));
nsize(41) = length(app.rcut4a(:));
nsize(42) = length(app.rcut4b(:));
nsize(43) = length(app.rcutsqmax(:));
nsize(44) = length(app.atom1b(:));
nsize(45) = length(app.atom2b(:));
nsize(46) = length(app.atom2c(:));
nsize(47) = length(app.atom3b(:));
nsize(48) = length(app.atom3c(:));
nsize(49) = length(app.atom4b(:));
nsize(51) = length(app.traininglist(:));
nsize(52) = length(app.validatelist(:));

endian = 'native';
app.nsize = nsize;
fileID = fopen(filename,'w');
fwrite(fileID,length(app.nsize(:)),'double',endian);
fwrite(fileID,app.nsize(:),'double',endian);
fwrite(fileID,app.ndims(:),'double',endian);
fwrite(fileID,app.flag(:),'double',endian);
fwrite(fileID,app.bcs(:),'double',endian);
fwrite(fileID,app.pbc(:),'double',endian);
fwrite(fileID,app.boxoffset(:),'double',endian);
fwrite(fileID,app.atomnumbers(:),'double',endian);
fwrite(fileID,app.atommasses(:),'double',endian);
fwrite(fileID,app.atomcharges(:),'double',endian);
fwrite(fileID,app.simparam(:),'double',endian);
fwrite(fileID,app.solparam(:),'double',endian);
fwrite(fileID,app.eta(:),'double',endian);
fwrite(fileID,app.kappa(:),'double',endian);
fwrite(fileID,app.muml(:),'double',endian);
fwrite(fileID,app.mu1a(:),'double',endian);
fwrite(fileID,app.mu1b(:),'double',endian);
fwrite(fileID,app.mu2a(:),'double',endian);
fwrite(fileID,app.mu2b(:),'double',endian);
fwrite(fileID,app.mu2c(:),'double',endian);
fwrite(fileID,app.mu3a(:),'double',endian);
fwrite(fileID,app.mu3b(:),'double',endian);
fwrite(fileID,app.mu3c(:),'double',endian);
fwrite(fileID,app.mu4a(:),'double',endian);
fwrite(fileID,app.mu4b(:),'double',endian);
fwrite(fileID,app.pot1a(:),'double',endian);
fwrite(fileID,app.pot1b(:),'double',endian);
fwrite(fileID,app.pot2a(:),'double',endian);
fwrite(fileID,app.pot2b(:),'double',endian);
fwrite(fileID,app.pot2c(:),'double',endian);
fwrite(fileID,app.pot3a(:),'double',endian);
fwrite(fileID,app.pot3b(:),'double',endian);
fwrite(fileID,app.pot3c(:),'double',endian);
fwrite(fileID,app.pot4a(:),'double',endian);
fwrite(fileID,app.pot4b(:),'double',endian);
fwrite(fileID,app.rcutml(:).^2,'double',endian);
fwrite(fileID,app.rcut2a(:).^2,'double',endian);
fwrite(fileID,app.rcut2b(:).^2,'double',endian);
fwrite(fileID,app.rcut2c(:).^2,'double',endian);
fwrite(fileID,app.rcut3a(:).^2,'double',endian);
fwrite(fileID,app.rcut3b(:).^2,'double',endian);
fwrite(fileID,app.rcut3c(:).^2,'double',endian);
fwrite(fileID,app.rcut4a(:).^2,'double',endian);
fwrite(fileID,app.rcut4b(:).^2,'double',endian);
fwrite(fileID,app.rcutsqmax(:),'double',endian);
fwrite(fileID,app.atom1b(:),'double',endian);
fwrite(fileID,app.atom2b(:),'double',endian);
fwrite(fileID,app.atom2c(:),'double',endian);
fwrite(fileID,app.atom3b(:),'double',endian);
fwrite(fileID,app.atom3c(:),'double',endian);
fwrite(fileID,app.atom4b(:),'double',endian);
fwrite(fileID,app.traininglist(:),'double',endian);
fwrite(fileID,app.validatelist(:),'double',endian);

fclose(fileID);


