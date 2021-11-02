%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function writeconfig(config, filename)

% disp('Writing configurations into a binary file...'); 

nconfigs = config.nconfigs; % number of configurations
dim = config.dim; % physical dimension
ncv = config.ncv; % number of compoments of v
nce = config.nce; % number of compoments of e
ncf = config.ncf; % number of compoments of f
ncq = config.ncq; % number of compoments of q

endian = 'native';
tmp = [nconfigs dim ncq ncv nce ncf];

nsize = zeros(20,1);
nsize(1) = length(tmp(:));
nsize(2) = length(config.natom(:));  
nsize(3) = length(config.a(:)); 
nsize(4) = length(config.b(:)); 
nsize(5) = length(config.c(:)); 
nsize(6) = length(config.e(:));
nsize(7) = length(config.t(:));
nsize(8) = length(config.x(:));
nsize(9) = length(config.q(:)); 
nsize(10) = length(config.v(:)); 
nsize(11) = length(config.f(:));    
nsize(12) = length(config.we(:)); 
nsize(13) = length(config.wf(:));    
nsize(14) = length(config.ws(:)); 
nsize(15) = length(config.stress(:)); 
nsize(16) = length(config.pbc(:)); 
nsize(17) = length(config.move(:)); 
nsize(18) = length(config.group(:)); 
%nsize(19) = length(config.tags(:)); 

fileID = fopen(filename,'w');
fwrite(fileID,length(nsize(:)),'double',endian);
fwrite(fileID,nsize,'double',endian);
fwrite(fileID,tmp,'double',endian);
fwrite(fileID,config.natom,'double',endian);
fwrite(fileID,config.a,'double',endian);    
fwrite(fileID,config.b,'double',endian);    
fwrite(fileID,config.c,'double',endian); 
fwrite(fileID,config.e,'double',endian);    
fwrite(fileID,config.t,'double',endian);    
fwrite(fileID,config.x,'double',endian);    
fwrite(fileID,config.q,'double',endian);    
fwrite(fileID,config.v,'double',endian);    
fwrite(fileID,config.f,'double',endian);
fwrite(fileID,config.we,'double',endian);
fwrite(fileID,config.wf,'double',endian);
fwrite(fileID,config.ws,'double',endian);
fwrite(fileID,config.stress,'double',endian);
fwrite(fileID,config.pbc,'double',endian);
fwrite(fileID,config.move,'double',endian);
fwrite(fileID,config.group,'double',endian);
%fwrite(fileID,config.tags,'double',endian);

fclose(fileID);

% config.natom = ones(1, nconfigs);   % number of atoms per each configuration
% config.lattice = zeros(ncl, nconfigs); % stresses for all configurations
% config.pbc = zeros(ncp, nconfigs); % periodic boundary conditions
% config.e = ones(nce, nconfigs);  % potential energies for all configurations
% config.stress = zeros(ncs, nconfigs); % stresses for all configurations
% 
% config.natomall = sum(config.natom);% number of atoms for all configurations
% config.we = ones(1, nconfigs);      % energy weight per each configuration
% config.wf = ones(1, nconfigs);      % force weight per each configuration
% config.ws = ones(1, nconfigs);      % stress weight per each configuration
% 
% % simulation box for each configuration
% config.a = zeros(3, nconfigs); % the 1st principal vector of the simulation box
% config.b = zeros(3, nconfigs); % the 2nd principal vector of the simulation box
% config.c = zeros(3, nconfigs); % the 3rd principal vector of the simulation box
% 
% config.Z = zeros(ncz, config.natomall);    % atom numbers for all configurations
% config.mass = zeros(ncm, config.natomall); % atom masses for all configurations
% config.move = zeros(nco, config.natomall);   % atom move masks for all configurations
% % config.eatom = zeros(1, config.natomall); % atom energies for all configurations
% % config.vatom = zeros(6, config.natomall); % atom virial for all configurations
% config.tags = zeros(nci, config.natomall); % atom tags for all configurations
% config.t = zeros(nct, config.natomall);   % atom types for all configurations
% config.g = zeros(ncg, config.natomall);   % atom groups for all configurations
% config.x = zeros(ncx, config.natomall); % atom positions for all configurations
% config.q = zeros(ncq, config.natomall); % atom charges for all configurations
% config.v = zeros(ncv, config.natomall); % atom velocities for all configurations
% config.f = zeros(ncf, config.natomall); % atom forces for all configurations
% 

% fileID = fopen(filename,'w');
% if mode==0 % binary    
%     fwrite(fileID,tmp,'double',endian);
%     fwrite(fileID,[config.natom(:) config.a' config.b' config.c' config.e(:) config.we(:) config.wf(:)],'double',endian);
%     fwrite(fileID,[config.t(:) config.x' config.q' config.v' config.f'],'double',endian);    
% elseif mode==1 % text    
%     fprintf(fileID,'%d %d %d %d %d %d\n',tmp);
%     
%     nc = 1 + dim*dim + nce;
%     mystr = "%d";
%     for d = 1:(nc-1)
%         mystr = mystr + " %.16f";
%     end
%     mystr = mystr + " \n";    
%     fprintf(fileID,char(mystr),[config.natom(:) config.a' config.b' config.c' config.e(:)]);        
%     
%     nc = 1 + ncx + ncq + ncv + ncf;
%     mystr = "%d";
%     for d = 1:(nc-1)
%         mystr = mystr + " %.16f";
%     end
%     mystr = mystr + " \n";    
%     fprintf(fileID,char(mystr),[config.t(:) config.x' config.q' config.v' config.f']);    
% else
%     nsize = zeros(20,1);
%     nsize(1) = length(tmp(:));
%     nsize(2) = length(config.natom(:));  
%     nsize(3) = length(config.a(:)); 
%     nsize(4) = length(config.b(:)); 
%     nsize(5) = length(config.c(:)); 
%     nsize(6) = length(config.e(:));
%     nsize(7) = length(config.t(:));
%     nsize(8) = length(config.x(:));
%     nsize(9) = length(config.q(:)); 
%     nsize(10) = length(config.v(:)); 
%     nsize(11) = length(config.f(:));    
%     
%     fwrite(fileID,length(nsize(:)),'double',endian);
%     fwrite(fileID,nsize,'double',endian);
%     fwrite(fileID,tmp,'double',endian);
%     fwrite(fileID,config.natom,'double',endian);
%     fwrite(fileID,config.a,'double',endian);    
%     fwrite(fileID,config.b,'double',endian);    
%     fwrite(fileID,config.c,'double',endian); 
%     fwrite(fileID,config.e,'double',endian);    
%     fwrite(fileID,config.t,'double',endian);    
%     fwrite(fileID,config.x,'double',endian);    
%     fwrite(fileID,config.q,'double',endian);    
%     fwrite(fileID,config.v,'double',endian);    
%     fwrite(fileID,config.f,'double',endian);
%     fwrite(fileID,config.we,'double',endian);
%     fwrite(fileID,config.wf,'double',endian);
% end
% fclose(fileID);
% 
% 
