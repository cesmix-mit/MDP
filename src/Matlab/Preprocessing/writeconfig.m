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
fclose(fileID);

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
