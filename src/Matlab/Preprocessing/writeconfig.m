function writeconfig(app, config, filename, mode)

% nconfigs = app.nconfigs; % number of configurations
% dim = app.dim; % physical dimension
% ncx = app.ncx; % number of compoments of x
% ncv = app.ncv; % number of compoments of v
% ncf = app.ncf; % number of compoments of f
% ncq = app.ncq; % number of compoments of q

nconfigs = app.nconfigs; % number of configurations
dim = app.dim; % physical dimension
ncx = app.ncx; % number of compoments of x
ncv = app.ncv; % number of compoments of v
ncf = app.ncf; % number of compoments of f
ncq = app.ncq; % number of compoments of q

endian = 'native';
tmp = [nconfigs dim ncx ncq ncv ncf];

fileID = fopen(filename,'w');
if mode==0 % binary    
    fwrite(fileID,tmp,'double',endian);
    fwrite(fileID,[config.natom(:) config.a' config.b' config.c'],'double',endian);
    fwrite(fileID,[config.t(:) config.x' config.q' config.v' config.f'],'double',endian);    
%     fwrite(fileID,config.natom,'double',endian);
%     fwrite(fileID,config.a,'double',endian);    
%     fwrite(fileID,config.b,'double',endian);    
%     fwrite(fileID,config.c,'double',endian);    
%     fwrite(fileID,config.t,'double',endian);    
%     fwrite(fileID,config.x,'double',endian);    
%     fwrite(fileID,config.q,'double',endian);    
%     fwrite(fileID,config.v,'double',endian);    
%     fwrite(fileID,config.f,'double',endian);    
else % text    
    fprintf(fileID,'%d %d %d %d %d %d\n',tmp);
    
    nc = 1 + dim*dim;
    mystr = "%d";
    for d = 1:(nc-1)
        mystr = mystr + " %.16f";
    end
    mystr = mystr + " \n";    
    fprintf(fileID,char(mystr),[config.natom(:) config.a' config.b' config.c']');        
    
    nc = 1 + ncx + ncq + ncv + ncf;
    mystr = "%d";
    for d = 1:(nc-1)
        mystr = mystr + " %.16f";
    end
    mystr = mystr + " \n";    
    fprintf(fileID,char(mystr),[config.t(:) config.x' config.q' config.v' config.f']');    
end
fclose(fileID);


