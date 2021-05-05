function config = readweight(app, config)

if app.weightmode ==0
    tmp = ones(2*config.nconfigs,1);
    tmp(1:config.nconfigs) = app.we;
    tmp((config.nconfigs+1):2*config.nconfigs) = app.wf;
elseif app.weightmode == 1 % binary    
    filename = app.weightfile;
    fileID = fopen(filename,'r');        
    tmp = fread(fileID,'double');    
    fclose(fileID);    
elseif app.weightmode == 2 % text    
    filename = app.weightfile;
    fileID = fopen(filename,'r');        
    tmp = fscanf(fileID, '%f');
    fclose(fileID);    
end

% determine weights
if (app.dftdata == 1)
    config.we = tmp(1:config.nconfigs);
    config.wf = [];
elseif (app.dftdata == 2)
    config.we = [];
    config.wf = tmp(1:config.nconfigs);
elseif (app.dftdata == 3)
    config.we = tmp(1:config.nconfigs);
    config.wf = tmp((config.nconfigs+1):2*config.nconfigs);
end

