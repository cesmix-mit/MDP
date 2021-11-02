%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [app, config] = fitting(app)

% create exec folder if it does not exist
bindir = "exec";
cd(char(app.sourcepath));
if exist(bindir, "file") == 0
    mkdir(char(bindir));    
end
cd(char(app.currentdir));

if lower(app.fitmethod) == "lsq"
    app.training = 1;
end

ndata = length(app.trainingfiles);
for i = 1:ndata        
    app.configfile = app.trainingfiles(i);    
    app.configweights = [app.eweight(i) app.fweight(i) app.vweight(i)];

    [app, config] = preprocessing(app);     
    if i==1
        gencode(app);                     
        compilecode(app);                     
    end        
    runcode(app);
    
    [Aei, Afi, Avi] = descriptors(app);
    
    Ae(:,:,i) = Aei'*Aei;
    Af(:,:,i) = Afi'*Afi;
    Av(:,:,i) = Avi'*Avi;
    
%    pause
%         
%     [Ai, bi] = leastsquarefit(app, config);
%     if i==1
%         A = Ai;
%         b = bi;
%     else
%         A = A + Ai;
%         b = b + bi;
%     end        
end    

