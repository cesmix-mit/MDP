%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [c, A, b, Ae, Af, Av, be, bf, bv, Aei, Afi, Avi, app, config] = leastsquarefit(app)

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
    app.nconfigs
    if i==1
        %gencode(app);                     
        compilecode(app);                     
    end        
    
    runcode(app);        
    [Aei, Afi, Avi] = descriptors(app, config);
    
    [size(Aei) size(Afi) size(Avi)]
    if (i==1)
        M = size(Aei,2);
        Ae = zeros(M,M,ndata);
        Af = zeros(M,M,ndata);
        Av = zeros(M,M,ndata);
        be = zeros(M,ndata);
        bf = zeros(M,ndata);
        bv = zeros(M,ndata);
    end
    
    Ae(:,:,i) = Aei'*Aei;
    Af(:,:,i) = Afi'*Afi;
    Av(:,:,i) = Avi'*Avi;
    if ~isempty(config.e)
        be(:,i) = Aei'*config.e(:);
    end
    if ~isempty(config.f)
        bf(:,i) = Afi'*config.f(:);
    end
    if ~isempty(config.stress)
        bv(:,i) = Avi'*reshape(config.stress([1 2 3 5 6 9],:),[],1);
    end
end    

M = size(Ae,1);
A = zeros(M,M);
b = zeros(M,1);
if app.dftdata == 1
    for i = 1:ndata
        A = A + app.eweight(i)*Ae(:,:,i);
        b = b + app.eweight(i)*be(:,i);
    end
elseif app.dftdata == 2
    for i = 1:ndata
        A = A + app.fweight(i)*Af(:,:,i);
        b = b + app.fweight(i)*bf(:,i);
    end    
elseif app.dftdata == 3
    for i = 1:ndata
        A = A + app.eweight(i)*Ae(:,:,i) + app.fweight(i)*Af(:,:,i);
        b = b + app.eweight(i)*be(:,i) + app.fweight(i)*bf(:,i);
    end        
elseif app.dftdata == 4
    for i = 1:ndata
        A = A + app.eweight(i)*Ae(:,:,i) + app.fweight(i)*Af(:,:,i) + app.vweight(i)*Av(:,:,i);
        b = b + app.eweight(i)*be(:,i) + app.fweight(i)*bf(:,i) + app.vweight(i)*bv(:,i);
    end            
end

c = pcg(A , b, 1e-14, 10000);

end



