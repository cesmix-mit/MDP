%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function lat = setlattice(style, scale, orientx, orienty, orientz, spacing, type, a1, a2, a3, basis)

%  NONE = 0, SC = 1, BCC = 2, FCC = 3, HCP = 4, DIAMOND = 5, SQ = 6, SQ2 = 7, HEX = 8, CUSTOM = 9;

if nargin<2
    error("Invalid input arguments");
end

lat.scale = scale;
lat.spaceflag = 0;

if (strcmpi(style,"none")) 
    lat.style = 0;
elseif (strcmpi(style,"sc")) 
    lat.style = 1;    
elseif (strcmpi(style,"bcc")) 
    lat.style = 2;
elseif (strcmpi(style,"fcc")) 
    lat.style = 3;
elseif (strcmpi(style,"hcp")) 
    lat.style = 4;
elseif (strcmpi(style,"diamond")) 
    lat.style = 5;
elseif (strcmpi(style,"sq")) 
    lat.style = 6;
elseif (strcmpi(style,"sq2")) 
    lat.style = 7;
elseif (strcmpi(style,"hex")) 
    lat.style = 8;        
elseif (strcmpi(style,"custom")) 
    lat.style = 9;        
else
    error("Invalid lattice style");
end

%set defaults 
lat.origin = [0 0 0];
if nargin >=3 && ~isempty(orientx)
    if length(orientx) ~= 3
        error("length of orientx must be 3");
    end    
    lat.orientx = orientx;
else
    lat.orientx = [1 0 0];
end
if nargin >=4 && ~isempty(orienty)
    if length(orienty) ~= 3
        error("length of orienty must be 3");
    end        
    lat.orienty = orienty;
else
    lat.orienty = [0 1 0];
end
if nargin >=5 && ~isempty(orientz)
    if length(orientz) ~= 3
        error("length of orientz must be 3");
    end        
    lat.orientz = orientz;
else
    lat.orientz = [0 0 1];
end
if nargin >=6 && ~isempty(spacing)
    if length(spacing) ~= 3
        error("length of spacing must be 3");
    end        
    lat.spaceflag = 1;
    lat.spacing = spacing;
else
    lat.spacing = [1 1 1];    
end

if (strcmpi(style,"custom"))
    if nargin<11
        error("Invalid input arguments");
    end
    if length(a1) ~= 3
        error("length of a1 must be 3");
    end        
    if length(a2) ~= 3
        error("length of a2 must be 3");
    end        
    if length(a3) ~= 3
        error("length of a3 must be 3");
    end        
    if size(basis,2) ~= 3
        error("the 2nd dimension of basis must be 3");
    end        
    lat.a1 = a1;
    lat.a2 = a2;
    lat.a3 = a3;
    lat.basis = basis;    
else
    lat.a1 = [1 0 0];
    lat.a2 = [0 1 0];
    lat.a3 = [0 0 1];
    
    if (strcmpi(style,"hex")) 
        lat.a2(2) = sqrt(3.0);
    end
    if (strcmpi(style,"hcp")) 
        lat.a2(2) = sqrt(3.0);
        lat.a3(3) = sqrt(8.0/3.0);
    end
    
    if (strcmpi(style,"sc")) 
        lat.basis = [0 0 0];        
    elseif (strcmpi(style,"bcc")) 
        lat.basis = [0 0 0; 0.5 0.5 0.5];
    elseif (strcmpi(style,"fcc")) 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5];
    elseif (strcmpi(style,"hcp")) 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 5.0/6.0 0.5; 0.0 1.0/3.0 0.5];
    elseif (strcmpi(style,"diamond")) 
        lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5; 0.25 0.25 0.25; 0.25 0.75 0.75; 0.75 0.25 0.75; 0.75 0.75 0.25];        
    elseif (strcmpi(style,"sq")) 
        lat.basis = [0 0 0];        
    elseif (strcmpi(style,"sq2")) 
        lat.basis = [0 0 0; 0.5 0.5 0.0];
    elseif (strcmpi(style,"hex")) 
        lat.basis = [0 0 0; 0.5 0.5 0.0];    
    end
end

nbasis = size(lat.basis,1);
if nargin >=7 && ~isempty(type)
    if length(type) ~= nbasis
        error("length of type must be equal to # basis atoms");
    end            
    lat.type = type;
else
    lat.type = ones(1,nbasis);    
end


