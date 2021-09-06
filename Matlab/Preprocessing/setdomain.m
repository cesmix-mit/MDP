%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function dom = setdomain(boxhi, boxtilt, pbc, bcs)

if nargin<1
    error("Invalid input arguments");
end
if length(boxhi) ~= 3
    error("length of boxhi must be 3");
end    

dom.boxlo = [0.0 0.0 0.0];
dom.boxhi = boxhi;

if nargin>=2 && ~isempty(boxtilt)
    if length(boxtilt) ~= 3
        error("length of boxtilt must be 3");
    end        
    dom.boxtilt = boxtilt;
    dom.triclinic = 1;
else
    dom.boxtilt = [0 0 0];
    dom.triclinic = 0;
end

if nargin>=3 && ~isempty(pbc)
    if length(pbc) ~= 3
        error("length of pbc must be 3");
    end        
    dom.pbc = pbc;
else
    dom.pbc = [1 1 1];
end

if nargin>=4 && ~isempty(bcs)
    if length(bcs) ~= 3
        error("length of bcs must be 6");
    end        
    dom.bcs = bcs;
else
    dom.bcs = [1 1 1 1 1 1];
end

end


