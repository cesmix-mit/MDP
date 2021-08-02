%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function config = readLAMMPS(tmp)
    
config.dim = 3;
config.ncq = 0;
config.ncv = 0;
config.nce = 1;
config.ncf = 3;
config.ncx = 3;
config.natom = [];
config.t = [];
config.x = [];
config.f = [];
config.e = [];
config.a = [];
config.b = [];
config.c = [];
config.q = [];
config.v = [];

i = 1;
k = 1;
while (k < length(tmp))
    n = tmp(k); % number of atoms
    config.natom = [config.natom n];
    b = reshape(tmp((k+1):(k + 7*n)), [7 n]);
    config.t = [config.t b(1,:)]; 
    config.x = [config.x b(2:4,:)]; 
    config.f = [config.f b(5:7,:)]; 
    b = reshape(tmp((k + 1 + 7*n):(k + 4 + 7*n)), [1 4]);
    config.e = [config.e b(1)];
    config.a = [config.a [b(2) 0 0]'];
    config.b = [config.b [0 b(3) 0]'];
    config.c = [config.c [0 0 b(4)]'];
    i = i + 1;
    k = k + 7*n + 5;
end

end
