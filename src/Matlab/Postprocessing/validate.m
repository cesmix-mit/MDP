%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function valid = validate(app)

if ~isempty(app.validatelist)
    filename = app.sourcepath + "C++/Main/validation.bin";
    fileID = fopen(filename,'r');
    tmp = fread(fileID,'double');    
    d = app.dim;
    nc = tmp(1);
    valid.conf = zeros(1, nc);
    valid.inum = zeros(1, nc);
    valid.e = zeros(1, nc);
    m = 1;    
    for i = 1:nc
        c = tmp(m+1);
        n = tmp(m+2);
        valid.conf(i) = c;
        valid.inum(i) = n/d;
        valid.x{i} = reshape(tmp((m+3):(m+2+n)),[d n/d]);
        valid.f{i} = reshape(tmp((m+3+n):(m+2+2*n)),[d n/d]);
        valid.e(i) = tmp(m+3+2*n);
        m = m + 3 + 2*n;
    end
    fclose(fileID);
end

