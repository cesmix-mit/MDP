function jlist = createjlist(xi, ilist, r, inum)

n = inum;
dim = size(xi,1);
rsq = r*r;
jlist = cell(n,1);
for i = 1:n
    dx = xi(:,i) - xi;
    ri = dx(1,:).^2;
    for j=2:dim
        ri = ri + dx(j,:).^2;
    end
    %jlist{i} = ilist((ri <= rsq));     
    jlist{i} = find(ri <= rsq);
    jlist{i} = setdiff(jlist{i},i);
end


