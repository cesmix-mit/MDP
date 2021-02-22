function neighlist = neighborlist(xi, atomtype, jlist, ellipsoid, inum)

neighlist = cell(inum,1);

for i = 1:inum  % loop over each i atom in the simulation domain    
    neighlist{i} = [];
    t = atomtype(i); % element of atom i 
    A = ellipsoid(:,:,t); % ellipsoid for element t        
    xt = xi(:,i);   % position of atom i
    m = length(jlist{i}); % number of atoms around i 
    for k = 1:m
        j = jlist{i}(k);  % atom j            
        % distance between atom i and atom j 
        rij = (xi(:,j)-xt)'*(A*(xi(:,j)-xt));
        if rij<=1 % atom j is 
            % add atom j into the list
            neighlist{i} = [neighlist{i} j];
        end            
    end        
end

