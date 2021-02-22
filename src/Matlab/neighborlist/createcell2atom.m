function [cellcounts, cell2atom] = createcell2atom(clist, cellnum)


natom = length(clist); % number of atoms
ncell = prod(cellnum); % number of cells

% sort atom-to-cell list
[cell, cell2atom] = sort(clist);

% determine number of atoms for every cell
cellcounts = zeros(1,ncell+1);
for i = 0:(natom-1)
    i1 = i + 1;
    cid = cell(i1); % cell id of atom i
    
    % fill with zeros   
    if (i==0)
        for j = 1:cid
            cellcounts(j) = 0;
        end
    end
    
    % fill with natom
    if (i==natom-1)
        for j = (cid+1):(ncell+1)
            cellcounts(j) = natom;
        end
    end
    
    if ((i>0) && (i<natom))
        lid = cell(i1-1); % cell id of atom i-1
        if (lid ~= cid)
            for j = (lid+1):cid
                cellcounts(j) = i;
            end 
        end
    end    
end

    
