%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function config = readEXTXYZ(filename, species)

fid = fopen(filename);

config.nconfigs = 0;
config.dim = 0;
config.ncx = 0;
config.ncv = 0;
config.ncq = 0;
config.nce = 0;
config.ncf = 0;
config.ncs = 0;
config.nci = 0;
config.nct = 0;
config.ncg = 0;
config.ncz = 0;
config.ncm = 0;
config.nco = 0;
config.ncl = 0;
config.ncp = 0;
config = initializeconfig(config);

natomall = 0;
nframe = 0;
while (1)     
    tline = fgetl(fid); 
    if ischar(tline) == 0
        break;
    end
    
    nframe = nframe + 1;
    natom = str2double(tline);
    
    config.natom(nframe) = natom;
    
    tline = fgetl(fid);  
    [keys, values, indp] = parsecomment(tline);
    for i = 1:length(keys)
        mystr = lower(string(keys{i}));
        if mystr == "lattice"
            config.lattice(:,nframe) = values{i};
            config.ncl = length(values{i});
        elseif mystr == "energy"
            config.e(nframe) = values{i};
            config.nce = 1;
        elseif mystr == "stress"
            config.stress(:,nframe) = values{i};
            config.ncs = length(values{i});
        elseif mystr == "pbc"
            config.pbc(:,nframe) = values{i};
            config.ncp = length(values{i});
        end
    end
            
    % handle properties
    prop = reshape(values{indp},3,[]);
    nfield = size(prop,2);        
    if nframe == 1
        fieldtype = prop(2,:);
        fieldlength = ones(nfield,1);        
        fmt = "";
        for i = 1:nfield            
            fieldlength(i) = str2double(prop(3,i));
            if lower(fieldtype(i)) == "s" || lower(fieldtype(i)) == "l"  % string        
                fmt = fmt + "%s";
            else                
                if fieldlength(i) == 3
                    fmt = fmt + "%f %f %f";
                elseif fieldlength(i) == 2
                    fmt = fmt + "%f %f";
                elseif fieldlength(i) == 1
                    fmt = fmt + "%f";                    
                end
            end
            if i < nfield
                fmt = fmt + " ";
            end
        end                
        fieldlength = cumsum(fieldlength);
    end        
        
    tm = textscan(fid, fmt, natom);       
    fgetl(fid);
    
    frame = cell(nfield,1);
    for j = 1:nfield        
        if j == 1
            k1 = 1;
            k2 = fieldlength(1);
        else
            k1 = fieldlength(j-1)+1;
            k2 = fieldlength(j);
        end
        if lower(fieldtype(j)) == "s" || lower(fieldtype(j)) == "l"  % string     
           frame{j} = string(tm{k1}); 
        else
           frame{j} = tm{k1}; 
           for k = (k1+1):k2
               frame{j} = [frame{j} tm{k}];
           end
        end        
    end
    
%     frame = cell(nfield,1);
%     for i = 1:natom
%         tline = split(fgetl(fid));          
%         for j = 1:nfield
%             if j == 1
%                 k1 = 1;
%                 k2 = fieldlength(1);
%             else
%                 k1 = fieldlength(j-1)+1;
%                 k2 = fieldlength(j);
%             end
%             if lower(fieldtype(j)) == "s" || lower(fieldtype(j)) == "l"  % string        
%                 for k = k1:k2                
%                     frame{j}(i,k+1-k1) = string(tline{k});                
%                 end
%             else
%                 for k = k1:k2                
%                     frame{j}(i,k+1-k1) = str2double(tline{k});                
%                 end
%             end
%         end        
%     end                   
    
    for j = 1:nfield
        b = ones(natom,1);
        if lower(fieldtype(j)) == "s" % string        
            for k=1:length(species)
               ind = frame{j} == species(k); 
               b(ind) = k;
            end
            frame{j} = b;
        end                
        b = zeros(natom,1);
        if lower(fieldtype(j)) == "l" % boolean   
            ind = frame{j} == 'T';
            b(ind) = 1;            
            frame{j} = b;
        end
        mystr = lower(prop(1,j));
        if mystr == "species"
            config.t((natomall+1):(natomall+natom)) = frame{j};
            config.nct = 1;
        elseif mystr == "pos"
            config.x(:,(natomall+1):(natomall+natom)) = frame{j}';
            config.ncx = size(config.x,1);
            config.dim = size(config.x,1);
        elseif mystr == "vel"
            config.v(:,(natomall+1):(natomall+natom)) = frame{j}';
            config.ncv = size(config.v,1);
        elseif mystr == "forces"
            config.f(:,(natomall+1):(natomall+natom)) = frame{j}';
            config.ncf = size(config.f,1);
        elseif mystr == "move_mask"
            config.move((natomall+1):(natomall+natom)) = frame{j};
            config.nco = 1;
        elseif mystr == "tags"
            config.tags((natomall+1):(natomall+natom)) = frame{j};
            config.nci = 1;
        elseif mystr == "groups"
            config.group((natomall+1):(natomall+natom)) = frame{j};    
            config.ncg = 1;
        elseif mystr == "charge"
            config.q((natomall+1):(natomall+natom)) = frame{j};
            config.ncq = 1;
        elseif mystr == "mass"
            config.mass((natomall+1):(natomall+natom)) = frame{j};
            config.ncm = 1;
        elseif mystr == "z"
            config.Z((natomall+1):(natomall+natom)) = frame{j};
            config.ncz = 1;
        end
    end               
    
    natomall = natomall + natom;
end

config.natomall = natomall; 
config.nconfigs = nframe;
fclose(fid);

end

function [keys, values, indp] = parsecomment(com)

dl = ['"' char("'") '[' '(' '{']; % left delimiters
dr = ['"' char("'") ']' ')' '}']; % right delimiters

com = char(com);         % convert to character vector
nchars = length(com);    % # of characters
inde = strfind(com, '=');% # of equal signs 
npairs = length(inde);   % # of (key, value) pairs

% find keys
keys = cell(npairs,1);
indk = zeros(1,npairs);
for i = 1:npairs
    nr = inde(i);
    keys{i} = [];
    for j = (nr-1):-1:1
        if com(j) == ' '
            break;
        else
            keys{i} = [keys{i} com(j)];
        end
    end
    keys{i} = keys{i}(end:-1:1);
    indk(i) = strfind(com, keys{i});
end

% find values
values = cell(npairs,1);
for i = 1:npairs    
    nl = inde(i)+1;
    if i<npairs        
        nr = indk(i+1)-2;
    else
        nr = nchars;
    end
    values{i} = com(nl:nr);    
    for j = 1:length(dl)        
        indl = strfind(values{i}, dl(j));        
        if isempty(indl) == 0            
            indr = strfind(values{i}, dr(j));            
            if isempty(indr) == 1
                error("Right delimiter is not found in value");
            else
                indl = indl(1);
                indr = indr(end);
                values{i} = values{i}((indl+1):(indr-1));
            end
        end
    end
    tm = textscan(values{i}, "%f");
    if isempty(tm{1}) == 0
        values{i} = tm{1};
    end
end

% handle pbc
for i = 1:npairs    
    if lower(string(keys{i})) == "pbc"        
        values{i} = strrep(values{i}, 'T', '1');
        values{i} = strrep(values{i}, 'F', '0');
        tm = textscan(values{i}, "%f");
        if isempty(tm{1}) == 0
            values{i} = tm{1};
        end    
    end
end

% handle properties
for i = 1:npairs    
    if lower(string(keys{i})) == "properties"        
        values{i} = split(string(values{i}),":");
        indp = i;
        break;        
    end
end

end

