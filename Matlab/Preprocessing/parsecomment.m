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
%         value = values{i};        
%         indc = strfind(value, ':');        
%         nl = [1 indc+1];
%         nr = [indc-1 length(value)];                
%         nc = length(nr);
%         c = cell(nc,1);
%         for j = 1:nc            
%             c{j} = value(nl(j):nr(j));
%         end 
        indp = i;
        break;        
    end
end

end

