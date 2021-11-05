%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function config = readJSONfromFolder(pathdir, species)

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

cdir = pwd;
cd(pathdir);
jsonfiles = dir('**/*.json');
n = length(jsonfiles);

nframe = 0;
natomall = 0;
for i = 1:n
    
    filename = jsonfiles(i).name;

    % read JSON file to string
    fid = fopen(filename); 
    raw = fread(fid,inf); 
    mystr = char(raw'); 
    fclose(fid); 

    % remove comment line
    i1 = strfind(mystr,"#");
    i2 = strfind(mystr,"{");
    mystr(i1:i2(1)-1) = [];

    % convert JSON string to struct
    Dataset = jsondecode(mystr);
    Dataset = Dataset.Dataset.Data;

    for j = 1:length(Dataset)
        nframe = nframe + 1;
        Data = Dataset(j);

        if isfield(Data,'NumAtoms')
            config.natom(nframe) = Data.NumAtoms;        
            natom = Data.NumAtoms;
        end    
        if isfield(Data,'Lattice')            
            config.lattice(:,nframe) = Data.Lattice(:);    
            config.ncl = length(Data.Lattice(:));
        end    
        if isfield(Data,'Energy') 
            config.e(nframe) = Data.Energy;
            config.nce = 1;
        end    
        if isfield(Data,'Stress')
            config.stress(:,nframe) = Data.Stress(:);
            config.ncs = length(Data.Stress(:));
        end        
        if isfield(Data,'Positions')
            config.x(:,(natomall+1):(natomall+natom)) = Data.Positions';
            config.ncx = size(config.x,1);
            config.dim = size(config.x,1);
        end            
        if isfield(Data,'Forces')
            config.f(:,(natomall+1):(natomall+natom)) = Data.Forces';
            config.ncf = size(config.f,1);            
        end            
        if isfield(Data,'Velocities')
            config.v(:,(natomall+1):(natomall+natom)) = Data.Velocities';            
            config.ncv = size(config.v,1);            
        end            
        if isfield(Data,'Charges')
            config.q((natomall+1):(natomall+natom)) = Data.Charges;            
            config.ncq = 1;
        end            
        if isfield(Data,'Masses')
            config.mass((natomall+1):(natomall+natom)) = Data.Masses;            
            config.ncm = 1;
        end                
        if isfield(Data,'Groups')
            config.group((natomall+1):(natomall+natom)) = Data.Groups;            
            config.ncg = 1;
        end                
        if isfield(Data,'Moves')
            config.move((natomall+1):(natomall+natom)) = Data.Moves;            
            config.nco = 1;
        end                
        if isfield(Data,'Tags')
            config.tags((natomall+1):(natomall+natom)) = Data.Tags;            
            config.nci = 1;
        end                        
        if isfield(Data,'Z')
            config.Z((natomall+1):(natomall+natom)) = Data.Z;            
            config.ncz = 1;
        end                        
        if isfield(Data,'AtomTypes')        
            types = cell2mat(Data.AtomTypes);
            b = ones(1,natom);        
            for k=1:length(species)
               ind = string(types) == string(species(k)); 
               b(ind) = k;
            end                                
            config.t((natomall+1):(natomall+natom)) = b;
            config.nct = 1;
        end

        natomall = natomall + natom;
    end
end

config.natomall = natomall; 
config.nconfigs = nframe;

cd(cdir);

