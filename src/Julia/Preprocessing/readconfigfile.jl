#***************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
#  
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#***************************************************************************

function readconfigfile(datapath, datafile, dataformat, datamode, atomspecies, transposelattice=nothing)

    #print("Read configuration from " * datapath * "\n");
    filename = datapath * "/" *  datafile;
    style = lowercase(dataformat);
    
    if style == "lammps"     
        config = readLAMMPS(filename,datamode);
    elseif style == "mdp" 
        config = readMDP(filename, datamode);
    elseif style == "extxyz" 
        config = readEXTXYZ(filename, atomspecies);    
    elseif style == "json" 
        config = readJSON(filename, atomspecies);    
    end
    
    config.nconfigs = length(config.natom); # number of configurations
    config.natomall = sum(config.natom);# number of atoms for all configurations    

    dim = config.dim;  
    if transposelattice == 1
        n = size(config.lattice,2)
        if n > 0
            lattice = reshape(config.lattice, (dim, dim, n))
            lattice = permutedims(lattice, (2,1,3))
            config.lattice = reshape(lattice, (dim*dim, n))
        end         
    end
      
    if size(config.lattice,1) == dim*dim    
        config.a = config.lattice[1:dim,:];
        config.b = config.lattice[(dim+1):2*dim,:];
        if dim==3
            config.c = config.lattice[(2*dim+1):3*dim,:];
        end
        n = size(config.lattice,2)
        natom = [0; cumsum(config.natom[:])];

        # rotate lattice and positions 
        for i = 1:n
            A = config.a[:,i]
            B = config.b[:,i]
            C = config.c[:,i]
            a, b, c = rotatelattice(A, B, C)
            config.a[:,i] = a
            config.b[:,i] = b
            config.c[:,i] = c
            X = config.x[:,(natom[i]+1):natom[i+1]]
            x, R = transformcoords(a,b,c,A,B,C,X)
            config.x[:,(natom[i]+1):natom[i+1]] = x
        end
    end

    return config    
end
        
function readconfigfile(app, transposelattice=nothing)

#print("Read configuration from " * app.datapath);
# filename = app.datapath * "/" *  app.datafile;
# style = lowercase(app.dataformat);

# if style == "lammps"     
#     config = readLAMMPS(filename,app.datamode);
# elseif style == "mdp" 
#     config = readMDP(filename, app.datamode);
# elseif style == "extxyz" 
#     config = readEXTXYZ(filename, app.atomspecies);    
# elseif style == "json" 
#     config = readJSON(filename, app.atomspecies);    
# end

# config.nconfigs = length(config.natom); # number of configurations
# config.natomall = sum(config.natom);# number of atoms for all configurations    

# dim = config.dim;    
# if size(config.lattice,1) == dim*dim    
#     config.a = config.lattice[1:dim,:];
#     config.b = config.lattice[(dim+1):2*dim,:];
#     if dim==3
#         config.c = config.lattice[(2*dim+1):3*dim,:];
#     end
# end

config = readconfigfile(app.datapath, app.datafile, app.dataformat, app.datamode, 
            app.atomspecies, transposelattice)

return config

end

function readconfigpath(datapath, dataformat, fileextension, atomspecies, transposelattice=nothing)

    fileextension = lowercase(fileextension)
    datamode = (fileextension == "bin") ? 0 : 1
        
    files = readdir(datapath)

    j = 0;
    for i = 1:length(files)
        datafile = files[i]
        ii = findlast(".", datafile)
        ext = datafile[(ii[end]+1):end]    
        if (fileextension == lowercase(ext)) | (fileextension == "all")
            j = i
            break
        end
    end
    
    datafile = files[j]
    config = readconfigfile(datapath, datafile, dataformat, datamode, atomspecies, transposelattice)

    for i = (j+1):length(files)
        datafile = files[i]
        ii = findlast(".", datafile)
        ext = datafile[(ii[end]+1):end]    
        if (fileextension == lowercase(ext)) | (fileextension == "all")
            config2 = readconfigfile(datapath, datafile, dataformat, datamode, atomspecies, transposelattice)
            config = catconfig(config, config2)
        end
    end
    
    return config
end

function readconfigdata(data, pbc=nothing, a=nothing, b=nothing, c=nothing)

    datapath, dataformat, fileextension, percentage, randomize, atomspecies, 
            weight, translation, rotation, transposelattice  = getdata(data)
    config = readconfigpath(datapath, dataformat, fileextension, atomspecies, transposelattice)
    config.nconfigs = length(config.natom)

    if randomize==1
        ind = randperm(config.nconfigs)
    else
        ind = Array(1:config.nconfigs)
    end
    ntrain = Int64(round(percentage*config.nconfigs/100.0))        
    indices = ind[1:ntrain]    
    
    config = extractconfig(config, indices)
    config.we = reshape([weight[1]],(1,1))
    config.wf = reshape([weight[2]],(1,1))
    config.ws = reshape([weight[3]],(1,1))
    config.nconfigs = length(config.natom)

    if pbc !== nothing
        config.pbc = reshape(pbc, (3,1))
    end
    if a !== nothing        
        config.a = reshape(a, (3,1))
    end
    if b !== nothing
        config.b = reshape(b, (3,1))
    end
    if c !== nothing
        config.c = reshape(c, (3,1))
    end

    # apply rotation and translation to position vectors
    dim = size(config.x,1)
    for i = 1:dim
        config.x[i,:] = config.x[i,:] .+ translation[i]
    end
    config.x = rotation[1:dim,1:dim]*config.x

    return config, indices
end


# function readconfigpaths(paths, dataformat, datamode, atomspecies)

#     for j = 1:length(paths)
#         datapath = paths[j]
#         files = readdir(datapath)
#         for i = 1:length(files)
#             datafile = files[i]
#             config2 = readconfigfile(datapath, datafile, dataformat[j], datamode[j], atomspecies[j])
#             if (i == 1) & (j==1)
#                 config = congig2 
#             else
#                 config = catconfig(config, config2)
#             end
#         end
#     end
    
#     return config
# end

