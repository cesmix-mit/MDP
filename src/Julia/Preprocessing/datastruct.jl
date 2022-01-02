mutable struct DataStruct

    datapath::String 
    dataformat::String 
    fileextension::String     
    percentage::Float64
    randomize::Bool 
    atomspecies::Vector{String}
    weight::Vector{Float64}    
    translation::Vector{Float64}    
    rotation::Array{Float64,2}
    transposelattice::Bool 

end


function adddata(datapath::String, dataformat::String, fileextension::String, percentage::Float64, 
    randomize::Bool, atomspecies::Vector{String}, weight=nothing, translation=nothing, 
    rotation=nothing, transposelattice=nothing)

    if weight === nothing
        weight = [1.0; 1.0; 1.0]
    end
    if translation === nothing
        translation = [0.0; 0.0; 0.0]
    end
    if rotation === nothing
        rotation = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    end
    if transposelattice === nothing
        transposelattice = false
    end

    data = DataStruct(datapath, dataformat, fileextension, percentage, 
                randomize, atomspecies, weight, translation, rotation, transposelattice)

    return data
end

function getdata(data::DataStruct)

    datapath = data.datapath
    dataformat = data.dataformat
    fileextension = data.fileextension
    atomspecies = data.atomspecies 
    percentage = data.percentage
    randomize = data.randomize
    weight = data.weight
    translation = data.translation
    rotation = data.rotation
    transposelattice = data.transposelattice

    return datapath, dataformat, fileextension, percentage, randomize, atomspecies, weight, translation, rotation, transposelattice
end

