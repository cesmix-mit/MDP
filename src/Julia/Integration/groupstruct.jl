
mutable struct GroupStruct
    groupnum::Int32        
    iparam::Vector{Int32} 
    fparam::Array{Float64}    
end

function setgroup(groupnum::Int32, iparam::Vector{Int32}, fparam::Array{Float64})    
    group = GroupStruct(groupnum, iparam, fparam)
    return group  
end

function getgroup(group)    
    groupnum = group.groupnum
    iparam = group.iparam
    fparam = group.fparam
    return groupnum, iparam, fparam  
end


