#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Potential

using Revise, LinearAlgebra

export empiricalpotential, snappotential, lammpspotential, addpotential, getpotential, initsna, PotentialStruct
export empiricaldescriptors, snapdescriptors, lammpssnapdescriptors, emldescriptors
export boundingbox, domainmaps, atomlist, neighborlist, fullneighborlist, neighsingles, accumarray 
export neighpairlist, neighpairs, neightripletlist, neightriplets, neighquadletlist, neighquadlets 
export tallysingle, tallypair, tallytriplet, tallyquadlet, findatomtype, periodicimages

function latticevolume(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64})
    vol = dot(a, cross(b, c))    
    return vol
end

# function getpotentials(potentials)

#     for i = length(potentials):-1:1    
#         if potentials[i] === nothing
#             deleteat!(potentials, i)            
#         end
#     end   

#     n = length(potentials) 
#     empotentials = Array{Any}(nothing, n)
#     mlpotentials = nothing;

#     n = 0
#     for i = 1:length(potentials)    
#         if typeof(potentials[i]) == Potential.PotentialStruct
#             n = n + 1
#             empotentials[n] = potentials[i]                         
#         else
#             mlpotentials = potentials[i]                         
#         end
#     end

#     if n>0
#         if length(potentials) > n
#             for i = length(potentials):-1:(n+1)    
#                 deleteat!(empotentials, i)            
#             end   
#         end     
#     else
#         empotentials = []
#     end

#     return empotentials, mlpotentials
# end

include("accumarray.jl");
include("atomlist.jl");
include("boundingbox.jl");
include("domainmaps.jl");
include("findatomtype.jl");
include("periodicimages.jl");
include("prebox.jl");
include("potentialstruct.jl");
include("latticecoords.jl");
include("neighborlist.jl");
include("fullneighborlist.jl");
include("neighsingles.jl");
include("neighpairlist.jl");
include("neighpairs.jl");
include("neightripletlist.jl");
include("neightriplets.jl");
include("neighquadletlist.jl");
include("neighquadlets.jl");
include("empiricalpotential.jl");
include("empiricaldescriptors.jl");
include("tally.jl");
include("snastruct.jl");
include("snapcpp.jl");
include("lammpspotential.jl");
include("lammpssnapdescriptors.jl");
include("emlpotential.jl");
include("emldescriptors.jl");

end




