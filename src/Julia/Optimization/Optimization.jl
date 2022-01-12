#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

module Optimization

using Revise

export setoptim, printerrors, linearfit, linearfit2, optimize, optimize2, validate, validate2
export mpolyeval, mpolyinterp, gradientdescent, num2string, replacestring

function deletenothing(a)
    if a !== nothing
        n = length(a)
        if n > 0
            m = 0
            for i = n:-1:1    
                if a[i] === nothing
                    m = m + 1
                    deleteat!(a, i)            
                end
            end   
            if m==n 
                a = nothing
            end
        end
    end
    return a
end

include("printerrors.jl");
include("OptimStruct.jl");
include("tensorpolyfit.jl");
include("tensorpolyopt.jl");
include("linearfit.jl");
include("maespace.jl");

end





