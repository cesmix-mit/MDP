module GOLFF

using StaticArrays
using SpecialFunctions

include("coordinatetransform.jl")
include("sphericalharmonics.jl")

"""
    sphericalbessel(r, x0, L, K)

# Arguments
- `r`: Radius of spherical coordinates
- `x0`: Zeros of the bessely function
- `L`: Order
- `K`: Number of zeros
"""
function sphericalbessel(r,x0,L,K)
    # roots of spherical Bessel functions
    N = L+1;
    n = (1:N)' - 1/2;
    
    nx = length(r)
    
    # Spherical Bessel functions
    g = zeros(nx,K,N)
    for n in 1:N
        for k in 1:K
            x = x0[n,k] * r
            @. g[:,k,n] = -sqrt(pi/(2x)) * bessely(n-1/2,x) # - SpecialFunctions.sphericalbessel(nu-1, x)
        end
    end
end

function random(N)

end

# 1. Calculate potential
# 2. Calculate derivative potential
# 3. Optimization of potential
# 4. Profit! Call lamps


end # module
