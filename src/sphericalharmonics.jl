"""
    sphericalharmonics(coords, L)
"""
function sphericalharmonics(coords::Vector{Spherical{T}}, L) where T
    Ylm = Vector{Matrix{Complex{T}}}(undef, L+1)

    N = length(coords)
    P = zeros(T, N, L+1)

    temp = zeros(T, N, L+1)
    temp[:, 1] = 1

    Ylm[1] = sqrt(1/(4*pi)) .* ones(N,1) .+ zeros(N, 1)*im

    the = getfield.(coords, :θ)
    phi = getfield.(coords, :ϕ)
    costhe = cos.(the)
    a = -(sin.(the))

    l = 1
    n = 2*l+1    
    m = 0:l
    C = @. sqrt(n*factorial(l - m + 1)./(4*pi*factorial(l + m + 1)));

    P[:, 1] = costhe
    P[:, 2] = a

    out = zeros(Complex{T}, N, 2)
    @. out[:,1] = C[1] * P[:, 1] * (cos(0)   + sin(0)im);
    @. out[:,2] = C[2] * P[:, 2] * (cos(phi) + sin(phi)im);
    Ylm[2] = out

    for l in 2:L
        out = zeros(Complex{T}, N, l+1)
        n = 2*l+1    
        m = 0:l
        C = @. sqrt(n*factorial(l - m + 1)./(4*pi*factorial(l + m + 1)));

        Pll = P[:, l] 
        tmp[:, l] = Pll 
   
        P[:,l+1] .= (2*(l-1)+1) .* (a     .* Pll)
        P[:,l]   .= (2*(l-1)+1) .* (costhe.* Pll)
        for m in 1:(l-1)
            Pll    = P[:,m]
            tmp[:,m] = Pll;
            P[:,m] = ((2*(l-1)+1).*(costhe.*Pll) - (l+m-2) .* tmp[:,m]) ./ (l-m+1);
        end
       
        for m in 1:(l+1)
            @. out[:, m] = C[m] * P[:, m] * (cos((m-1) * phi) + sin((m-1)*phi)im)
        end
        Ylm[l+1] = out
    end 
    return Ylm
end
