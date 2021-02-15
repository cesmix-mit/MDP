"""
    P(l, m, x) Associated Legendre Polynomials
    Eq. 8-11 in summary. It is not present in the original manuscript.
    
"""
function P(l, m, x)
    comb(n, p) = factorial(n) / (factorial(p) * factorial(n - p))
    res = 0.0
    if m == 0
        res =   1.0 / 2.0^l * 
                sum([((-1.0)^(l - k) * comb(l, k) * comb(2.0 * k, l)
                      * x^(2.0 * k - l))
                      for k = ceil(l / 2.0):l])
    elseif m > 0
        res =   (-1.0)^m * (1.0 - x^2)^(m / 2.0) * 1.0 / 2.0^l * 
                sum([((-1.0)^(l - k) * comb(l, k) * comb(2.0 * k, l) 
                    * factorial(2.0 * k - l) / factorial(2.0 * k - l - m)
                    * x^(2.0 * k - l - m))
                    for k = ceil((l + m) / 2.0):l])
    else
        m = -m
        res =   (-1.0)^m * (1.0 - x^2)^(m / 2.0) * 1.0 / 2.0^l * 
                sum([((-1.0)^(l - k) * comb(l, k) * comb(2.0 * k, l) 
                    * factorial(2.0 * k - l) / factorial(2.0 * k - l - m)
                    * x^(2.0 * k - l - m))
                    for k = ceil((l + m) / 2.0):l])
        res = -1.0^m * factorial(l - m) / factorial(l + m) * res
    end
    return res
end

"""
    `Y(l, m, θ, ϕ)`: spherical harmonics of degree l and order m.
    Eq. 7 in summary. Eq. 12 in original manuscript.
    
"""
Y(l, m, θ, ϕ) = √((2.0 * l + 1.0) * factorial(l - m) /
                  (4.0 * π * factorial(l + m)) * 
                   P(l, m, cos(θ)) * exp(im * m * ϕ))

"""
    `g(l, k, r)`: radial basis function
    It is not present in original manuscript.
    
"""
g(l, k, r) = sphericalbessely(l, r)

"""
    `u(k, l, m, r))`: three-dimensional basic functions.
    Eq. 6 in summary. Eq. 11 in original manuscript.
    
"""
u(k, l, m, r) = g(l, k, norm(r)) * Y(l, m, convert(Spherical, r).θ, convert(Spherical, r).ϕ)

"""
    `deriv_u(k, l, m, r)`: derivative of the 3D basic functions.
    Eq. 5 in summary. It is not present in the original manuscript.
    
"""
deriv_u(k, l, m, r) = 
     [ u(k, l, m, r + Δ['x']) - u(k, l, m, r - Δ['x']) / (2.0 * norm(Δ['x'])),
       u(k, l, m, r + Δ['y']) - u(k, l, m, r - Δ['y']) / (2.0 * norm(Δ['y'])),
       u(k, l, m, r + Δ['z']) - u(k, l, m, r - Δ['z']) / (2.0 * norm(Δ['z']))]

"""
    `p(i0, i1, k, k′, l, r, j)`: partial derivatives of the power spectrum components.
    Eq. 4 in summary. Eq. 23 and 24 in original manuscript.
    
"""
p(i0, i1, k, k′, l, r, j) = 
     (  sum([( deriv_u(k, l, m, r[i0] - r[i1])
             * sum([u(k′, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
             for m = -l:l])
      + sum([( deriv_u(k′, l, m, r[i0] - r[i1])
             * sum([u(k, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
             for m = -l:l]))

"""
    `deriv_d(t, k, k′, l, r, i, j)`: partial derivatives of the basis function.
    Eq. 3 in summary. Eq. 23 and 24 in original manuscript.
"""
deriv_d(t, k, k′, l, r, i, j) =
    (  sum([ p(i, s, k, k′, l, r, j) for s in Ω′[j][i][t]])
     - sum([ p(s, i, k, k′, l, r, j) for s in Ω′′[j][i][t]]))

"""
    `f(i, j, c, r)`: atomic forces. The `c` vector will be optimized.
    Eq. 2 in summary. Eq. 4 and 5 in original manuscript.
"""
f(i, j, c, r, m = 0) = Cartesian((sum([c[m+=1] * deriv_d(t, k, k′, l, r, i, j)
                         for t = 1:NZ for k = 1:K for k′ = k:K for l = 0:L ]))...)

