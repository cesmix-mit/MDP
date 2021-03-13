
################################################################################
# Force calculated based on the power spectrum functions 
################################################################################

"""
    P(l, m, x) Associated Legendre Polynomials
    Eq. 9-12 in summary. It is not present in the original manuscript.
    
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
    Eq. 8 in summary. Eq. 12 in original manuscript.
    
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
    Eq. 7 in summary. Eq. 11 in original manuscript.
    
"""
u(k, l, m, r) = g(l, k, norm(r)) * Y(l, m, convert(Spherical, r).θ, convert(Spherical, r).ϕ)

"""
    `deriv_u(k, l, m, r)`: derivative of the 3D basic functions.
    Eq. 6 in summary. It is not present in the original manuscript.
    
"""
deriv_u(k, l, m, r, Δ) = 
     [ u(k, l, m, r + Δ['x']) - u(k, l, m, r - Δ['x']) / (2.0 * norm(Δ['x'])),
       u(k, l, m, r + Δ['y']) - u(k, l, m, r - Δ['y']) / (2.0 * norm(Δ['y'])),
       u(k, l, m, r + Δ['z']) - u(k, l, m, r - Δ['z']) / (2.0 * norm(Δ['z']))]

"""
    `p(i0, i1, k, k′, l, r, j)`: partial derivatives of the power spectrum components.
    Eq. 5 in summary. Eq. 23 and 24 in original manuscript.
    
"""
p(i0, i1, k, k′, l, r, j, Ω, Δ) = 
     (  sum([( deriv_u(k, l, m, r[i0] - r[i1], Δ)
             * sum([u(k′, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
             for m = -l:l])
      + sum([( deriv_u(k′, l, m, r[i0] - r[i1], Δ)
             * sum([u(k, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
             for m = -l:l]))

"""
    `deriv_d_ps(t, k, k′, l, r, i, j)`: partial derivatives of the basis function.
    Eq. 4 in summary. Eq. 23 and 24 in original manuscript.
"""
deriv_d_ps(t, k, k′, l, r, i, j, Ω, Ω′, Ω′′, Δ) =
    (  sum([ p(i, s, k, k′, l, r, j, Ω, Δ) for s in Ω′[j][i][t]])
     - sum([ p(s, i, k, k′, l, r, j, Ω, Δ) for s in Ω′′[j][i][t]]))

"""
    `f_ps(i, j, c, r, NZ, K, L, Ω, Ω′, Ω′′, Ω′′′, Δ, m = 0)`:
     Atomic force based on power spectrum basis functions.
     The `c` vector will be optimized.
     Eq. 2 in summary. Eq. 4 and 5 in original manuscript.
"""
f_ps(i, j, c, r, NZ, K, L, Ω, Ω′, Ω′′, Ω′′′, Δ, m = 0) =
     Cartesian((sum(c[m+=1] * deriv_d_ps(t, k, k′, l, r, i, j, Ω, Ω′, Ω′′, Δ)
                     for t = 1:NZ for k = 1:K for k′ = k:K for l = 0:L))...)

################################################################################
# Force calculated based on the bispectrum spectrum functions 
################################################################################

"""
    `a(s, k', l, m, r, j)`
    Eq. 17 in summary. Eq. 10 in original manuscript.
"""
a_bs(i, k, l, m, r, j, Ω) =
    sum([u(k, l, m, r[s] - r[i]) for s in Ω[j][i]])
    


"""
    `b(t, k, k′, l, r, i, j)`
    Eq. 16 in summary. Eq. 29 in original manuscript.
"""
b_bs(s, k, k′, l, l1, l2, r, j, i, Ω) =
    sum([conj(a_bs(s, k, l, m, r, j, Ω)) * CG(l1,l2,m1,m2,l,m) *
         a_bs(s, k', l1, m1, r, j, Ω) * a_bs(s, k', l2, m2, r, j, Ω)
         for m = -l:l for m1 = -l1:l1 for m2 = -l2:l2])

"""
    `d_bs(k, k′, l, l1, l2, r, j, i, Ω′′′, r)`: bispectrum basis functions
    Eq. 15 in summary. Eq. 30 in original manuscript.
"""
d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′) =
    sum([b_bs(s, k, k′, l, l1, l2, r, j, i, Ω) for s in Ω′′′[j][t]])

"""
    `deriv_d_bs(t, k, k′, l, r, j, i)`: partial derivatives of the basis function.
    Eq. 14 in summary. Eq. 23 and 24 in original manuscript.
"""
function deriv_d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′, Δ)
    res = [0.0+0.0im, 0.0+0.0im, 0.0+0.0im]
    # d(d_bs)/dx
    h = norm(Δ['x'])
    r[i].x = r[i].x + h
    a = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    r[i].x = r[i].x - 2.0 * h
    b = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    res[1] = (a - b) / (2.0 * h)
    # d(d_bs)/dy
    h = norm(Δ['y'])
    r[i].y = r[i].y + h
    a = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    r[i].y = r[i].y - 2.0 * h
    b = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    res[2] = (a - b) / (2.0 * h)
    # d(d_bs)/dz
    h = norm(Δ['z'])
    r[i].z = r[i].z + h
    a = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    r[i].z = r[i].z - 2.0 * h
    b = d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′)
    res[3] = (a - b) / (2.0 * h)
    # res = (d(d_bs)/dx, d(d_bs)/dy, d(d_bs)/dz)
    return res
end

"""
    `f_bs(i, j, c, r, K, L, Ω, Ω′, Ω′′, m = 0)`:
     Atomic force based on bispectrum basis functions.
     The `c` vector will be optimized.
     Eq. 3 in summary. Eq. 4 and 5 in original manuscript.
"""
f_bs(i, j, c, r, NZ, K, L, Ω, Ω′, Ω′′, Ω′′′, Δ, m = 0) =
     Cartesian((sum([c[m+=1] * deriv_d_bs(t, k, k′, l, l1, l2, r, j, i, Ω, Ω′′′, Δ)
                     for t = 1:NZ for k = 1:K for k′ = k:K for l = 0:L for l1 = 0:L for l2 = 0:L]))...)


