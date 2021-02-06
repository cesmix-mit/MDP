
"""
    P(l, m, x) Associated Legendre Polynomials
    See https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    
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
    calculate_forces(M, Ω, Ω′, Ω′′, Δ)
    
"""
function calculate_forces(M, Ω, Ω′, Ω′′, Δ)

    # `Y(l, m, theta, phi)`: spherical harmonics of degree l and order m (Eq. 12)
    Y(l, m, θ, ϕ) = √((2 * l + 1) * factorial(l - m) /
                      (4 * π * factorial(l + m))
                       * P(l, m, cos(θ)) * exp(im * m * ϕ))

    # `g(l, k, r)`: radial basis function
    g(l, k, r) = sphericalbessely(l, r)

    # `u(m, l, k, r)`: 3D basic functions (Eq. 11)
    u(k, l, m, r) = g(l, k, norm(r)) * Y(l, m, convert(Spherical, r).θ, convert(Spherical, r).ϕ)

    # `deriv_u(k, l, m, r)`: derivative of the 3D basic functions
    deriv_u(k, l, m, r) = 
         [ u(k, l, m, r + Δ['x']) - u(k, l, m, r - Δ['x']) / (2.0 * norm(Δ['x'])),
           u(k, l, m, r + Δ['y']) - u(k, l, m, r - Δ['y']) / (2.0 * norm(Δ['y'])),
           u(k, l, m, r + Δ['z']) - u(k, l, m, r - Δ['z']) / (2.0 * norm(Δ['z']))]

    # `p(n, i, k, k′, l, r, j)`: partial derivatives of the power spectrum components (Eq. 24)
    p(i0, i1, k, k′, l, r, j) = 
         (  sum([( deriv_u(k, l, m, r[i0] - r[i1])
                 * sum([u(k′, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
                 for m = -l:l])
          + sum([(  deriv_u(k′, l, m, r[i0] - r[i1])
                 * sum([u(k, l, m, r[s] - r[i1]) for s in Ω[j][i1]]))
                 for m = -l:l]))
    
    # `deriv_d(t, k, k′, l, r, i)`: partial derivatives of the basis function (Eq. 28)
    deriv_d(t, k, k′, l, r, i, j) =
        (  sum([ p(i, s, k, k′, l, r, j) for s in Ω′[j][i][t]])
         - sum([ p(s, i, k, k′, l, r, j) for s in Ω′′[j][i][t]]))

    # pars(m) = [t, k, k′, l]
    pars(m) = [1, 1, 1, 1]

    # `f(i, j, c, r)`: atomic forces. The `c` vector will be optimized. (Eq. 3).
    # TODO: verify the lines below: 
    # Original code line: f(i, j, c, r) = sum([c[m] * deriv_d(pars(m)..., r, i, j) for m = 1:M]) 
    # I added Cartesian(real.( ...))
    f(i, j, c, r) = Cartesian(real.(sum([c[m] * deriv_d(pars(m)..., r, i, j) for m = 1:M]))...)

    return f
end

