module GOLFF

using StaticArrays
using SpecialFunctions
using GalacticOptim, Optim

include("coordinatetransform.jl")
include("sphericalharmonics.jl")

"""
   get_input_parameters()
   
"""
function get_input_parameters()
    #TODO: open input files with actual input information
    
    # Input parameters #

    # Configuration 1: H, H
    # Configuration 2: O, O
    # Configuration 3: H, O, Na, Cl

    # `J`: number of configurations
    J = 3

    # `N[j]`: number of atoms in configuration j
    N = @SArray [2, 2, 4]
    
    # `Z[j][i]`: atomic number of atom i in configuration j
    Z = @SArray [[1, 1], [8, 8], [1, 8, 11, 17]]

    # `NZ`: number of unique atomic numbers in all configurations
    NZ = length(unique(vcat(unique.(Z)...)))

    # `r_N[j][i]`: position (Cartesian) of atom i in the configuration j
    r_N = Array{Array}(undef, J)
    for j = 1:J
        positions_j = []
        for k = 1:N[j]
            push!(positions_j, Cartesian(rand(),rand(),rand()))
        end
        r_N[j] = positions_j
    end

    # `r_cut`: cut radius needed to calculate the neighbors of each atom 
    r_cut = rand()

    # `Ω[j][i]`: neighbors of the atom i in conf. j
    # `Ω′[j][i][t]`: neighbors of the atom i in conf. j, whose atomic number is t
    # `Ω′′[j][i][t]`: neighbors of the atom i in conf. j, if the atomic number of i is t, else it returns empty
    Ω   = [[ [] for i=1:N[j]] for j=1:J]
    Ω′  = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    Ω′′ = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    for j = 1:J
        for i0 = 1:r_N[j]
            for i1 = 1:r_N[j]
                if i0 != i1 && norm(r_N[i0] - r_N[i1]) < r_cut
                    push!(Ω[j][i0], i1)
                    push!(Ω′[j][i0][Z[j][i1]], i1)
                    push!(Ω′′[j][i0][Z[j][i0]], i1)
                end
            end
        end
    end
    
    # `f_qm[j][i]` quantum force associated to the atom i in the configuration j
    f_qm = Array{Array}(undef, J)
    for j = 1:J
        forces_j = []
        for k = 1:N[j]
            push!(forces_j, Cartesian(rand(),rand(),rand()))
        end
        f_qm[j] = forces_j
    end
    
    # `M`: number of basis functions. M must be divisible by NZ.
    M = 12
    

    # `c[m]`: coefficient needed to calculate the potential/force.
    c =  @SArray [10.0, 2.0, 30.0, 20.0, 1.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 12.0]
    
    # `w[j]`: weight associated to the configuration j
    w = @SArray [1.0, 1.0, 1.0]
    
    # `Δ`: finite difference Delta
    Δ = 0.001
    
    return  J, N, Z, NZ, r_N, Ω, Ω′, Ω′′, f_qm, M, c, w, Δ
    
end

"""
    P(l, m, x) Associated Legendre Polynomials
    See https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    
"""
function P(l, m, x)
    comb(n, p) = factorial(n) / (factorial(p) * factorial(n - p))
    res = 0.0
    if m == 0
        res =   1.0 / 2.0^l * 
                sum([(-1.0)^(l - k) * comb(l , k) * comb(2 * k, l)
                      * x^(2 * k - l)
                      for k = ceil(l / 2.0):l])
    else
        res =   (-1.0)^m * (1.0 - x^2)^(m / 2.0) * 1.0 / 2.0^l * 
                sum([(-1.0)^(l - k) * comb(l, k) * comb(2 * k, l) 
                      * factorial(2.0 * k - l) / factorial(2.0 * k - l - m)
                      * x^(2.0 * k - l - m)
                      for k = ceil((l + m) / 2.0):l])
        if m < 0
            res = -1.0^m * factorial(l - m) / factorial(l + m) * res
        end
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
    u(k, l, m, r) = g(l, k, r) * Y(l, m, convert(Spherical, r).θ, convert(Spherical, r).ϕ)

    # `deriv_u(k, l, m, r)`: derivative of the 3D basic functions
    deriv_u(k, l, m, r) = 
         [ u(k, l, m, r + [Δ, 0.0, 0.0]) - u(k, l, m, r - [Δ, 0.0, 0.0]) / (2.0 * Δ),
           u(k, l, m, r + [0.0, Δ, 0.0]) - u(k, l, m, r - [0.0, Δ, 0.0]) / (2.0 * Δ),
           u(k, l, m, r + [0.0, 0.0, Δ]) - u(k, l, m, r - [0.0, 0.0, Δ]) / (2.0 * Δ)]
    
    # `p(n, i, k, k′, l, r, j)`: partial derivatives of the power spectrum components (Eq. 24)
    p(i0, i1, k, k′, l, r, j) = 
            sum([  deriv_u(k, l, m, r[i0] - r[i1])
                 * sum([u(k′, l, m, r[s] - r[i1]) for s = 1:Ω[j][i1]])
                 for m = -l:l])
          + sum([  deriv_u(k′, l, m, r[i0] - r[i1])
                 * sum([u(k, l, m, r[s] - r[i1]) for s = 1:Ω[j][i1]])
                 for m = -l:l])

    # `deriv_d(t, k, k′, l, r, i)`: partial derivatives of the basis function (Eq. 28)
    deriv_d(t, k, k′, l, r, i, j) = 
           sum([ p(i, s, k, k′, l, r, j) for s in Ω′[j][i][t]])
         - sum([ p(s, i, k, k′, l, r, j) for s in Ω′′[j][i][t]])
         
    # pars(m) = [t, k, k′, l]
    pars(m) = [1.0, 1.0, 1.0, 1.0]
    
    # `f(i, j, c, r)`: atomic forces. The `c` vector will be optimized. (Eq. 3).
    f(i, j, c, r) = sum([c[m] * deriv_d(pars(m)...,r, i, j) for m = 1:M])
    
    return f
end


"""
    optimize_coefficients(w, f, f_qm, r_N, M, N, J)
    
"""
function optimize_coefficients(w, f, f_qm, r_N, M, N, J)
    # Eq. 4
    cost_function(c, p) = sum([w[j] * 
                               sum([normsq.(f(i, j, c, r_N[j]) .- f_qm[j]))
                                    for i=1:N[j]])
                               for j=1:J])

    c0 = zeros(M)
    prob = OptimizationProblem(cost_function, c0)
    sol = solve(prob, NelderMead())
    
    return sol.minimizer
end


"""
    GOLFF(): main function
    
"""
function GOLFF()
    # Input variables ##########################################################
    J, N, Z, NZ, r_N, Ω, Ω′, Ω′′, f_qm, M, c, w, Δ = get_input_parameters()

    # Force calculation ########################################################
    f = calculate_forces(M, Ω, Ω′, Ω′′, Δ)

    # Optimize coeffiecients ###################################################
    c_opt = optimize_coefficients(w, f, f_qm, r_N, M, N, J)

    # Print/plot/save results ##################################################
    println("Coefficients:", c_opt)

    # Call LAMMPS ##############################################################
    # TODO
end


end # module
