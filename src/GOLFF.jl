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
    
    # Degree
    l = ?
    
    # Order
    m = ?
    
    # k and k′
    k = ?
    k′ = ?

    # `N`: number of atoms
    N = 8
    
    # `J`: number of configurations
    J = 3

    # `Nj`: 1D array of number of atoms for each configuration j
    Nj = @SArray [2, 2, 4]
    
    # Configuration 1: H, H
    # Configuration 2: O, O
    # Configuration 3: H, O, Na, Cl
    
    # `w`: 1D array of weights for each configuration j
    w = @SArray [1.0, 1.0, 1.0]
    
    # `T`: is the number of atom types.
    T = 4

    # `M`: number of basis functions. M must be divisible by T.
    M = 12

    # `Z`: 1D array of size T, where Z[t] is the atomic number for atom type t.
    # Z =  @SArray [1, 8, 11, 17] # E.g.: H, O, Na, Cl
    
    # `c`: 1D array of M coefficients needed to calculate the potential/force.
    c =  @SArray [10.0, 2.0, 30.0, 20.0, 1.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 12.0]
    
    # `q`: 1D array indicating the type of atom i in configuration j
    # E.g. if we define H is type 1, O is type 2, Na is type 3, and Cl is type 4
    #      q(i=1,j=1)=1
    #      q(i=2,j=1)=1
    #      q(i=1,j=2)=2
    #      q(i=2,j=2)=2
    #      q(i=1,j=3)=1
    #      q(i=2,j=3)=2
    #      q(i=3,j=3)=3
    #      q(i=4,j=3)=4
    # q = @SArray [1, 1, 2, 2, 1, 2, 3, 4]


    # `Z`: 1D array indicating the type of atom i in configuration j
    # E.g. if we define H is type 1, O is type 2, Na is type 3, and Cl is type 4
    #      q(i=1,j=1)=1
    #      q(i=2,j=1)=1
    #      q(i=1,j=2)=2
    #      q(i=2,j=2)=2
    #      q(i=1,j=3)=1
    #      q(i=2,j=3)=2
    #      q(i=3,j=3)=3
    #      q(i=4,j=3)=4
    Z = @SArray [1, 1, 2, 2, 1, 2, 3, 4]

    NZ = length(unique(Z))

    # `r_N_italic`: array positions (Cartesian) of each atom sorted by configuration
    #        each position is associated to the atom i in configuration j
    #      r_N_italic[1] is the position of atom i=1 in conf. j=1
    #      r_N_italic[2] is the position of atom i=2 in conf. j=1
    #      r_N_italic[3] is the position of atom i=1 in conf. j=2
    #      r_N_italic[4] is the position of atom i=2 in conf. j=2
    #      r_N_italic[5] is the position of atom i=1 in conf. j=3
    #      r_N_italic[6] is the position of atom i=2 in conf. j=3
    #      r_N_italic[7] is the position of atom i=3 in conf. j=3
    #      r_N_italic[8] is the position of atom i=4 in conf. j=3
    #    r_N_italic = Array{Cartesian}(undef, N) 
    #    for i = 1:N
    #        r_N_italic[i] = Cartesian(rand(),rand(),rand())
    #    end
    
    r_Nj = Array{Array}(undef, J)
    for j = 1:J
        positions_j = []
        for k = 1:Nj[j]
            push!(positions_j, Cartesian(rand(),rand(),rand()))
        end
        r_Nj[j] = positions_j
    end

    # `r_cut`: cut radii needed to calculate the neighbors of each atom 
    r_cut = rand()

    # `Ω`   : Ω[i] returns the neighbors of the atom i
    # `Ω′`  : Ω′[i,t] returns the neighbors of the atom i whose type is t
    # `Ω′′` : Ω′′[i,t] returns the neighbors of the atom i if the type of i is t, else it returns empty
    Ω   = Array{Int32}(undef, N)
    Ω′  = Array{Array}(undef, N)
    Ω′′ = [[ [] for j=1:NZ ] for i=1:N]
    for i = 1:N
        for j = 1:N
            if j != i && norm(r_N[i] - r_N[j]) < r_cut
                push!(Ω[i], j)
                push!(Ω′[i,Z[j]], j)
                push!(Ω′′[i,Z[i]], j)
            end
        end
    end
    
    # `f_qm`: array of quantum forces (Cartesian) of each atom
    #         each position is associated to the atom i in configuration j
    # E.g. if we define H is type 1, O is type 2, Na is type 3, and Cl is type 4
    #      f_qm[1] is quantum force associated to atom i=1 in conf. j=1
    #      f_qm[2] is quantum force associated to atom i=2 in conf. j=1
    #      f_qm[3] is quantum force associated to atom i=1 in conf. j=2
    #      f_qm[4] is quantum force associated to atom i=2 in conf. j=2
    #      f_qm[5] is quantum force associated to atom i=1 in conf. j=3
    #      f_qm[6] is quantum force associated to atom i=2 in conf. j=3
    #      f_qm[7] is quantum force associated to atom i=3 in conf. j=3
    #      f_qm[8] is quantum force associated to atom i=4 in conf. j=3
    f_qm = Array{Cartesian}(undef, N) 
    for i = 1:N
        f_qm[i] = Cartesian(rand(),rand(),rand())
    end
    
    return N, J, Nj, w, T, M, Z, c, q, r_N, Ω, f_qm
    
end

"""
    P(l, m, x) Associated Legendre Polynomials
    
    `l`: degree
    `m`: order
    `x`: cos(θ)
    
    see https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
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
    calculate_forces()

"""
function calculate_forces()
    # TODO: complete function deriv_d and p

    # `Y(l, m, theta, phi)`: spherical harmonics of degree l and order m (Eq. 12).
    Y(l, m, θ, ϕ) = sqrt(
                         (2 * l + 1) * factorial(l - m) /
                         (4 * π * factorial(l + m))
                         * P(l, m, cos(θ)) * exp(i * m * ϕ)
                        )

    # `g(l, k, r)`: radial basis function
    g(l, k, r) = sphericalbessely(l, r)

    # `u(m, l, k, r)`: 3D basic functions (Eq. 11).
    u(k, l, m, r) = g(l, k, r) * Y(l, m, convert(Spherical, r).θ, convert(Spherical, r).ϕ) 
    
    # `p(n, i, k, k′, l, r_Nj_j)`: partial derivatives of the power spectrum components. (Eq. 24).
    h = 0.001
    p(n, i, k, k′, l, r) = 
        sum([[ u(k, l, m, r[n] - r[i] + [h, 0.0, 0.0]) - u(k, l, m, r[n] - r[i] - [h, 0.0, 0.0]) / (2.0 * h),
               u(k, l, m, r[n] - r[i] + [0.0, h, 0.0]) - u(k, l, m, r[n] - r[i] - [0.0, h, 0.0]) / (2.0 * h),
               u(k, l, m, r[n] - r[i] + [0.0, 0.0, h]) - u(k, l, m, r[n] - r[i] - [0.0, 0.0, h]) / (2.0 * h)]
             * sum([u(k′, l, m, r[j] - r[i]) for j = 1:Ω[i]])
             for m = -l:l])
      + sum([[ u(k′, l, m, r[n] - r[i] + [h, 0.0, 0.0]) - u(k, l, m, r[n] - r[i] - [h, 0.0, 0.0]) / (2.0 * h),
               u(k′, l, m, r[n] - r[i] + [0.0, h, 0.0]) - u(k, l, m, r[n] - r[i] - [0.0, h, 0.0]) / (2.0 * h),
               u(k′, l, m, r[n] - r[i] + [0.0, 0.0, h]) - u(k, l, m, r[n] - r[i] - [0.0, 0.0, h]) / (2.0 * h)]
             * sum([u(k, l, m, r[j] - r[i]) for j = 1:Ω[i]])
             for m = -l:l])

    # `deriv_d(t, k, k′, l, r_Nj_j, n)`: partial derivatives of the basis function.
    #                                   (Eq. 28)
    deriv_d(t, k, k′, l, r_Nj_j, n) = 
                           sum([ p(n, i, k, k′, l, r_Nj_j) for i in Ω′[n,t]])
                         - sum([ p(j, n, k, k′, l, r_Nj_j) for j in Ω′′[n,t]])
    
    # `f(i, j, c, r_Nj)`: atomic forces. The `c` vector will be optimized. (Eq. 3).
    f(i, j, c, r_Nj) = sum([c[m] * deriv_d(t, k, k′, l, r_Nj[j], i) for m = 1:M])
    
    return f
end


"""
    optimize_coefficients(f,f_qm)

    `w`: 1D array of weights for each configuration j
    `f`: forces (they depend on a coeff. vector c, which will be minimized)
    `f_qm`: quantum mechanical forces
    `M`: number of basis functions. M must be divisible by T.
    `J`: number of configurations 

"""
function optimize_coefficients(w, f, f_qm, M, J)
    # Eq. 4
    cost_function(c, p) = sum([w[j] * sum(normsq.(f(c, j) .- f_qm(j))) for j=1:J])
    
    c0 = zeros(M)
    prob = OptimizationProblem(cost_function,c0)
    sol = solve(prob, NelderMead())
    
    return sol.minimizer
end


"""
    GOLFF(): main function

"""
function GOLFF()
    # Input variables ##########################################################
    N, J, Nj, w, T, M, Z, c, q, r_N, Ω, f_qm = get_input_parameters()

    # Force calculation ########################################################
    f = calculate_forces()

    # Optimize coeffiecients ###################################################
    c_opt = optimize_coefficients(w, f, f_qm, M, J)

    # Print/plot/save results ##################################################
    println("Coefficients:", c_opt)

    # Call LAMMPS ##############################################################
    # TODO
end


end # module
