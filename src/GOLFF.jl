module GOLFF

using StaticArrays
using SpecialFunctions
using GalacticOptim, Optim
using LinearAlgebra

include("coordinatetransform.jl")
include("sphericalharmonics.jl")
include("forcecalculation.jl")


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
    
    # `I[i]`: index of the atomic number
    I = Dict(1 => 1, 8 => 2, 11 => 3, 17 => 3 )

    # `NZ`: number of unique atomic numbers in all configurations
    NZ = length(unique(vcat(unique.(Z)...)))

    # `r_N[j][i]`: position (Cartesian) of atom i in the configuration j
    r_N = Array{Array}(undef, J)
    for j = 1:J
        positions_j = []
        for k = 1:N[j]
            push!(positions_j, Cartesian(rand(), rand(), rand()))
        end
        r_N[j] = positions_j
    end

    # `r_cut`: cut radius needed to calculate the neighbors of each atom 
    r_cut = rand()

    # `Ω[j][i]`: neighbors of the atom i in conf. j
    # `Ω′[j][i][t]`: neighbors of the atom i in conf. j, whose atomic number is t
    # `Ω′′[j][i][t]`: neighbors of the atom i in conf. j, if the atomic number
    #                 of i is t, else it returns empty
    Ω   = [[ [] for i=1:N[j]] for j=1:J]
    Ω′  = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    Ω′′ = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    for j = 1:J
        for i0 = 1:N[j]
            for i1 = 1:N[j]
                if i0 != i1 && norm(r_N[j][i0] - r_N[j][i1]) < r_cut
                    push!(Ω[j][i0], i1)
                    t = I[Z[j][i1]]
                    push!(Ω′[j][i0][t], i1)
                    t = I[Z[j][i1]]
                    push!(Ω′′[j][i0][t], i1)
                end
            end
        end
    end
    
    # `f_qm[j][i]` quantum force associated to the atom i in the configuration j
    f_qm = Array{Array}(undef, J)
    for j = 1:J
        forces_j = []
        for k = 1:N[j]
            push!(forces_j, Cartesian(rand(), rand(), rand()))
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
    h = 0.001
    Δ = Dict()
    Δ['x'] = Cartesian(h, 0.0, 0.0)
    Δ['y'] = Cartesian(0.0, h, 0.0)
    Δ['z'] = Cartesian(0.0, 0.0, h)
    
    return  J, N, Z, NZ, r_N, Ω, Ω′, Ω′′, f_qm, M, c, w, Δ
    
end


"""
    optimize_coefficients(w, f, f_qm, r_N, M, N, J)
    
"""
function optimize_coefficients(w, f, f_qm, r_N, M, N, J)
    # Eq. 4
    cost_function(c, p) = sum([w[j] * 
                               sum([normsq(f(i, j, c, r_N[j]) - f_qm[j][i])
                                    for i=1:N[j]])
                               for j=1:J])

    c0 = zeros(M)
    prob = OptimizationProblem(cost_function, c0)
    sol = solve(prob, NelderMead())
    
    return sol.minimizer
end


"""
    compute(): main function
    
"""
function compute()
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

compute()


end # module
