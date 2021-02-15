################################################################################
#
#    Module MDP.jl
#
#    How to use it:
#        julia> include("MDP.jl")
#        julia> using .MDP
#        julia> MDP.compute()
#
################################################################################

module MDP

using StaticArrays
using SpecialFunctions
using GalacticOptim, Optim
using LinearAlgebra

include("coordinate-transform.jl")
include("input-loading.jl")
include("force-calculation.jl")

export compute


"""
    optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J)
    Eq. 1 in summary. Eq. 4 in original manuscript.
    
"""
function optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J)
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
    # Load input ###############################################################
    println("Loading input data..."); flush(stdout)
    J, N, Z, NZ, r_N, Ω, Ω′, Ω′′, f_qm, K, L, M, c, w, Δ = load_input()

    # Force calculation ########################################################
    println("Calculating force..."); flush(stdout)
    f = calc_force(NZ, Ω, Ω′, Ω′′, K, L,  Δ)

    # Optimize coeffiecients ###################################################
    println("Optimizing coefficients..."); flush(stdout)
    c_opt = optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J)

    # Print/plot/save results ##################################################
    println("Finished!"); flush(stdout)
    println("Optimized coefficients:", c_opt); flush(stdout)
    
    # Call LAMMPS ##############################################################
    # TODO
end

end # module
