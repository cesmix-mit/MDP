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
using PartialWaveFunctions

include("coordinate-transform.jl")
include("input-loading.jl")
include("force-calculation.jl")

"""
    optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J, Ω, Ω′, Ω′′, Δ)
    Eq. 1 in summary. Eq. 4 in original manuscript.
    
"""
function optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J, Ω, Ω′, Ω′′, Ω′′′, Δ)
    cost_function(c, p) =
        sum([w[j] * 
             sum([normsq(f(i, j, c, r_N[j], NZ, K, L, Ω, Ω′, Ω′′, Ω′′′, Δ) - f_qm[j][i])
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

    println("Calculation of the force field:"); flush(stdout)

    # Load input ###############################################################
    println("Loading input data..."); flush(stdout)
    J, N, NZ, r_N, Ω, Ω′, Ω′′, Ω′′′, f_qm, K, L, M_ps, M_bs, w, Δ = load_input()

    # Optimize coeffiecients ###################################################
    println("Power spectrum functions case. Optimizing coefficients..."); flush(stdout)
    c_opt_ps = optimize_coefficients(w, f_ps, f_qm, r_N, NZ, K, L, M_ps, N, J, Ω, Ω′, Ω′′, Ω′′′, Δ)
    println("Finished!"); flush(stdout)
    println("Optimized coefficients:", c_opt_ps); flush(stdout)

#    println("Bispectrum functions case. Optimizing coefficients..."); flush(stdout)
#    c_opt_bs = optimize_coefficients(w, f_bs, f_qm, r_N, NZ, K, L, M_bs, N, J, Ω, Ω′, Ω′′, Ω′′′, Δ)
#    println("Finished!"); flush(stdout)
#    println("Optimized coefficients:", c_opt_bs); flush(stdout)

    # Call LAMMPS ##############################################################
    # TODO
end

end # module
