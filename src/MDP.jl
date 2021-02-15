module MDP

using StaticArrays
using SpecialFunctions
using GalacticOptim, Optim
using LinearAlgebra

include("coordinatetransform.jl")
include("forcecalculation.jl")
include("loadinput.jl")


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
    # Input variables ##########################################################
    J, N, Z, NZ, r_N, Ω, Ω′, Ω′′, f_qm, K, L, M, c, w, Δ = load_input()

    # Optimize coeffiecients ###################################################
    c_opt = optimize_coefficients(w, f, f_qm, r_N, NZ, K, L, M, N, J)

    # Print/plot/save results ##################################################
    println("Coefficients:", c_opt)

    # Call LAMMPS ##############################################################
    # TODO
end

end # module
