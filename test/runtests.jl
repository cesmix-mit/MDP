include("../src/MDP.jl")
using .MDP
using StaticArrays
using LinearAlgebra
using Test


@testset "1-Atom-Basis-Function" begin

    # Configuration 1: 101 H atoms

    # `J`: Number of configurations.
    J = 1

    # `N[j]`: Number of atoms in configuration j.
    N = @SArray [101]
    
    # `Z[j][i]`: Atomic number of atom i in configuration j.
    Z = @SArray [ones(101)]
    
    # `T[z]`: Type of the atomic number z.
    #         Each atomic number z is indexed in T.
    T = Dict(1 => 1)

    # `NZ`: Number of different atomic numbers (or atomic types)
    #       present in all configurations.
    NZ = length(unique(vcat(unique.(Z)...)))
    
    # `Δ`: Finite difference Delta.
    h = 0.001
    Δ = Dict()
    Δ['x'] = MDP.Cartesian(h, 0.0, 0.0)
    Δ['y'] = MDP.Cartesian(0.0, h, 0.0)
    Δ['z'] = MDP.Cartesian(0.0, 0.0, h)

    for L = 1:8
        K = L + 1
        
        # Rotation matrix
        euler_angles = [2.0 * π * rand(), π * rand(), 2.0 * π * rand()]
        rot_mat = MDP.euler2rotm(euler_angles...)
        
        # `r_N[j][i]`: Position (Cartesian) of atom i in the configuration j.
        r_N = Array{Array}(undef, J)
        r_N_rot = Array{Array}(undef, J)
        for j = 1:J
            positions_j = []
            positions_j_rot = []
            c = MDP.Cartesian(rand(), rand(), rand())
            push!(positions_j, c)
            push!(positions_j_rot, MDP.rotc(c, rot_mat))
            for k = 2:N[j]
                c = MDP.Cartesian(rand(), rand(), rand())
                push!(positions_j, c)
                push!(positions_j_rot, MDP.rotc(c, rot_mat))
            end
            r_N[j] = positions_j
            r_N_rot[j] = positions_j_rot
        end
        
        # `r_cut`: Cut radius needed to calculate the neighbors of each atom. 
        r_cut = rand()
        
        # Calc. neighbors
        Ω, Ω′, Ω′′, Ω′′′ = MDP.calc_neighbors(J, N, NZ, Z, T, r_N, r_cut)
        
        # Verify power spectrum and bispectrum basis function invariance regarding
        # neighbors position rotations 
        j = 1
        i = N[j]
        t = 1
        for k = 1:K
            for k′ = k:K
                for l = 0:L
                    d1 = MDP.deriv_d_ps(t, k, k′, l, r_N[j], i, j, Ω, Ω′, Ω′′, Δ)
                    d2 = MDP.deriv_d_ps(t, k, k′, l, r_N_rot[j], i, j, Ω, Ω′, Ω′′, Δ)
                    @test d1 == d2
                    
                    for l1 = 0:L
                        for l2 = 0:L
                            d1 = MDP.deriv_d_bs(t, k, k′, l, l1, l2, r_N[j], j, i, Ω, Ω′′′, Δ)
                            d2 = MDP.deriv_d_bs(t, k, k′, l, l1, l2, r_N_rot[j], j, i, Ω, Ω′′′, Δ)
                            @test d1 == d2
                        end
                    end
                end
            end
        end
        
    end

end

@testset "Rotation" begin
    I =  @SMatrix [1.0 0.0 0.0
                   0.0 1.0 0.0
                   0.0 0.0 1.0]

    c = MDP.Cartesian(1.0, 2.0, 3.0)
    c′ = MDP.rotc(c, I)
    @test c == c′

    m = @SMatrix [0.0 0.0 1.0
                  0.0 1.0 0.0
                  1.0 0.0 0.0]

    c′  = MDP.rotc(c, I)
    c′′  = MDP.rotc(c′, I)
    @test c == c′′
end
