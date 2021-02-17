include("../src/MDP.jl")
using .MDP
using StaticArrays
using LinearAlgebra
using Test

@testset "1-Atom-Power-Spectrum-Basis-Function" begin

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
    
    I =  @SMatrix [1.0 0.0 0.0
                   0.0 1.0 0.0
                   0.0 0.0 1.0]

    # `r_N[j][i]`: Position (Cartesian) of atom i in the configuration j.
    r_N = Array{Array}(undef, J)
    r_N_rot = Array{Array}(undef, J)
    for j = 1:J
        positions_j = []
        positions_j_rot = []
        c = MDP.Cartesian(rand(), rand(), rand())
        push!(positions_j, c)
        push!(positions_j_rot, MDP.rotc(c, I))
        for k = 2:N[j]
            c = MDP.Cartesian(rand(), rand(), rand())
            push!(positions_j, c)
            push!(positions_j_rot, MDP.rotc(c, I))
        end
        r_N[j] = positions_j
        r_N_rot[j] = positions_j_rot
    end
    
    # `r_cut`: Cut radius needed to calculate the neighbors of each atom. 
    r_cut = rand()

    # `Ω[j][i]`: Neighbors of the atom i in conf. j.
    # `Ω′[j][i][t]`: Neighbors of the atom i in conf. j, whose atomic number type is t.
    # `Ω′′[j][i][t]`: if the atomic number type is t, 
    #                     it returns the neighbors of the atom i in conf. j,
    #                else it returns empty
    Ω   = [[ [] for i=1:N[j]] for j=1:J]
    Ω′  = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    Ω′′ = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    for j = 1:J
        for i0 = 1:N[j]
            for i1 = 1:N[j]
                # TODO: add periodic boundaries
                if i0 != i1 && norm(r_N[j][i0] - r_N[j][i1]) < r_cut
                    push!(Ω[j][i0], i1)
                    t = T[Z[j][i1]]
                    push!(Ω′[j][i0][t], i1)
                    t = T[Z[j][i1]]
                    push!(Ω′′[j][i0][t], i1)
                end
            end
        end
    end
    
    # `Δ`: Finite difference Delta.
    h = 0.001
    Δ = Dict()
    Δ['x'] = MDP.Cartesian(h, 0.0, 0.0)
    Δ['y'] = MDP.Cartesian(0.0, h, 0.0)
    Δ['z'] = MDP.Cartesian(0.0, 0.0, h)
    
    # `K`: ?
    K = 3
    
    # `L`: Degree.
    L = 3

    j = 1
    i = N[j]
    t = 1
    for k = 1:K
        for k′ = k:K
            for l = 0:L
                d1 = MDP.deriv_d(t, k, k′, l, r_N[j], i, j, Ω, Ω′, Ω′′, Δ)
                d2 = MDP.deriv_d(t, k, k′, l, r_N_rot[j], i, j, Ω, Ω′, Ω′′, Δ)
                @test d1 == d2
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
