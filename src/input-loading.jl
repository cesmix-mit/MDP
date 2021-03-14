"""
   calc_neighbors()
  `Ω[j][i]`: Neighbors of the atom i in conf. j.
  `Ω′[j][i][t]`: Neighbors of the atom i in conf. j, whose atomic number type is t.
  `Ω′′[j][i][t]`: if the atomic number type is t, 
                     it returns the neighbors of the atom i in conf. j,
                 else it returns empty
  `Ω′′′[j][t]`: Atoms of the configuration j whose atomic number types is t.
   
"""
function calc_neighbors(J, N, NZ, Z, T, r_N, r_cut)
    Ω   = [[ [] for i=1:N[j]] for j=1:J]
    Ω′  = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    Ω′′ = [[[ [] for t=1:NZ] for i=1:N[j]] for j=1:J]
    Ω′′′ = [[ [] for t=1:NZ] for j=1:J]
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
    
    for j = 1:J
        for i = 1:N[j]
            t = T[Z[j][i]]
            push!(Ω′′′[j][t], i)
        end
    end
    
    return Ω, Ω′, Ω′′, Ω′′′
end

"""
   load_input()
   
"""
function load_input()
    #TODO: open input files with actual input information

    # Input parameters #

    # Configuration 1: H, H
    # Configuration 2: O, O
    # Configuration 3: H, O, Na, Cl

    # `J`: Number of configurations.
    J = 3

    # `N[j]`: Number of atoms in configuration j.
    N = @SArray [2, 2, 4]
    
    # `Z[j][i]`: Atomic number of atom i in configuration j.
    Z = @SArray [[1, 1], [8, 8], [1, 8, 11, 17]]
    
    # `T[z]`: Type of the atomic number z.
    #         Each atomic number z is indexed in T.
    T = Dict(1 => 1, 8 => 2, 11 => 3, 17 => 3 )

    # `NZ`: Number of different atomic numbers (or atomic types)
    #       present in all configurations.
    NZ = length(unique(vcat(unique.(Z)...)))

    # `r_N[j][i]`: Position (Cartesian) of atom i in the configuration j.
    r_N = Array{Array}(undef, J)
    for j = 1:J
        positions_j = []
        for k = 1:N[j]
            push!(positions_j, Cartesian(rand(), rand(), rand()))
        end
        r_N[j] = positions_j
    end
    
    # `f_qm[j][i]`: Quantum force associated to the atom i in the configuration j.
    f_qm = Array{Array}(undef, J)
    for j = 1:J
        forces_j = []
        for k = 1:N[j]
            push!(forces_j, Cartesian(rand() + im * rand(), rand() + im * rand(),
                                      rand() + im * rand()))
        end
        f_qm[j] = forces_j
    end

    # `r_cut`: Cut radius needed to calculate the neighbors of each atom. 
    r_cut = 10.0

    # Calc. neighbors
    Ω, Ω′, Ω′′, Ω′′′ = calc_neighbors(J, N, NZ, Z, T, r_N, r_cut)
    
    # `L`: Degree.
    L = 1
    
    # `K`: ?
    K = L + 1
    
    # `M`: Number of power spectrum basis functions.
    M_ps = ceil(Int, NZ * K * (K + 1) / 2.0 * (L + 1))
    
    # `M`: Number of power spectrum basis functions.
    M_bs = ceil(Int, NZ * K * (K + 1) / 2.0 * (L + 1)^3)
    
    # `w[j]`: Weight associated to the configuration j.
    w = @SArray [1.0, 1.0, 1.0]
    
    # `Δ`: Finite difference Delta.
    h = 0.001
    Δ = Dict()
    Δ['x'] = Cartesian(h, 0.0, 0.0)
    Δ['y'] = Cartesian(0.0, h, 0.0)
    Δ['z'] = Cartesian(0.0, 0.0, h)
    
    return  J, N, NZ, r_N, Ω, Ω′, Ω′′, Ω′′′, f_qm, K, L, M_ps, M_bs, w, Δ
    
end

