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

    # `N`: number of atoms
    N = 8
    
    # `J`: number of configurations
    J = 3

    # `Nj`: 1D array of number of atoms for each configuration j
    Nj = @SArray [ 2, 2, 4 ]
    
    # Configuration 1: H, H
    # Configuration 2: O, O
    # Configuration 3: H, O, Na, Cl
    
    # `w`: 1D array of weights for each configuration j
    w = @SArray [ 1.0, 1.0, 1.0 ]
    
    # `T`: is the number of atom types.
    T = 4

    # `M`: number of basis functions. M must be divisible by T.
    M = 12

    # `Z`: 1D array of size T, where Z[t] is the atomic number for atom type t.
    Z =  @SArray [1, 8, 11, 17] # E.g.: H, O, Na, Cl
    
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
    q = @SArray [1, 1, 2, 2, 1, 2, 3, 4]

    # `x`: 1D array indicating the x-component of atom i in configuration j
    # E.g. if we define H is type 1, O is type 2, Na is type 3, and Cl is type 4
    #      x(i=1,j=1)=10.0
    #      x(i=2,j=1)=30.0
    #      x(i=1,j=2)=20.0
    #      x(i=2,j=2)=50.0
    #      x(i=1,j=3)=30.0
    #      x(i=2,j=3)=20.0
    #      x(i=3,j=3)=70.0
    #      x(i=4,j=3)=40.0
    x = @SArray [10.0, 30.0, 20.0, 50.0, 30.0, 20.0, 70.0, 40.0]
    
    # `y`: analogous to `x`, for the y-component
    y = @SArray [0.0, 10.0, 40.0, 10.0, 30.0, 80.0, 100.0, 70.0]
    
    # `x`: analogous to `x`, for the z-component
    z = @SArray [10.0, 30.0, 20.0, 50.0, 30.0, 20.0, 80.0, 100.0]
    
    positions = @SArray Cartesian(x, y, z)
    
    # `f_qm_x`: 1D array indicating the x-component of quantum mechanical force
    #           of atom i in configuration j
    # E.g. if we define H is type 1, O is type 2, Na is type 3, and Cl is type 4
    #      f_qm_x(i=1,j=1)=15.0
    #      f_qm_x(i=2,j=1)=35.0
    #      f_qm_x(i=1,j=2)=25.0
    #      f_qm_x(i=2,j=2)=55.0
    #      f_qm_x(i=1,j=3)=35.0
    #      f_qm_x(i=2,j=3)=25.0
    #      f_qm_x(i=3,j=3)=75.0
    #      f_qm_x(i=4,j=3)=45.0
    f_qm_x = @SArray [15.0, 35.0, 25.0, 55.0, 35.0, 25.0, 75.0, 45.0]
    
    # `f_qm_y`: analogous to `x`, for the y-component
    f_qm_y = @SArray [5.0, 15.0, 45.0, 15.0, 35.0, 85.0, 95.0, 75.0]
    
    # `f_qm_z`: analogous to `x`, for the z-component
    f_qm_z = @SArray [15.0, 35.0, 25.0, 55.0, 35.0, 25.0, 85.0, 90.0]
    
end


"""
    calculate_forces()

"""
function calculate_forces()
    # TODO: complete functions deriv and deriv_d

    # `P(l,m,x)`: Associated Legendre Polynomials
    #             https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    P(l,x) = 1.0 / ( 2^l * factorial(l) ) * deriv(l, (x^2-1)^l )
    deriv_P(l,x) = deriv(m,P(l,x))
    P(l,m,x) = (-1)^m * sin(θ)^m * deriv_P(l,cos(θ))

    # `Y(l,m,theta,phi)`: spherical harmonics of degree l and order m (Eq. 12).
    Y(l,m,θ,ϕ) = sqrt(
                        ( 2 * l + 1 ) * factorial( l - m ) /
                        ( 4 * π * factorial( l + m ) )
                        * P(l,m,cos(θ)) * exp(i * m * ϕ)
                     )

    # `g(l,k,r)`: radial basis function
    g(l,k,r) = sphericalbessely(l,r)

    # `u(m,l,k,r)`: 3D basic functions (Eq. 11).
    u(k,l,m,r) = g(l,k,r) * Y(l,m,convert(Spherical,r).θ,convert(Spherical,r).ϕ) 
    
    # `p(q_p,i_p,k,k_p,l,r_n)`: partial derivatives of the power spectrum components
    p(q_p,i_p,k,k_p,l,r_n) = 

    # `deriv_d(t,k,k_p,l,r_N,r_n)`: partial derivatives of the basis function.
    #                                 (Eq. 28)
    deriv_d(m,r_Nj,r_i) = ?
    deriv_d(t,k,k_p,l,r_N,r_n) = ?
    
    # `f(c,j)`: atomic forces. The `c` vector will be optimized. (Eq. 3).
    f(i,j,c) = sum([ c[m] * deriv_d(m,r_Nj[j],r[i]) for m = 1:M ])
    return f
end


"""
    optimize_coefficients(f,f_qm)

# Arguments
- `w`: 1D array of weights for each configuration j
- `f`: forces (they depend on a coeff. vector c, which will be minimized)
- `f_qm`: quantum mechanical forces
- `M`: number of basis functions. M must be divisible by T.
- `J`: number of configurations 

"""
function optimize_coefficients(w,f,f_qm,M,J)
    # Eq. 4
    cost_function(c,p) = sum([w[j] * sum((abs.(f(c,j) .- f_qm(j)).^2)) for j=1:J])
    
    c0 = zeros(M)
    prob = OptimizationProblem(cost_function,c0)
    sol = solve(prob,NelderMead())
    
    return sol.minimizer
end



"""
    GOLFF(): main function

"""
function GOLFF()
    # Input variables ##########################################################
    N, J, Nj, w, T, M, Z, c, q, x, y, z, f_qm, g, h 
        = get_input_parameters()

    # Force calculation ########################################################
    f = calculate_forces()

    # Optimize coeffiecients ###################################################
    c_opt = optimize_coefficients(w,f,f_qm,M,J)

    # Print/plot/save results ##################################################
    println("Coefficients:",c_opt)

    # Call LAMMPS ##############################################################
    # TODO
end


end # module
