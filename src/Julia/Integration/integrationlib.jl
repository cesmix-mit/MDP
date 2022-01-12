#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

__precompile__()

using Revise, Libdl

global integrationlib

function loadintegration(mdppath::String)    
    libexist = 0
    if Sys.isapple()
        libpath = mdppath * "src/Julia/Integration/cpuIntegration.dylib";
    elseif Sys.isunix()
        libpath = mdppath * "src/Julia/Integration/cpuIntegration.so";
    elseif Sys.windows()
        libpath = mdppath * "src/Julia/Integration/cpuIntegration.so";
    end
    if isfile(libpath)
        libexist = 1;
    end        
    if libexist==0            
        cdir = pwd();
        cd(mdppath * "src/Julia/Integration")        
        cstr = `g++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuIntegration.cpp -o cpuintegration.o`                
        run(cstr)
        if Sys.isapple()            
            cstr = `g++ -shared cpuintegration.o -o cpuIntegration.dylib`            
        else
            cstr = `g++ -shared cpuintegration.o -o cpuIntegration.so`            
        end 
        run(cstr)
        cd(cdir)                   
    end
    global integrationlib = Libdl.dlopen(libpath)
end

function closeintegration()
    Libdl.dlclose(integrationlib)
end

function SetGlobalBox(dom::Domain)
    ccall(Libdl.dlsym(integrationlib, :cpuSetGlobalBox), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), 
    dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic)
    return dom     
end

function SetLocalOrthBox(dom::Domain)
    ccall(Libdl.dlsym(integrationlib, :cpuSetLocalOrthBox), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), 
    dom.subhi, dom.sublo, dom.boxhi, dom.boxlo, dom.subhi_lamda, dom.sublo_lamda, dom.dimension)
    return dom     
end

function ShiftedSubbox(dom::Domain)
    epsilon = zeros(3)
    for i=1:3
        epsilon[i] = 1e-6*(dom.boxhi[i]-dom.boxlo[i]);
    end
    ccall(Libdl.dlsym(integrationlib, :cpuShiftedSubbox), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint), 
    dom.ssublo, dom.ssubhi, dom.boxlo, dom.boxhi, dom.boxlo_lamda, dom.boxhi_lamda, 
    dom.sublo, dom.subhi, dom.sublo_lamda, dom.subhi_lamda, epsilon, dom.pbc, dom.triclinic)
    return dom     
end

function BoundingSubbox(dom::Domain)
    ccall(Libdl.dlsym(integrationlib, :cpuBoundingSubbox), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint), 
    dom.bsublo, dom.bsubhi, dom.sublo, dom.subhi, dom.sublo_lamda, dom.subhi_lamda, 
    dom.boxlo, dom.h, dom.triclinic)
    return dom     
end

function SetupLattice(lat::Lattice, unitstyle::Int32, dim::Int32)
    basis = zeros(12*3)    

    nbasis = ccall(Libdl.dlsym(integrationlib, :cpuLattice), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    basis, lat.primitive, lat.rotaterow, lat.primitinv, lat.rotatecol, lat.origin, 
    lat.spacing, lat.a1, lat.a2, lat.a3, lat.scale, lat.orientx, lat.orienty, lat.orientz, 
    lat.style, unitstyle, lat.spaceflag, dim)
    
    lat.nbasis = nbasis
    lat.basis = reshape(basis, (3,12))
    return lat     
end

function LatticeBoundingBox(lat::Lattice, dom::Domain)
    ccall(Libdl.dlsym(integrationlib, :cpuLatticeBoundingBox), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble), 
    lat.sublo, lat.subhi, dom.bsublo, dom.bsubhi, lat.primitive, lat.rotaterow, lat.primitinv, 
    lat.rotatecol, lat.origin, lat.spacing, lat.scale)
    return lat     
end

function AtomLattice(x::Vector{Float64}, t::Vector{Int32}, lat::Lattice, nlocal::Int32, dim::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuAtomLattice), Cvoid, (Ptr{Cdouble}, Ptr{Cint},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Cdouble, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), x, t,
    lat.atombasis, lat.primitive, lat.rotaterow, lat.origin, lat.spacing, lat.scale, 
    lat.atomtype, lat.nbasis, nlocal, lat.ilo, lat.ihi, lat.jlo, lat.jhi, lat.klo, lat.khi, dim)
    return x, t     
end

function SetAtomType(x::Vector{Float64}, fraction::Float64, t::Vector{Int32}, seed::Vector{Int32}, 
    save::Vector{Int32}, seed0::Int32, newtype::Int32, dim::Int32, nlocal::Int32)
    count = ccall(Libdl.dlsym(integrationlib, :cpuSetAtomType), Cint, (Ptr{Cdouble}, 
    Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, fraction, t, seed, save, seed0, newtype, dim, nlocal)
    return t, seed, save, count     
end

function VelocityCreate(x::Vector{Float64}, v::Vector{Float64}, mass::Vector{Float64},
    second::Vector{Float64}, seed::Vector{Int32}, save::Vector{Int32}, map::Vector{Int32}, 
    type::Vector{Int32}, seed0::Int32, sum_flag::Int32, dist_flag::Int32, loop_flag::Int32, 
    dim::Int32, mpiRank::Int32,  nlocal::Int32, natoms::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuVelocityCreate), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, 
    Cint, Cint, Cint, Cint, Cint), x, v, mass, second, seed, save, map, type, seed0, sum_flag, 
    dist_flag, loop_flag, dim, mpiRank, nlocal, natoms)
    return x, v, second, seed, save     
end

function ComputeMass(amass::Vector{Float64}, mass::Vector{Float64}, type::Vector{Int32}, 
    ilist::Vector{Int32}, inum::Int32)
    masstotal = ccall(Libdl.dlsym(integrationlib, :cpuComputeMass), Cdouble, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint), amass, mass, type, ilist, inum)
    return masstotal, amass     
end

function ComputeVCM(vcm::Vector{Float64}, v::Vector{Float64}, mass::Vector{Float64}, 
    masstotal::Float64, ilist::Vector{Int32}, type::Vector{Int32}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuComputeVCM), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    vcm, v, mass, masstotal, ilist, type, dim, inum)
    return vcm      
end

function VelocityZeroMomentum(v::Vector{Float64}, vcm::Vector{Float64}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuVelocityZeroMomentum), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Cint, Cint), v, vcm, dim, inum)
    return v, vcm       
end

function ComputeTempScalar(v::Vector{Float64}, mass::Vector{Float64}, tfactor::Float64, 
    type::Vector{Int32}, ilist::Vector{Int32}, dim::Int32, inum::Int32)
    temp = ccall(Libdl.dlsym(integrationlib, :cpuComputeTempScalar), Cdouble, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    v, mass, tfactor, type, ilist, dim, inum)
    return temp  
end

function ComputeKEAtom(ke::Vector{Float64}, mass::Vector{Float64}, v::Vector{Float64}, 
    mvv2e::Float64, type::Vector{Int32}, ilist::Vector{Int32}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuComputeKEAtom), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    ke, mass, v, mvv2e, type, ilist, dim, inum)
end

function InitialIntegrate(x, v, f, mass::Vector{Float64}, dtarray::Vector{Float64}, tarray::Vector{Float64}, eta_mass::Vector{Float64}, 
    eta::Vector{Float64}, eta_dot::Vector{Float64}, eta_dotdot::Vector{Float64}, vlimitsq::Float64, 
    type::Vector{Int32}, ilist::Vector{Int32}, eta_mass_flag::Int32, biasflag::Int32, 
    mtchain::Int32, nc_tchain::Int32, mode::Int32, dim::Int32, inum::Int32)

    ccall(Libdl.dlsym(integrationlib, :cpuInitialIntegrate), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
    x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, vlimitsq,
    type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum)
end

function FinalIntegrate(x, v, f, mass::Vector{Float64}, dtarray::Vector{Float64}, tarray::Vector{Float64}, eta_mass::Vector{Float64}, 
    eta::Vector{Float64}, eta_dot::Vector{Float64}, eta_dotdot::Vector{Float64}, vlimitsq::Float64, 
    type::Vector{Int32}, ilist::Vector{Int32}, eta_mass_flag::Int32, biasflag::Int32, 
    mtchain::Int32, nc_tchain::Int32, mode::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFinalIntegrate), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
    x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, vlimitsq,
    type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, mode, dim, inum)
end

function VelocityRescalingThermostat(v, mass::Vector{Float64}, dtarray::Vector{Float64}, 
    tarray::Vector{Float64}, second::Vector{Float64}, energy::Float64, type::Vector{Int32}, ilist::Vector{Int32}, 
    seed::Vector{Int32}, save::Vector{Int32}, biasflag::Int32, mode::Int32, dim::Int32, inum::Int32)
    energy = ccall(Libdl.dlsym(integrationlib, :cpuVelocityRescalingThermostat), Cdouble, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, 
    Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    v, mass, dtarray, tarray, second, energy, type, ilist, seed, save, biasflag, mode, dim, inum)
    return energy
end

function SetVelocityInitialIntegrate(v, mass::Vector{Float64}, tfactor::Float64, 
    type::Vector{Int32}, ilist::Vector{Int32}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuSetVelocityInitialIntegrate), Cdouble, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    v, mass, tfactor, type, ilist, dim, inum)
end

function SetVelocityInitialIntegrate(x, v, f, mass::Vector{Float64}, fparam::Vector{Float64}, 
    dtf::Float64, dtv::Float64, type::Vector{Int32}, ilist::Vector{Int32}, iparam::Vector{Int32}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuSetVelocityInitialIntegrate), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum)
end

function SetVelocityFinalIntegrate(x, v, f, mass::Vector{Float64}, dtf::Float64, 
    type::Vector{Int32}, ilist::Vector{Int32}, iparam::Vector{Int32}, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuSetVelocityFinalIntegrate), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint), 
    x, v, f, mass, dtf, type, ilist, iparam, dim, inum)
end

function FixWallReflect(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallReflect), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixSetForce(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixSetForce), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixLineForce(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixLineForce), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixPlaneForce(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixPlaneForce), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

# function FixAddForce(x, v, f, 
#     eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
#     ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
#     ccall(Libdl.dlsym(integrationlib, :cpuFixAddForce), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
#     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
#     x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
# end

# function FixDragForce(x, v, f, 
#     eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, box::Vector{Float64},
#     iparam::Vector{Int32}, ilist::Vector{Int32}, pbc::Vector{Int32}, triclinic::Int32, 
#     eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
#     ccall(Libdl.dlsym(integrationlib, :cpuFixDragForce), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
#     Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
#     x, v, f, eatom, vatom, fparam, box, iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum)
# end

function FixWallHarmonic(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallHarmonic), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixWallLJ93(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallLJ93), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixWallLJ126(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallLJ126), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixWallLJ1043(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallLJ1043), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function FixWallMorse(x, v, f, 
    eatom::Vector{Float64}, vatom::Vector{Float64}, fparam::Vector{Float64}, iparam::Vector{Int32}, 
    ilist::Vector{Int32}, eflag_atom::Int32, vflag_atom::Int32, dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuFixWallMorse), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
    x, v, f, eatom, vatom, fparam, iparam, ilist, eflag_atom, vflag_atom, dim, inum)
end

function ApplyPBC(x, v, image::Vector{Int32}, 
    boxhi::Vector{Float64}, boxlo::Vector{Float64}, hi_lambda::Vector{Float64}, 
    lo_lambda::Vector{Float64}, h::Vector{Float64}, h_inv::Vector{Float64},  
    h_rate::Vector{Float64}, pbc::Vector{Int32}, vdeform::Int32, triclinic::Int32, 
    dim::Int32, inum::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuPBC), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cint}, Cint, Cint, Cint, Cint), x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, 
    h_rate, pbc, vdeform, triclinic, dim, inum)
end

function ArrayDistSquareSum(y, x1, x2, m::Int32, n::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuArrayDistSquareSum), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Cint, Cint), y, x1, x2, m, n)
end

function ArrayPlusAtColumnIndex(A, B, colind::Vector{Int32}, 
    m::Int32, n::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuArrayPlusAtColumnIndex), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, 
    Ptr{Cdouble}, Cint, Cint), A, B, colind, m, n)
end

function ArraySumEveryColumn(A, B, m::Int32, n::Int32)
    ccall(Libdl.dlsym(integrationlib, :cpuArraySumEveryColumn), Cvoid, (Ptr{Cdouble}, 
    Ptr{Cdouble}, Cint, Cint), A, B, m, n)
end

