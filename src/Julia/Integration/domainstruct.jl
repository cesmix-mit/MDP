mutable struct Domain    
    box_exist::Int32                          # 0 = not yet created, 1 = exists
    dimension::Int32                          # 2 = 2d, 3 = 3d
    nonperiodic::Int32                        # 0 = periodic in all 3 dims
                                            # 1 = periodic or fixed in all 6
                                            # 2 = shrink-wrap in any of 6
    triclinic::Int32    # 0 = orthog box, 1 = triclinic
    tiltsmall::Int32    # 1 if limit tilt, else 0  
    box_change::Int32           # 1 if any of next 3 flags are set, else 0
    box_change_size::Int32      # 1 if box size changes, 0 if not
    box_change_shape::Int32     # 1 if box shape changes, 0 if not
    box_change_domain::Int32    # 1 if proc sub-domains change, 0 if not
    vdeform::Int32   
    deform_flag::Int32        # 1 if fix deform exist, else 0
    deform_vremap::Int32      # 1 if fix deform remaps v, else 0
    deform_groupbit::Int32    # atom group to perform v remap for    
    boundary::Vector{Int32} #[3*2];    # settings for 6 boundaries
                         # 0 = periodic
                         # 1 = fixed non-periodic
                         # 2 = shrink-wrap non-periodic
                         # 3 = shrink-wrap non-per w/ min
    pbc::Vector{Int32} #[3];          # xyz periodicity as array
    
    face::Array{Float64,2}#[6*3];   # unit normals of 6 box faces
    corners::Array{Float64,2}#[8*3];                     # 8 corner points
    faces::Array{Float64,2};#[6*4*3];   # 4 corner pts of 6 prism faces   
    boxhi::Vector{Float64}#[3];   # orthogonal box global bounds
    boxlo::Vector{Float64}#[3];   # orthogonal box global bounds
    boxtilt::Vector{Float64}#[3]; # 3 tilt factors
    boxhi_lamda::Vector{Float64}#[3]; # lamda box = (0,1)
    boxlo_lamda::Vector{Float64}#[3]; # lamda box = (0,1)
    boxhi_bound::Vector{Float64}#[3];  # bounding box of tilted domain
    boxlo_bound::Vector{Float64}#[3];  # bounding box of tilted domain
    subhi::Vector{Float64}#[3]; # sub-box bounds on this proc
    sublo::Vector{Float64}#[3]; # sub-box bounds on this proc
    ssubhi::Vector{Float64}#[3]; # shifted sub-box on this proc
    ssublo::Vector{Float64}#[3]; # shifted sub-box on this proc
    bsubhi::Vector{Float64}#[3]; # bounding for sub-box on this proc
    bsublo::Vector{Float64}#[3]; # bounding for sub-box on this proc
    subhi_lamda::Vector{Float64}#[3];  # bounds of subbox in lamda
    sublo_lamda::Vector{Float64}#[3];  # bounds of subbox in lamda
    h::Vector{Float64}#[6]               # shape matrix in Voigt ordering
                                      # Voigt = xx,yy,zz,yz,xz,xy
    h_inv::Vector{Float64}#[6];  # inverse of h
    h_rate::Vector{Float64}#[6]; # rate of box size/shape change    

    Domain() = new();
end

function printout(dom::Domain)    
    display("boxlo: "); display(dom.boxlo);
    display("boxhi: "); display(dom.boxhi);
    display("boxtilt: "); display(dom.boxtilt);
    display("sublo: "); display(dom.sublo);
    display("subhi: "); display(dom.subhi);
    display("sublo_lamda: "); display(dom.sublo_lamda);
    display("subhi_lamda: "); display(dom.subhi_lamda);
    display("bsublo: "); display(dom.bsublo);
    display("bsubhi: "); display(dom.bsubhi);
    display("ssublo: "); display(dom.ssublo);
    display("ssubhi: "); display(dom.ssubhi);
    display("h: "); display(dom.h);
    display("h_inv: "); display(dom.h_inv);
end       

function initdomain(boxhi, boxtilt, pbc=[1; 1; 1], bcs=[1; 1; 1; 1; 1; 1])
    
    dom = Domain();
    dim = 3
    dom.boxlo_lamda = zeros(3)
    dom.boxhi_lamda = zeros(3)
    dom.sublo_lamda = zeros(3)
    dom.subhi_lamda = zeros(3)
    dom.h = zeros(6)
    dom.h_inv = zeros(6)
    dom.h_rate = zeros(6)

    dom.box_change = 0
    dom.deform_flag = 0
    dom.vdeform = 0
        
    dom.dimension = dim
    dom.boxlo = [0.0; 0.0; 0.0];
    dom.boxhi = boxhi;
    dom.boxtilt = boxtilt;
    dom.pbc = pbc;
    dom.bcs = bcs;
    
    if (dom.boxtilt[1] == 0 && dom.boxtilt[2] == 0 && dom.boxtilt[3] == 0)
        dom.triclinic = 0;
    else    
        dom.triclinic = 1;
    end
    
    for i=1:3
        dom.boxlo_lamda[i] = 0.0;
        dom.boxhi_lamda[i] = 1.0;                   
        dom.sublo_lamda[i] = 0.0;
        dom.subhi_lamda[i] = 1.0;                            
    end

    return dom    
end
        
function setdomainstruct(pbc, a, b, c)
    
    dom = Domain();
    dim = 3
    dom.boxlo = zeros(dim)
    dom.boxhi = zeros(dim)
    dom.sublo = zeros(dim)
    dom.subhi = zeros(dim)
    dom.bsublo = zeros(dim)
    dom.bsubhi = zeros(dim)
    dom.ssublo = zeros(dim)
    dom.ssubhi = zeros(dim)
    dom.boxtilt = zeros(dim)
    dom.boxlo_bound = zeros(dim)
    dom.boxhi_bound = zeros(dim)
    dom.boxlo_lamda = zeros(dim)
    dom.boxhi_lamda = zeros(dim)
    dom.sublo_lamda = zeros(dim)
    dom.subhi_lamda = zeros(dim)
    dom.h = zeros(6)
    dom.h_inv = zeros(6)
    dom.h_rate = zeros(6)

    dom.box_change = 0
    dom.deform_flag = 0
    dom.vdeform = 0    
    
    dom.dimension = dim
    dom.pbc = pbc    
    dom.boxlo[1] = 0.0;
    dom.boxlo[2] = 0.0;
    dom.boxlo[3] = 0.0;    
    dom.boxhi[1] = a[1];
    dom.boxhi[2] = b[2];
    dom.boxhi[3] = c[3];        
    dom.boxtilt[1] = b[1];
    dom.boxtilt[2] = c[1];
    dom.boxtilt[3] = c[2];    
    for i=1:3
        dom.boxlo_lamda[i] = 0.0;
        dom.boxhi_lamda[i] = 1.0;                   
        dom.sublo_lamda[i] = 0.0;
        dom.subhi_lamda[i] = 1.0;                            
    end
    
    if ((abs(dom.boxtilt[1]) > 1e-10) || (abs(dom.boxtilt[2]) > 1e-10) || (abs(dom.boxtilt[3]) > 1e-10))
        dom.triclinic = 1;
    else
        dom.triclinic = 0;
    end
            
    #dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic)
    dom = SetGlobalBox(dom);  
    dom = SetLocalOrthBox(dom);
    dom = ShiftedSubbox(dom);
    dom = BoundingSubbox(dom);                                                  
    return dom            
end

function setlatticedomain(dom::Domain, lat::Lattice, reg::Region)    

    dim = 3
    dom.dimension = dim
    lat = SetupLattice(lat, unitstyle, dim);

    dom.triclinic = reg.triclinic;
    dom.boxtilt[0] = lat.spacing[0]*reg.boxtilt[0];                     
    dom.boxtilt[1] = lat.spacing[0]*reg.boxtilt[1];                     
    dom.boxtilt[2] = lat.spacing[1]*reg.boxtilt[2];                                     
    for i = 1:3 
        dom.boxhi[i] = lat.spacing[i]*reg.boxhi[i]; 
        dom.boxlo[i] = lat.spacing[i]*reg.boxlo[i];        
        dom.boxlo_lamda[i] = 0.0;
        dom.boxhi_lamda[i] = 1.0;         
        dom.sublo_lamda[i] = 0.0;
        dom.subhi_lamda[i] = 1.0;                            
    end                          
    
    dom = SetGlobalBox(dom);  
    dom = SetLocalOrthBox(dom);
    dom = ShiftedSubbox(dom);
    dom = BoundingSubbox(dom);                  
    lat = LatticeBoundingBox(lat, dom)                                    
    lat = latticebounds(lat)

    # cpuSetGlobalBox(dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, 
    #         dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic);                
    # cpuSetLocalOrthBox(dom.subhi, dom.sublo, dom.boxhi, dom.boxlo, 
    #         dom.subhi_lamda, dom.sublo_lamda, 3);            
    # epsilon = zeros(3)
    # for i=1:3
    #     epsilon[i] = 1e-6*(dom.boxhi[i]-dom.boxlo[i]);
    # end
    # cpuShiftedSubbox(dom.ssublo, dom.ssubhi, dom.boxlo, dom.boxhi, 
    #         dom.boxlo_lamda, dom.boxhi_lamda, dom.sublo, dom.subhi, 
    #         dom.sublo_lamda, dom.subhi_lamda, epsilon, pbc, dom.triclinic);
    # cpuBoundingSubbox(dom.bsublo, dom.bsubhi, dom.sublo, dom.subhi,   
    #         dom.sublo_lamda, dom.subhi_lamda, dom.boxlo, dom.h, dom.triclinic);               
    # cpuLatticeBoundingBox(lat.sublo, lat.subhi, dom.bsublo, dom.bsubhi, 
    #         lat.primitive, lat.rotaterow, lat.primitinv, lat.rotatecol, 
    #         lat.origin, lat.spacing, lat.scale);
    # lat.latticebounds();        

    return dom, lat
end

function CreateLatticeAtom(lat)
    dim = 3    
    natom = lat.natom
    x = zeros(dim*natom)
    t = zeros(natom)
    x, t = AtomLattice(x, t, lat, 0, dim)
    x = reshape(x, (dim, natom))
    return x, t
end

function SetAtomTypeFraction(x, t, atomtypefraction, seed, save)    
    dim, nlocal = size(x)    
    newtype = Int32(atomtypefraction[1]);
    fraction = atomtypefraction[2];
    seed0 = Int32(atomtypefraction[3]);                    
    t, seed, save, count = SetAtomType(x, fraction, t, seed, save, seed0, newtype, dim, nlocal)
    return t, seed, save, count  
end

function CreateInitialVelocity(x, t, createvelocity, atommass, second, seed, save)

    dim, nlocal = size(x)    
    t_desired = createvelocity[1];
    seed0 = Int32(createvelocity[2]);
    dist_flag = Int32(createvelocity[3]);
    sum_flag = Int32(createvelocity[4]);
    loop_flag = Int32(createvelocity[5]);
    #momentum_flag = Int32(createvelocity[6]);         
    #rotation_flag = Int32(createvelocity[7]);     
    mvv2e = createvelocity[8]
    boltz = createvelocity[9]

    mpiRank = 0
    natom=nlocal;
    tdof = (nlocal-1)*dim;    
    tfactor = mvv2e / (tdof * boltz);        

    ilist = zeros(nlocal)
    amass = zeros(nlocal)
    vcm = zeros(3)

    for i=0:nlocal
         ilist[i+1] = i;
    end    
    v = 0.0*x 

    VelocityCreate(x, v, atommass, second, seed, save, ilist, t, seed0, 
            sum_flag, dist_flag, loop_flag, dim, mpiRank, nlocal, natom);
    masstotal = ComputeMass(amass, atommass, t, ilist, nlocal);                                
    ComputeVCM(vcm, v, atommass, masstotal, ilist, t, dim, nlocal);                
    VelocityZeroMomentum(v, vcm, dim, nlocal);
    temp = ComputeTempScalar(v, atommass, tfactor, t, ilist, dim, nlocal);    
    cpuArrayMultiplyScalar(v, sqrt(t_desired/temp), nlocal*dim);          
    
    return v, vcm, masstotal, temp  
end


