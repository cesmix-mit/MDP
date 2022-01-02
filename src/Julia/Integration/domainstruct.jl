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

    function printout()    
        display("boxlo: "); display(boxlo);
        display("boxhi: "); display(boxhi);
        display("boxtilt: "); display(boxtilt);
        display("sublo: "); display(sublo);
        display("subhi: "); display(subhi);
        display("sublo_lamda: "); display(sublo_lamda);
        display("subhi_lamda: "); display(subhi_lamda);
        display("bsublo: "); display(bsublo);
        display("bsubhi: "); display(bsubhi);
        display("ssublo: "); display(ssublo);
        display("ssubhi: "); display(ssubhi);
        display("h: "); display(h);
        display("h_inv: "); display(h_inv);
    end       

    Domain() = new();
end

function initdomain(boxhi, boxtilt, pbc=[1; 1; 1], bcs=[1; 1; 1; 1; 1; 1])

    dom = Domain();
    
    reg.boxlo = [0.0; 0.0; 0.0];
    dom.boxhi = boxhi;
    dom.boxtilt = boxtilt;
    dom.pbc = pbc;
    dom.bcs = bcs;
    
    if (dom.boxtilt[1] == 0 && dom.boxtilt[2] == 0 && dom.boxtilt[3] == 0)
        dom.triclinic = 0;
    else    
        dom.triclinic = 1;
    end
    
    return dom
    
end
        
function setdomain(dom::Domain, config, ci)

    dim = config.dim;    
    dom.boxlo[1] = 0.0;
    dom.boxlo[2] = 0.0;
    dom.boxlo[3] = 0.0;    
    dom.boxhi[1] = config.a[dim*(ci-1)+1];
    dom.boxhi[2] = config.b[dim*(ci-1)+2];
    dom.boxhi[3] = config.c[dim*(ci-1)+3];        
    dom.boxtilt[1] = config.b[dim*(ci-1)+1];
    dom.boxtilt[2] = config.c[dim*(ci-1)+2];
    dom.boxtilt[3] = config.c[dim*(ci-1)+3];    
    
    if ((abs(dom.boxtilt[1]) > 1e-10) || (abs(dom.boxtilt[2]) > 1e-10) || (abs(dom.boxtilt[3]) > 1e-10))
        dom.triclinic = 1;
    else
        dom.triclinic = 0;
    end
            
    cpuSetGlobalBox(dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, 
            dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic);  

    cpuSetLocalOrthBox(dom.subhi, dom.sublo, dom.boxhi, dom.boxlo, 
            dom.subhi_lamda, dom.sublo_lamda, 3);

    epsilon = zeros(3)
    for i = 1:dim
        epsilon[i] = 1e-6*(dom.boxhi[i]-dom.boxlo[i]);
    end

    cpuShiftedSubbox(dom.ssublo, dom.ssubhi, dom.boxlo, dom.boxhi, 
            dom.boxlo_lamda, dom.boxhi_lamda, dom.sublo, dom.subhi, 
            dom.sublo_lamda, dom.subhi_lamda, epsilon, common.pbc, dom.triclinic);

    cpuBoundingSubbox(dom.bsublo, dom.bsubhi, dom.sublo, dom.subhi,   
            dom.sublo_lamda, dom.subhi_lamda, dom.boxlo, dom.h, dom.triclinic);    
            
    return dom            
end

function setdomain(dom::Domain, lat, reg, common)    
    if (common.readlattice==0) 
        dom.triclinic = reg.triclinic;
        for i = 1:3
            dom.boxhi[i] = reg.boxhi[i]; 
            dom.boxlo[i] = reg.boxlo[i]; 
            dom.boxtilt[i] = reg.boxtilt[i]; 
            dom.boxlo_lamda[i] = 0.0;
            dom.boxhi_lamda[i] = 1.0;                    
        end

        cpuSetGlobalBox(dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, 
                dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic);                
        cpuSetLocalOrthBox(dom.subhi, dom.sublo, dom.boxhi, dom.boxlo, 
                dom.subhi_lamda, dom.sublo_lamda, 3);
                
        epsilon = zeros(3)
        for i=1:3
            epsilon[i] = 1e-6*(dom.boxhi[i]-dom.boxlo[i]);
        end

        cpuShiftedSubbox(dom.ssublo, dom.ssubhi, dom.boxlo, dom.boxhi, 
                dom.boxlo_lamda, dom.boxhi_lamda, dom.sublo, dom.subhi, 
                dom.sublo_lamda, dom.subhi_lamda, epsilon, common.pbc, dom.triclinic);
        cpuBoundingSubbox(dom.bsublo, dom.bsubhi, dom.sublo, dom.subhi,   
                dom.sublo_lamda, dom.subhi_lamda, dom.boxlo, dom.h, dom.triclinic);                                                    
    else 
        cpuLattice(lat.basis, lat.primitive, lat.rotaterow, lat.primitinv, 
                lat.rotatecol, lat.origin, lat.spacing, lat.a1, lat.a2, 
                lat.a3, lat.scale, lat.orientx, lat.orienty, lat.orientz, 
                lat.style, common.unitstyle, lat.spaceflag, common.dim);

        dom.triclinic = reg.triclinic;
        dom.boxtilt[0] = lat.spacing[0]*reg.boxtilt[0];                     
        dom.boxtilt[1] = lat.spacing[0]*reg.boxtilt[1];                     
        dom.boxtilt[2] = lat.spacing[1]*reg.boxtilt[2];                                     
        for i = 1:3 
            dom.boxhi[i] = lat.spacing[i]*reg.boxhi[i]; 
            dom.boxlo[i] = lat.spacing[i]*reg.boxlo[i];        
            dom.boxlo_lamda[i] = 0.0;
            dom.boxhi_lamda[i] = 1.0;                    
        end                          
        
        cpuSetGlobalBox(dom.h, dom.h_inv, dom.boxlo_bound, dom.boxhi_bound, 
                dom.boxhi, dom.boxlo, dom.boxtilt, dom.triclinic);                
        cpuSetLocalOrthBox(dom.subhi, dom.sublo, dom.boxhi, dom.boxlo, 
                dom.subhi_lamda, dom.sublo_lamda, 3);
                
        epsilon = zeros(3)
        for i=1:3
            epsilon[i] = 1e-6*(dom.boxhi[i]-dom.boxlo[i]);
        end

        cpuShiftedSubbox(dom.ssublo, dom.ssubhi, dom.boxlo, dom.boxhi, 
                dom.boxlo_lamda, dom.boxhi_lamda, dom.sublo, dom.subhi, 
                dom.sublo_lamda, dom.subhi_lamda, epsilon, common.pbc, dom.triclinic);
        cpuBoundingSubbox(dom.bsublo, dom.bsubhi, dom.sublo, dom.subhi,   
                dom.sublo_lamda, dom.subhi_lamda, dom.boxlo, dom.h, dom.triclinic);                                        
        cpuLatticeBoundingBox(lat.sublo, lat.subhi, dom.bsublo, dom.bsubhi, 
                lat.primitive, lat.rotaterow, lat.primitinv, lat.rotatecol, 
                lat.origin, lat.spacing, lat.scale);
        lat.latticebounds();        
    end

    return dom, lat
end
