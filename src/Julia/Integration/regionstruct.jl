mutable struct Region 
    style::Int32
    triclinic::Int32 
    interior::Int32                     # 1 for interior, 0 for exterior
    scaleflag::Int32                    # 1 for lattice, 0 for box  
    bboxflag::Int32                # 1 if bounding box is computable
    varshape::Int32                # 1 if region shape changes over time
    dynamic::Int32                 # 1 if position/orient changes over time
    moveflag::Int32 
    rotateflag::Int32    # 1 if position/orientation changes
    openflag::Int32                # 1 if any face is open
    
    r::Vector{Float64}                   # distance between particle & surf, r > 0.0
    delx::Vector{Float64}
    dely::Vector{Float64}
    delz::Vector{Float64}    # vector from surface pt to particle
    radius::Vector{Float64}              # curvature of region at contact point
    iwall::Vector{Int32}                # unique id of wall for storing shear history    
    openfaces::Vector{Int32};#[6];           # flags for which faces are open
    tri::Vector{Int32}#[12*3];    # 3 corner pts of 12 triangles (2 per face)
    theta::Float64 # orientation
    rprev::Float64 # speed of time-dependent radius, if applicable
    
    face::Array{Float64,2}#[6*3];   # unit normals of 6 box faces
    corners::Array{Float64,2}#[8*3];                     # 8 corner points
    faces::Array{Float64,2};#[6*4*3];   # 4 corner pts of 6 prism faces   
    boxhi::Vector{Float64}#[3];   # orthogonal box global bounds
    boxlo::Vector{Float64}#[3];   # orthogonal box global bounds
    boxtilt::Vector{Float64}#[3]; # 3 tilt factors
    clo::Vector{Float64}#[3]; # opposite corners of prism
    chi::Vector{Float64}#[3]; # opposite corners of prism   
    scale::Vector{Float64}#[3];
    extent_lo::Vector{Float64}#[3];
    extent_hi::Vector{Float64}#[3];
    a::Vector{Float64}#[3];  # edge vectors of region
    b::Vector{Float64}#[3];  # edge vectors of region
    c::Vector{Float64}#[3];  # edge vectors of region
    h::Vector{Float64}#[6];
    h_inv::Vector{Float64}#[6];    
    dx::Vector{Float64}#[3];    # current displacement 
    v::Vector{Float64}#[3];                 # translational velocity
    rpoint::Vector{Float64}#[3];            # current origin of rotation axis
    omega::Vector{Float64}#[3];             # angular velocity    
    xcenter::Vector{Float64}#[3];    # translated/rotated center of cylinder/sphere (only used if varshape)
    prev::Vector{Float64}#[5];       # stores displacement (X3), angle and if
    
    function printout()    
        display("boxlo: "); display(boxlo);
        display("boxhi: "); display(boxhi);
        display("boxtilt: "); display(boxtilt);
        display("clo: "); display(clo);
        display("chi: "); display(chi);
        display("extent_lo: "); display(extent_lo);
        display("extent_hi: "); display(extent_hi);
        display("a: "); display(a);
        display("b: "); display(b);
        display("c: "); display(c);
        display("scale: "); display(scale);
        display("h: "); display(h);
        display("h_inv: "); display(h_inv);
        display("dx: "); display(dx);
        display("v: "); display(v);
        display("rpoint: "); display(rpoint);
        display("omega: "); display(omega);
        display("xcenter: "); display(xcenter);
        display("prev: "); display(prev);
    end    
   
    Region() = new()
end

function initregion(boxhi, boxtilt=[0.0; 0.0; 0.0])

    reg = Region();
    
    reg.boxlo = [0.0; 0.0; 0.0];
    reg.boxhi = boxhi;
    reg.boxtilt = boxtilt;
    
    if (reg.boxtilt[1] == 0.0 && reg.boxtilt[2] == 0.0 && reg.boxtilt[3] == 0.0)
        reg.triclinic = 0;
    else    
        reg.triclinic = 1;
    end
    
    return reg    
end
    
