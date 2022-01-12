mutable struct Lattice
    basis::Array{Float64,2} #[12*3]; # fractional coords of each basis atom within unit cell (0 <= coord < 1)
    primitive::Vector{Float64} #[9]; # lattice <-> box transform matrices
    primitinv::Vector{Float64} #[9];
    rotaterow::Vector{Float64} #[9];
    rotatecol::Vector{Float64} #[9];
    spacing::Vector{Float64} #[3]; # lattice scale factors in 3 dims
    origin::Vector{Float64} #[3]; # lattice origin
    sublo::Vector{Float64} #[3];  # sub-box bounds in lattice space on this proc
    subhi::Vector{Float64} #[3];  # sub-box bounds in lattice space on this proc
    a1::Vector{Float64} #[3]; # edge vectors of unit cell  
    a2::Vector{Float64} #[3]; # edge vectors of unit cell  
    a3::Vector{Float64} #[3]; # edge vectors of unit cell  
    scale::Float64 

    atomtype::Vector{Int32}   # type of basis atoms
    orientx::Vector{Int32} #[3]; # lattice orientation vecs
    orienty::Vector{Int32} #[3]; # orientx = what lattice dir lies
    orientz::Vector{Int32} #[3]; #           along x dim in box

    style::Int32  # NONE,SC,FCC,etc
    spaceflag::Int32    
    nbasis::Int32   # # of basis atoms in unit cell    
    ilo::Int32
    ihi::Int32
    jlo::Int32 
    jhi::Int32
    klo::Int32
    khi::Int32     # lattice bounds for sub-box on this proc
    natom::Int32                              # # of atoms in the sub-box on this proc

    Lattice() = new();
end

function printout(lat::Lattice)        
    display("style, spaceflag, nbasis, natom, ilo, ihi, jlo, jhi, klo, khi, scale")
    display([lat.style, lat.spaceflag, lat.nbasis, lat.natom, lat.ilo, lat.ihi, lat.jlo, lat.jhi, lat.klo, lat.khi, lat.scale])        
    display("origin: "); display(lat.origin);
    display("spacing: "); display(lat.spacing);
    display("orientx: "); display(lat.orientx);
    display("orienty: "); display(lat.orienty);
    display("orientz: "); display(lat.orientz);
    display("a1: "); display(lat.a1);
    display("a2: "); display(lat.a2);
    display("a3: "); display(lat.a3);
    display("sublo: "); display(lat.sublo);
    display("subhi: "); display(lat.subhi);                
    display("type: "); display(lat.atomtype);
    display("basis: "); display(lat.atombasis);
    display("primitive: "); display(lat.primitive);
    display("primitinv: "); display(lat.primitinv);
    display("rotaterow: "); display(lat.rotaterow);
    display("rotatecol: "); display(lat.rotatecol);
end

function initlattice(stylein, scale, orientx=[1; 0; 0], orienty=[0; 1; 0], orientz=[0; 0; 1], 
        spacing=[1.0; 1.0; 1.0], type=nothing, a1=[1.0; 0.0; 0.0], a2=[0.0; 1.0; 0.0], a3=[0.0; 0.0; 1.0], 
        basis=nothing)

    lat = Lattice();

    lat.scale = scale;
    lat.spaceflag = 0;

    if (lowercase(stylein) == "none") 
        lat.style = 0;
    elseif (lowercase(stylein) == "sc") 
        lat.style = 1;    
    elseif (lowercase(stylein) == "bcc") 
        lat.style = 2;
    elseif (lowercase(stylein) == "fcc") 
        lat.style = 3;
    elseif (lowercase(stylein) == "hcp") 
        lat.style = 4;
    elseif (lowercase(stylein) == "diamond") 
        lat.style = 5;
    elseif (lowercase(stylein) == "sq") 
        lat.style = 6;
    elseif (lowercase(stylein) == "sq2") 
        lat.style = 7;
    elseif (lowercase(stylein) == "hex") 
        lat.style = 8;        
    elseif (lowercase(stylein) == "custom") 
        lat.style = 9;        
    else
        error("Invalid lattice style");
    end

    lat.origin = [0.0; 0.0; 0.0];
    lat.orientx = orientx;
    lat.orienty = orienty;
    lat.orientz = orientz;
    lat.spacing = spacing;

    if (lowercase(stylein) == "custom") 
        lat.a1 = a1;
        lat.a2 = a2;
        lat.a3 = a3;
        lat.basis = basis;    
    else
        lat.a1 = [1.0; 0.0; 0.0];
        lat.a2 = [0.0; 1.0; 0.0];
        lat.a3 = [0.0; 0.0; 1.0];
        
        if (lowercase(stylein) === "hex") 
            lat.a2[2] = sqrt(3.0);
        end
        if (lowercase(stylein) === "hcp") 
            lat.a2[2] = sqrt(3.0);
            lat.a3[3] = sqrt(8.0/3.0);
        end
        
        if (lowercase(stylein) === "sc") 
            lat.basis = [0 0 0];        
        elseif (lowercase(stylein) === "bcc") 
            lat.basis = [0 0 0; 0.5 0.5 0.5];
        elseif (lowercase(stylein) === "fcc") 
            lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5];
        elseif (lowercase(stylein) === "hcp") 
            lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 5.0/6.0 0.5; 0.0 1.0/3.0 0.5];
        elseif (lowercase(stylein) === "diamond") 
            lat.basis = [0 0 0; 0.5 0.5 0.0; 0.5 0.0 0.5; 0.0 0.5 0.5; 0.25 0.25 0.25; 0.25 0.75 0.75; 0.75 0.25 0.75; 0.75 0.75 0.25];        
        elseif (lowercase(stylein) === "sq") 
            lat.basis = [0 0 0];        
        elseif (lowercase(stylein) === "sq2") 
            lat.basis = [0 0 0; 0.5 0.5 0.0];
        elseif (lowercase(stylein) === "hex") 
            lat.basis = [0 0 0; 0.5 0.5 0.0];    
        end
    end

    nbasis = size(lat.basis)[1];
    lat.basis = transpose(lat.basis)

    if type===nothing
        lat.atomtype = Int32.(ones(nbasis));
    else
        lat.atomtype = atomtype;
    end

    return lat
end

function latticebounds(lat::Lattice)    
    lat.ilo = Int32(lat.sublo[1]) - 1;
    lat.jlo = Int32(lat.sublo[2]) - 1;
    lat.klo = Int32(lat.sublo[3]) - 1;
    lat.ihi = Int32(lat.subhi[1]) + 1;
    lat.jhi = Int32(lat.subhi[2]) + 1;
    lat.khi = Int32(lat.subhi[3]) + 1;

    if (lat.sublo[1] < 0.0) 
        lat.ilo = lat.ilo - 1
    end
    if (lat.sublo[2] < 0.0) 
        lat.jlo = lat.jlo - 1
    end        
    if (lat.sublo[3] < 0.0) 
        lat.klo = lat.klo - 1
    end
      
    lat.natom = 0;
    for k = lat.klo:lat.khi
        for j = lat.jlo:lat.jhi
            for i = lat.ilo:lat.ihi
                for m = 0:(lat.nbasis-1) 
                    lat.natom = lat.natom + 1            
                end
            end
        end
    end
end    

#SetupLattice(lat::Lattice, unitstyle::Int32, dim::Int32)




