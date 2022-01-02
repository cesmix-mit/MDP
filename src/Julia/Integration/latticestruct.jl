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
    
    function latticebounds()    
        ilo = Int32(sublo[0]) - 1;
        jlo = Int32(sublo[1]) - 1;
        klo = Int32(sublo[2]) - 1;
        ihi = Int32(subhi[0]) + 1;
        jhi = Int32(subhi[1]) + 1;
        khi = Int32(subhi[2]) + 1;
    
        if (sublo[0] < 0.0) 
            ilo = ilo - 1
        end
        if (sublo[1] < 0.0) 
            jlo = jlo - 1
        end        
        if (sublo[2] < 0.0) 
            klo = klo - 1
        end
          
        natom = 0;
        for k = klo:khi
            for j = jlo:jhi
                for i = ilo:ihi
                    for m = 0:(nbasis-1) 
                        natom = natom + 1            
                    end
                end
            end
        end
    end    

    function printout()        
        display("style, spaceflag, nbasis, natom, ilo, ihi, jlo, jhi, klo, khi, scale")
        display([style, spaceflag, nbasis, natom, ilo, ihi, jlo, jhi, klo, khi, scale])        
        display("origin: "); display(origin);
        display("spacing: "); display(spacing);
        display("orientx: "); display(orientx);
        display("orienty: "); display(orienty);
        display("orientz: "); display(orientz);
        display("a1: "); display(a1);
        display("a2: "); display(a2);
        display("a3: "); display(a3);
        display("sublo: "); display(sublo);
        display("subhi: "); display(subhi);                
        display("type: "); display(atomtype);
        display("basis: "); display(atombasis);
        display("primitive: "); display(primitive);
        display("primitinv: "); display(primitinv);
        display("rotaterow: "); display(rotaterow);
        display("rotatecol: "); display(rotatecol);
    end

    Lattice() = new();
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

    if type===nothing
        lat.type = Int64.(ones(1,nbasis));
    else
        lat.type = type;
    end

    return lat

end



