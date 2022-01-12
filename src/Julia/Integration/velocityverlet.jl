using Potential, Preprocessing, Optimization 

function VelocityVerlet(x, v, q, t, a, b, c, pbc, app, potentials)

    dim = Int32(size(x,1))
    inum = Int32(length(t))

    ensemblemode = app.ensemblemode;
    time = app.time;             # initial time
    dt = app.dt
    ntimesteps =  app.ntimesteps;       # # time steps
    globalfreq =  app.globalfreq;    
    skin = app.neighskin;
    every = app.neighevery;      # perform neighbor rebuild list for "every" iterations
    delay = app.neighdelay;      # delay neighbor rebuild
    distcheck = app.neighcheck;  # 0 -> neighbor rebuild by delay and every settings,     
    outputflag = app.outputflag  
    eta = app.eta[:]
    kappa = app.kappa[:]
    atommass = [0; app.atommasses[:]]

    if (lowercase(ensemblemode) == "nve") 
        ensemblemodenum = 0;        
    elseif (lowercase(ensemblemode) == "nvelimit") 
        ensemblemodenum = 1;    
    elseif (lowercase(ensemblemode) == "nvt") 
        ensemblemodenum = 2;        
    else
        ensemblemodenum = -1;
    end
    ensemblemodenum = Int32(ensemblemodenum)

    units = Preprocessing.setunits(app.unitstyle) 
    boltz = Int32(units.boltz)
    mvv2e = Int32(units.mvv2e)
    nktv2p = Int32(units.nktv2p)
    
    dtarray, tarray, vlimitsq,  eta_mass, eta_thermostat, eta_dot, eta_dotdot, 
        eta_mass_flag, biasflag, mtchain, nc_tchain = ensemble(app, dim, inum)

    dom = setdomainstruct(pbc, a, b, c)
    boxlo = dom.boxlo
    boxhi = dom.boxhi
    hi_lambda = dom.boxhi_lamda
    lo_lambda = dom.boxlo_lamda
    h = dom.h
    h_inv = dom.h_inv
    h_rate = dom.h_rate 
    vdeform = dom.vdeform  
    triclinic = dom.triclinic
    inv_volume = 1.0 / (dom.h[1] * dom.h[2] * dom.h[3])
    
    rcutmax = skin + Potential.maximumcutoff(potentials)
    x, alist, neighlist, neighnum = Potential.fullneighborlist(x, a, b, c, pbc, rcutmax);    
    etot, eatom, f, vatom, vtot = Potential.emlpotential(x, q, t, alist, neighlist, neighnum, eta, kappa, potentials)
    pe, temp, pres, ke = thermooutput(v, t, alist, atommass, vtot, etot, mvv2e, boltz, nktv2p, inv_volume, outputflag)
    print("Total number of atoms = $inum\n");
    print("Time step     Temperature        Potential energy        Total energy            Pressure \n"); 
    printoutput(0, ntimesteps, temp, pe, ke, pres) 

    xhold = 1.0*x 
    for istep=1:ntimesteps
        time = time + dt        
        gnum = Int32(length(alist))                
        dtarray[6] = istep      
        xstart = 1.0*x

        InitialIntegrate(x, v, f, atommass, dtarray, tarray, eta_mass, eta_thermostat, 
            eta_dot, eta_dotdot, vlimitsq, Int32.(t[:]), Int32.(alist.-1), eta_mass_flag, biasflag, mtchain, 
            nc_tchain, ensemblemodenum, dim, inum);        

        PostInitialIntegration(x, v, f, t, eatom, vatom, atommass, dtarray, app.setvelocity, app.wallreflect)            

        build = checkbuild(x, xhold, skin, delay, every, distcheck, istep, dim, inum)
        if build==1
            print("Rebuild neighbor list at timestep $istep...\n");
            if sum(pbc) > 0
                ApplyPBC(x, v, image, boxhi, boxlo, hi_lambda, lo_lambda, h, h_inv, 
                    h_rate, pbc, vdeform, triclinic, dim, inum)
            end
            x, alist, neighlist, neighnum = Potential.fullneighborlist(x[:,1:inum], a, b, c, pbc, rcutmax);
            xhold = 1.0*x         
        else            
            if sum(pbc) > 0
                # calculate the distance between x and xstart: xstart = x - xstart
                xstart = x - xstart  
                # shift periodic images by the above distance        
                xtm = x[:,(inum+1):end]
                atm = alist[(inum+1):end]            
                ArrayPlusAtColumnIndex(xtm, xstart, Int32.(atm.-1), dim, gnum);                            
                x[:,(inum+1):end] = xtm         
            end
        end

        etot, eatom, f, vatom, vtot = Potential.emlpotential(x, q, t, alist, neighlist, neighnum, eta, kappa, potentials)
        
        PostForceComputation(x, v, f, eatom, vatom, app.setforce, app.lineforce, app.planeforce, 
                    app.wallharmonic, app.walllj126, app.wallmorse)

        FinalIntegrate(x, v, f, atommass, dtarray, tarray, eta_mass, eta_thermostat, eta_dot, 
            eta_dotdot, vlimitsq, Int32.(t[:]), Int32.(alist.-1), eta_mass_flag, biasflag, 
            mtchain, nc_tchain,  ensemblemodenum, dim, inum);        

        PostFinalIntegration(x, v, f, t, atommass, dtarray, tarray, 
                biasflag, app.setvelocity, app.velocityrescaling)

        if (istep % globalfreq == 0)            
            pe, temp, pres, ke = thermooutput(v, t, alist, atommass, vtot, etot, 
                mvv2e, boltz, nktv2p, inv_volume, outputflag)
            printoutput(istep, ntimesteps, temp, pe, ke, pres) 
        end
    end
end

function VelocityVerlet(config, app, potentials, ci=1)

    natom = [0; cumsum(config.natom[:])];
    
    x = config.x[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    v = config.v[:,(natom[ci]+1):natom[ci+1]] # atom positions of configuration ci
    t = config.t[:,(natom[ci]+1):natom[ci+1]] # atom types of configuration ci
    if length(config.q)>0
        q = config.q[:,(natom[ci]+1):natom[ci+1]]
    else
        q = []
    end
    if size(config.a,2) > 1
        a = config.a[:,ci]                        
        b = config.b[:,ci]                        
        c = config.c[:,ci]             
    else
        a = config.a[:,1]                        
        b = config.b[:,1]                        
        c = config.c[:,1]                       
    end
    if size(config.pbc,2) > 1
        pbc = config.pbc[:,ci]           
    else
        pbc = config.pbc[:,1]         
    end

    VelocityVerlet(x, v, q, t, a, b, c, pbc, app, potentials)
end

function printoutput(istep, ntimesteps, temp, pe, ke, pres) 
    n = length(digits(ntimesteps))
    s0 = repeat(' ', n)
    s1 = replacestring(s0, string(istep)) * "        "
    s2 = num2string(temp, 16) * "      "
    s3 = num2string(pe, 16) * "      "
    s4 = num2string(pe+ke, 16) * "      "
    s5 = num2string(pres, 16) * "      "
    print(s1 * s2 * s3 * s4 * s5 * "\n");                    
end

function ensemble(app, dim, inum)

    ensemblemode = app.ensemblemode
    dtarray = zeros(10)
    tarray = zeros(10)
    eta_mass = zeros(10)
    eta_thermostat = zeros(10)
    eta_dot = zeros(10)
    eta_dotdot = zeros(10)

    units = Preprocessing.setunits(app.unitstyle) 
    boltz = units.boltz
    mvv2e = units.mvv2e
    ftm2v = units.ftm2v
    tdof = (inum-1)*dim;    

    dt = app.dt
    dtarray[1] = dt;     # dt
    dtarray[2] = 0.5*dt*ftm2v; # dtf
    dtarray[3] = dt;     # dtv
    dtarray[4] = 0;             # beginstep 
    dtarray[5] = app.ntimesteps; # endstep        
    vlimitsq = 0.0;    
    eta_mass_flag = 1;
    biasflag = 0;
    mtchain = 1;
    nc_tchain = 1;
    if ensemblemode != "nve"        
        if (ensemblemode == "nvelimit")
            nveparam = app.nveparam
            vlimitsq = (nveparam[1]/dt) * (nveparam[1]/dt);
        end
        if (ensemblemode == "nvt")  # NVT ensemble                
            nvtparam = app.nvtparam
            tarray[1] = nvtparam[0]; # temp start
            tarray[2] = nvtparam[1]; # temp stop
            tarray[3] = 1/nvtparam[2]; # temp frequency    
            #delta = (ntimestep - beginstep)/(endstep - beginstep)  
            #tarray[4] = t_start + delta * (t_stop - t_start); # t_target      
            #t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
            #tarray[5] = t_current;       
            tarray[6] = tdof; 
            tarray[7] = boltz; 
            tarray[8] = (length(nvtparam)>3) ? nvtparam[3] : 0.0; # drag factor added to barostat/thermostat         
            tarray[9] = mvv2e; 
            mtchain   = (length(nvtparam)>4) ? Int32(nvtparam[4]) : 3;         
            nc_tchain = (length(nvtparam)>5) ? Int32(nvtparam[5]) : 1;                 
        end
    end

    eta_mass_flag = Int32(eta_mass_flag)
    biasflag = Int32(biasflag)
    mtchain = Int32(mtchain)
    nc_tchain = Int32(nc_tchain)
    return dtarray, tarray, vlimitsq,  eta_mass, eta_thermostat, eta_dot, eta_dotdot, eta_mass_flag, biasflag, mtchain, nc_tchain
end

function thermooutput(v, t, ilist, atommass, vtot, etot, mvv2e, boltz, nktv2p, inv_volume, flag=0)

    dim, inum = size(v)
    tdof = (inum-1)*dim;    
    tfactor = mvv2e/(tdof * boltz)

    # function ComputeKEAtom(ke::Vector{Float64}, mass::Vector{Float64}, v::Vector{Float64}, 
    #     mvv2e::Float64, type::Vector{Int32}, ilist::Vector{Int32}, dim::Int32, inum::Int32)
    
    keatom = zeros(inum)
    ComputeKEAtom(keatom, atommass, v[:], 2.0, Int32.(t[:]), Int32.(ilist .- 1), Int32(dim), Int32(inum));    
    kesum = sum(keatom)

    pe = etot     
    temp = tfactor*kesum
    ke = 0.5 * tdof * boltz * temp
    pres = (tdof * boltz * temp + vtot[1] + vtot[2] + vtot[3]) / 3.0 * inv_volume * nktv2p

    if (flag == 0) 
        pe = pe/inum                
        ke = ke/inum                          
    end

    return pe, temp, pres, ke 
end

function PostInitialIntegration(x, v, f, t, eatom, vatom, mass, 
            dtarray, setvelocity=nothing, wallreflect=nothing)

    dim, anum = size(v)
    alist = Int32.(0:(anum-1))  
    eflag_atom = 0;
    vflag_atom = 0;            
    dtf = dtarray[2]
    dtv = dtarray[3]
    
    x = x[:]
    v = v[:]
    f = f[:]
    if (setvelocity !== nothing) 
        nsetvelocity = length(setvelocity)    
        for  i = 1:nsetvelocity        
            fparam = setvelocity[i].fparam 
            iparam = setvelocity[i].iparam 
            ind =  findall(t == setvelocity[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            SetVelocityInitialIntegrate(x, v, f, mass, fparam, dtf, 
                    dtv, t, ilist, iparam, dim, inum);        
        end
    end
            
    if (wallreflect !== nothing) 
        nwallreflect = length(wallreflect)
        for i = 1:nwallreflect
            fparam = wallreflect[i].fparam 
            iparam = wallreflect[i].iparam 
            ind =  findall(t == wallreflect[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixWallReflect(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end
end

function checkbuild(x, xhold, skin, delay, every, distcheck, istep, dim, inum)
    build = 0;    
    ago = istep + 1;
            
    if ((ago >= delay) && (ago % every == 0)) 
        if (distcheck == 0) 
            build = 1;
        else             
            ds = zeros(inum)
            ArrayDistSquareSum(ds, x[:], xhold[:], dim, inum);
            maxdistsquare = maximum(ds)
            if (maxdistsquare > 0.25*skin*skin)
                build = 1;            
            end
        end  
    end

    return build
end


function PostForceComputation(x, v, f, t, eatom, vatom, setforce=nothing, lineforce=nothing, 
    planeforce=nothing, wallharmonic=nothing, walllj126=nothing, wallmorse=nothing)

    dim, anum = size(v)
    alist = Int32.(0:(anum-1))  

    eflag_atom = 0;
    vflag_atom = 0;            
    
    if (setforce !== nothing) 
        nsetforce = length(setforce)
        for i = 1:nsetforce
            fparam = setforce[i].fparam 
            iparam = setforce[i].iparam 
            ind =  findall(t == setforce[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixSetForce(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end

    if (lineforce !== nothing) 
        nlineforce = length(lineforce)
        for i = 1:nlineforce
            fparam = lineforce[i].fparam 
            iparam = lineforce[i].iparam 
            ind =  findall(t == lineforce[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixLineForce(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end

    if (planeforce !== nothing) 
        nplaneforce = length(planeforce)
        for i = 1:nplaneforce
            fparam = planeforce[i].fparam 
            iparam = planeforce[i].iparam 
            ind =  findall(t == planeforce[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixPlaneForce(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end

    if (wallharmonic !== nothing) 
        nwallharmonic = length(wallharmonic)
        for i = 1:nwallharmonic
            fparam = wallharmonic[i].fparam 
            iparam = wallharmonic[i].iparam 
            ind =  findall(t == wallharmonic[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixWallHarmonic(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end
    
    if (walllj126 !== nothing) 
        nwalllj126 = length(walllj126)
        for i = 1:nwalllj126
            fparam = walllj126[i].fparam 
            iparam = walllj126[i].iparam 
            ind =  findall(t == walllj126[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixWallLJ126(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end

    if (wallmorse !== nothing)
        nwallmorse = length(wallmorse)
        for i = 1:nwallmorse
            fparam = wallmorse[i].fparam 
            iparam = wallmorse[i].iparam 
            ind =  findall(t == wallmorse[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            FixWallMorse(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
        end
    end
end

function PostFinalIntegration(x, v, f, t, mass, dtarray, tarray, 
        biasflag, setvelocity=nothing, velocityrescaling=nothing)

    dim, anum = size(v)
    alist = Int32.(0:(anum-1))  
    energy = 0.0 
    dtf = dtarray[2]
    #dtv = dtarray[3]

    if (setvelocity !== nothing) 
        nsetvelocity = length(setvelocity)    
        for  i = 1:nsetvelocity        
            #fparam = setvelocity[i].fparam 
            iparam = setvelocity[i].iparam 
            ind =  findall(t == setvelocity[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            SetVelocityFinalIntegrate(x, v, f, mass, dtf, t, ilist, iparam, dim, inum);        
        end
    end

    if (velocityrescaling !== nothing)        
        nvelocityrescaling = length(velocityrescaling)      
        second = zeros(1)                                
        seed = Int32.(ones(1))
        save = Int32.(zeros(1))    
        for  i = 1:nvelocityrescaling       
            fparam = velocityrescaling[i].fparam 
            iparam = velocityrescaling[i].iparam 
            energy = fparam[1];
            mode = iparam[1];    
            ind =  findall(t == velocityrescaling[i].groupnum) 
            ilist = alist[ind]
            inum = length(ilist)
            energy = VelocityRescalingThermostat(v, mass, dtarray, tarray, second, energy, 
                    t, ilist, seed, save, biasflag, mode, dim, inum);      
        end      
    end        

    return energy 
end
