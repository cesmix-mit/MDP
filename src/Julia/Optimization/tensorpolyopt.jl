function gradientdescent(c, yspace, y0=nothing)

    dim = size(yspace,1)
    ymin = yspace[:,1]
    ymax = yspace[:,2]
    N = Int64.(yspace[:,3])
    
    if y0 === nothing
        M = 100;
        y0 = rand(M, dim)
        for i = 1:dim
            y0[:,i] = ymin[i] .+ y0[:,i]*(ymax[i] - ymin[i])
        end
    end
    
    # estimate maximum stepsize
    f, df = mpolyder(y0, c, ymin, ymax, N)
    y = sum(y0.*y0, dims=2) 
    gamma = sum(df.*df, dims=2) 
    alphamax = minimum(sqrt.(y./gamma))
    alphamax = max(alphamax, 10.0)

    beta = 0.8;    
    y = 1.0*y0
    M = size(y,1) 
    y1 = 1.0*y
    r1 = zeros(M)
    maxstep = zeros(M)    
    alpha = alphamax*ones(M)
    a = zeros(dim)
    iter = 0
    while (true)
        iter = iter + 1
        f, df = mpolyder(y, c, ymin, ymax, N)        
        gamma = 0.5*sum(df.*df, dims=2) 
        
        if minimum(gamma[:]) <= 1e-12
            fmin, imin = findmin(f)
            sol = y[imin,:]
            return sol, fmin, iter, y, f, df, y0  
        end

        for m = 1:M
            for i = 1:dim
                a1 = (y[m,i]-ymax[i])/df[m,i]
                a2 = (y[m,i]-ymin[i])/df[m,i]    
                a[i] = (abs(df[m,i])<=1e-12) ? alphamax : max(a1, a2)    
            end
            maxstep[m] = min(minimum(a), alphamax)
            alpha[m] = min(alpha[m], maxstep[m])
            y1[m,:] = y[m,:] - alpha[m]*df[m,:]    
            r1[m] = f[m] - alpha[m]*gamma[m]
        end
                        
        f1 = mpolyeval(y1, c, ymin, ymax, N)
        for m = 1:M
            fm = f1[m]
            rm = r1[m]    
            
            while (fm > rm)
                alpha[m] = beta*alpha[m]
                rm = f[m] - alpha[m]*gamma[m]
                y1[m,:] = y[m,:] - alpha[m]*df[m,:]                    
                tm = mpolyeval(reshape(y1[m,:],(1,dim)), c, ymin, ymax, N)                
                fm = tm[1]
            end            
            y[m,:] = y1[m,:]                         
            alpha[m] = (1.0/beta)*alpha[m]
        end        
    end        
end



    
    