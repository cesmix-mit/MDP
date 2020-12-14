#ifndef __OPUELEMFACENODE
#define __OPUELEMFACENODE

template <typename T> void opuGetElemFaceNodes(T *un, T *u, int *perm, int npe, int npf, int nfe, int nc, int e1, int e2)
{        
    int ndf = npf*nfe;
    int nn = ndf*(e2-e1);    
    int N = nn*nc;
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn; // [0, npf*nfe*(e2-e1)]
        int j = (idx-i)/nn; // [0, nc]
        int k = i%ndf;   //[0, ndf]
        int e = (i-k)/ndf+e1; // [e1, e2]
        int m = k%npf; //[0, npf]
        int n = (k-m)/npf; // [0, nfe]
        un[idx] = u[perm[m+n*npf]+j*npe+e*npe*nc];               
    }        
    // un: npf*nfe*ne*nc
}

template <typename T> void opuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int nn = np*(e2-e1);
    int ncu = nc2-nc1;
    int N = nn*ncu;
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn; // [0, npe*ne]
        int j = (idx-i)/nn; // [0, ncu]
        int k = i%np;   //[0, npe]
        int e = (i-k)/np+e1;        
        un[idx] = u[k+(j+nc1)*np+e*np*nc];        
    }        
}

template <typename T> void opuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2)
{
    int nn = np*(e2-e1);
    int ncu = nc2-nc1;
    int N = nn*ncu;
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;
        int j = (idx-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        u[k+(j+nc1)*np+e*np*nc] = un[idx];        
    }            
}

template <typename T> void opuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int ncu = nc2-nc1;
    int M = np*ncu;
    int N = np*ncu*(e2-e1);
    for (int idx=0; idx<N; idx++)
    {
        int i = idx%M;        
        int k = i%np;
        int j = (i-k)/np+nc1;
        int e = (idx-i)/M+e1;
        un[idx] = u[k+j*np+e*np*nc];        
    }            
}

template <typename T> void opuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int ncu = nc2-nc1;
    int M = np*ncu;
    int N = np*ncu*(e2-e1);
    for (int idx=0; idx<N; idx++)
    {
        int i = idx%M;        
        int k = i%np;
        int j = (i-k)/np+nc1;
        int e = (idx-i)/M+e1;
        u[k+j*np+e*np*nc] = un[idx];        
    }            
}

template <typename T> void opuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        for (int idx = 0; idx<N; idx++)
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int k2 = facecon[2*m+1];
            int m1 = k1%npe;
            int m2 = k2%npe;
            int n1 = (k1-m1)/npe;
            int n2 = (k2-m2)/npe;          
            uh[idx] = 0.5*(udg[m1+j*npe+n1*npe*nc]+udg[m2+j*npe+n2*npe*nc]);
        }                        
    }
    else if (opts==1) {
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;
            uh[idx] = udg[m1+j*npe+n1*npe*nc];
        }                        
    }
    else if (opts==2) {
        for (int idx = 0; idx<N; idx++)       
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k2 = facecon[2*m+1];
            int m2 = k2%npe;
            int n2 = (k2-m2)/npe;            
            uh[idx] = udg[m2+j*npe+n2*npe*nc];
        }                        
    }
}


template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        for (int idx = 0; idx<N; idx++)
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int k2 = facecon[2*m+1];
            int m1 = k1%npe;
            int m2 = k2%npe;
            int n1 = (k1-m1)/npe;
            int n2 = (k2-m2)/npe;                                  
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
            udg[m2+j*npe+n2*npe*nc] = udg[m2+j*npe+n2*npe*nc] + uh[idx];            
        }                        
    }
    else {
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;            
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
        }                        
    }
}

template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts)
{    
    int ne = e2-e1;
    int K = npf*nc;
    int M = npe*nc;
    int N = M*ne;        
    int I = M*e1;
    if (opts==0) {
        for (int idx = 0; idx<N; idx++) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            int q, p, s;
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];
            for (j=0; j<n; j++) {
                q = cole2f1[rowe2f1[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
            }            
            
            i = ent2ind2[l+npe*k];
            e = (i > 0) ? i : 0;
            n = rowe2f2[i+1] - rowe2f2[e];
            for (j=0; j<n; j++) {
                q = cole2f2[rowe2f2[i]+j];
                p = q%npf;              // [1, npf]
                s = (q-p)/npf;          // [1, nf]      
                udg[I+idx] = udg[I+idx] + uh[p+npf*m+K*s]; 
            }                        
        }        
    }
    else {
        for (int idx = 0; idx<N; idx++) {
            int j = idx%M;              //[1, npe*nc]
            int k = (idx-j)/M+e1;       //[1, ne]      
            int l = j%npe;              //[1, npe]
            int m = (j-l)/npe;          //[1, nc] 
            
            int i = ent2ind1[l+npe*k];
            int e = (i > 0) ? i : 0;
            int n = rowe2f1[i+1] - rowe2f1[e];
            for (j=0; j<n; j++) {
                int q = cole2f1[rowe2f1[i]+j];
                int p = q%npf;          // [1, npf]
                int s = (q-p)/npf;          // [1, nf]    
                udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
            }            
        }
    }
}

// % for i=1:npe*nc*ne
// %     j = rem(i-1,npe*nc)+1; % [1,npe*nc]
// %     k = (i-j)/(npe*nc);    % [1,ne]  
// %     l = rem(j-1,npe)+1;    % [1,npe]
// %     p = (j-l)/(npe);       % [1,nc]          
// %     ind = ent2ind1(l+npe*k);
// %     ine = max(ind,1);
// %     n = rowe2f1(ind+1)-rowe2f1(ine);
// %     for m = 1:n
// %         k = cole2f1(rowe2f1(ind)+m);
// %         s = rem(k-1,npf)+1; % [1,npf]
// %         q = (k-s)/npf;      % [1,nf]            
// %         a = s + npf*p + npf*nc*q;
// %         ze(i) = ze(i) - uf(a);            
// %     end    
// % end


template void opuGetElemFaceNodes(double*, double*, int*, int, int, int, int, int, int);
template void opuGetElemNodes(double*, double*, int, int, int, int, int, int);
template void opuPutElemNodes(double*, double*, int, int, int, int, int, int);
template void opuGetElemNodes2(double*, double*, int, int, int, int, int, int);
template void opuPutElemNodes2(double*, double*, int, int, int, int, int, int);
template void opuGetFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void opuPutFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void opuPutFaceNodes(double*, double*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template void opuGetElemFaceNodes(float*, float*, int*, int, int, int, int, int, int);
template void opuGetElemNodes(float*, float*, int, int, int, int, int, int);
template void opuPutElemNodes(float*, float*, int, int, int, int, int, int);
template void opuGetElemNodes2(float*, float*, int, int, int, int, int, int);
template void opuPutElemNodes2(float*, float*, int, int, int, int, int, int);
template void opuGetFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void opuPutFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void opuPutFaceNodes(float*, float*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

#endif



