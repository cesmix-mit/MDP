/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDPSNAP
#define MDPSNAP

void SetSnap(snastruct &sna, int twojmax, int ntypes, int backend)
{
    sna.twojmax = twojmax;
    sna.ntypes = ntypes;
    
    int jdim = twojmax + 1;    
    int jdimpq = twojmax + 2;  
    
    TemplateMalloc(&sna.map, ntypes+1, backend);
    TemplateMalloc(&sna.idxcg_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxz_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxb_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxu_block, jdim, backend);   
    TemplateMalloc(&sna.idx_max, 5, backend);     
    
    int idxb_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= min(twojmax, j1 + j2); j += 2)
          if (j >= j1) idxb_count++;
    
    int idxz_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= min(twojmax, j1 + j2); j += 2)
          for (int mb = 0; 2*mb <= j; mb++)
            for (int ma = 0; ma <= j; ma++)
              idxz_count++;
    
    int idxcg_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= min(twojmax, j1 + j2); j += 2) {
                for (int m1 = 0; m1 <= j1; m1++)
                  for (int m2 = 0; m2 <= j2; m2++)
                    idxcg_count++;
            }
   
    TemplateMalloc(&sna.idxz, idxz_count*10, backend);        
    TemplateMalloc(&sna.idxb, idxb_count*3, backend);        
    
    TemplateMalloc(&sna.rcutsq, (ntypes+1)*(ntypes+1), backend);
    TemplateMalloc(&sna.radelem, ntypes+1, backend);
    TemplateMalloc(&sna.wjelem, ntypes+1, backend);
        
    TemplateMalloc(&sna.rootpqarray, jdimpq*jdimpq, backend);      
    TemplateMalloc(&sna.cglist, idxcg_count, backend);        
    TemplateMalloc(&sna.bzero, jdim, backend); 
    TemplateMalloc(&sna.fac, 168, backend); 
    
    if (backend <= 1) {
        for (int i=0; i<168; i++)
            sna.fac[i] = (dstype) factable[i];
        cpuInitSna(sna.rootpqarray, sna.cglist, sna.fac, sna.idx_max, sna.idxz, 
            sna.idxz_block, sna.idxb, sna.idxb_block, sna.idxu_block, sna.idxcg_block, sna.twojmax);    
    }
    else {
        //printf("184 \n");
        dstype *rootpqarray, *cglist, *fac;
        int *idxcg_block, *idxz_block, *idxb_block, *idxu_block, *idx_max, *idxz, *idxb;
                
        TemplateMalloc(&rootpqarray, jdimpq*jdimpq, 0);      
        TemplateMalloc(&cglist, idxcg_count, 0);        
        TemplateMalloc(&idxcg_block, jdim*jdim*jdim, 0);
        TemplateMalloc(&idxz_block, jdim*jdim*jdim, 0);
        TemplateMalloc(&idxb_block, jdim*jdim*jdim, 0);
        TemplateMalloc(&idxu_block, jdim, 0);   
        TemplateMalloc(&idx_max, 5, 0);     
        TemplateMalloc(&idxz, idxz_count*10, 0);        
        TemplateMalloc(&idxb, idxb_count*3, 0);        
        
        TemplateMalloc(&fac, 168, 0); 
        for (int i=0; i<168; i++)
            fac[i] = (dstype) factable[i];
        cpuInitSna(rootpqarray, cglist, fac, idx_max, idxz, 
          idxz_block, idxb, idxb_block, idxu_block, idxcg_block, twojmax);    
        
        TemplateCopytoDevice(sna.rootpqarray, rootpqarray, jdimpq*jdimpq, backend);        
        TemplateCopytoDevice(sna.cglist, cglist, idxcg_count, backend);              
        TemplateCopytoDevice(sna.idxcg_block, idxcg_block, jdim*jdim*jdim, backend);
        TemplateCopytoDevice(sna.idxz_block, idxz_block, jdim*jdim*jdim, backend);
        TemplateCopytoDevice(sna.idxb_block, idxb_block, jdim*jdim*jdim, backend);
        TemplateCopytoDevice(sna.idxu_block, idxu_block, jdim, backend);   
        TemplateCopytoDevice(sna.idx_max, idx_max, 5, backend);         
        TemplateCopytoDevice(sna.idxz, idxz, idxz_count*10, backend);     
        TemplateCopytoDevice(sna.idxb, idxb, idxb_count*3, backend);              
    }
}

void ComputePairSnap(dstype *eatom, dstype *fatom, dstype *vatom, snastruct &sna, commonstruct &common, 
        sysstruct &sys, neighborstruct &nb, tempstruct &tmp)
{    
    int inum = common.inum;
    //int anum = common.anum;
    int neighmax = common.neighmax;    
    int dim = common.dim;
    int backend = common.backend;
    
    int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int ncoeffall = sna.ncoeffall;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int quadraticflag = sna.quadraticflag;
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    dstype wself = sna.wself;
    dstype rmin0 = sna.rmin0;
    dstype rfac0 = sna.rfac0;
    dstype rcutfac = sna.rcutfac;
    dstype rcutmax = sna.rcutmax;        
    dstype *bzero = sna.bzero;
    dstype *rootpqarray = sna.rootpqarray;
    dstype *cglist = sna.cglist;
    dstype *rcutsq = sna.rcutsq;    
    dstype *radelem = sna.radelem;
    dstype *wjelem = sna.wjelem; 
    dstype *coeffelem = sna.coeffelem;           // element bispectrum coefficients
    //dstype *fac = sna.fac;
    
    int *atomtype = nb.atomtype;       // type of each atom i      
    int *alist = nb.alist;
    int *neighnum = nb.neighnum;  // numbers of neighbors for each atom i 
    int *neighlist = nb.neighlist; // list of neighbors for each atom i    
    int *mask = nb.mask;        
  
    int na = inum;
    int ne = 0; 
    int *pairnum = &tmp.intmem[ne]; // na
    int *pairlist = &tmp.intmem[na]; // na*neighmax
    ne = na + na*neighmax;
    int *pairnumsum = &tmp.intmem[ne]; // na+1         
    
    dstype rcutsq1 = ArrayGetValueAtIndex(rcutsq, 8, backend);
    
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutsq1, alist, neighlist, 
        neighnum, inum, neighmax, dim, backend);
    
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    
    ne += na+1;
    int *ai = &tmp.intmem[ne]; // ineighmax        
    int *aj = &tmp.intmem[ne+ineighmax]; // ineighmax        
    int *ti = &tmp.intmem[ne+2*ineighmax]; // ineighmax        
    int *tj = &tmp.intmem[ne+3*ineighmax]; // ineighmax        
    int *aii = &tmp.intmem[ne+4*ineighmax]; // ineighmax    
    // 1 + 2*na + na*neighmax + 5*ineighmax
    
    ne = 0;
    dstype *rij = &tmp.tmpmem[ne]; // ineighmax*dim    
    ne += ineighmax*dim; 
    dstype *ulist_r = &tmp.tmpmem[ne]; 
    dstype *zlist_r = &tmp.tmpmem[ne]; 
    dstype *ylist_r = &tmp.tmpmem[ne]; 
    ne += max(max(idxu_max*ineighmax, idxz_max*ndoubles*na), idxu_max*nelements*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    dstype *ylist_i = &tmp.tmpmem[ne]; 
    ne += max(max(idxu_max*ineighmax, idxz_max*ndoubles*na), idxu_max*nelements*na); 
    dstype *dulist_r = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;
    dstype *dulist_i = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;    
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    dstype *dedr = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;    
    dstype *beta = &tmp.tmpmem[ne];
    ne += ncoeff*na;    
    dstype *blist = &tmp.tmpmem[ne]; // dxb_max*ntriples*na          
              
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
   
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
    
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
          inum, switchflag, chemflag, backend);
     
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
      ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
    
    ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);

    SnapTallyEnergyFull(eatom, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);
    
    ComputeBeta2(beta, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);
    
    ComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, beta, idxz, idxb_block, 
          idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bnormflag, ncoeff, inum, backend);
        
    ComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, atomtype, map, 
          ai, aj, nelem, twojmax, idxu_max, chemflag, ineighmax, backend);
    
    SnapTallyForceFull(fatom, dedr, ai, aj, ineighmax, backend);

    SnapTallyVirialFull(vatom, dedr, rij, ai, aj, ineighmax, backend);  
}

void ComputeMDPSna(dstype *snaarray, snastruct &sna, commonstruct &common, 
        sysstruct &sys, neighborstruct &nb, tempstruct &tmp)
{    
    int inum = common.inum;
    //int anum = common.anum;
    int neighmax = common.neighmax;    
    int dim = common.dim;
    int backend = common.backend;
    
    
    int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int nperdim = sna.nperdim;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int quadraticflag = sna.quadraticflag;
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    dstype wself = sna.wself;
    dstype rmin0 = sna.rmin0;
    dstype rfac0 = sna.rfac0;
    dstype rcutfac = sna.rcutfac;
    dstype rcutmax = sna.rcutmax;        
    dstype *bzero = sna.bzero;
    dstype *rootpqarray = sna.rootpqarray;
    dstype *cglist = sna.cglist;
    dstype *rcutsq = sna.rcutsq;    
    dstype *radelem = sna.radelem;
    dstype *wjelem = sna.wjelem; 
    
    int *atomtype = nb.atomtype;       // type of each atom i      
    int *alist = nb.alist;
    int *neighnum = nb.neighnum;  // numbers of neighbors for each atom i 
    int *neighlist = nb.neighlist; // list of neighbors for each atom i    
    int *mask = nb.mask;        
  
    int na = inum;
    int ne = 0; 
    int *pairnum = &tmp.intmem[ne]; // na
    int *pairlist = &tmp.intmem[na]; // na*neighmax
    ne = na + na*neighmax;
    int *pairnumsum = &tmp.intmem[ne]; // na+1         
    
    dstype rcutsq1 = ArrayGetValueAtIndex(rcutsq, 8, backend);
    
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutsq1, alist, neighlist, 
        neighnum, inum, neighmax, dim, backend);
    
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    
    ne += na+1;
    int *ai = &tmp.intmem[ne]; // ineighmax        
    int *aj = &tmp.intmem[ne+ineighmax]; // ineighmax        
    int *ti = &tmp.intmem[ne+2*ineighmax]; // ineighmax        
    int *tj = &tmp.intmem[ne+3*ineighmax]; // ineighmax        
    int *aii = &tmp.intmem[ne+4*ineighmax]; // ineighmax    
    // 1 + 2*na + na*neighmax + 5*ineighmax
    
    ne = 0;
    dstype *rij = &tmp.tmpmem[ne]; // ineighmax*dim    
    ne += ineighmax*dim; 
    dstype *ulist_r = &tmp.tmpmem[ne]; 
    dstype *zlist_r = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na); 
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;        
    dstype *blist = &tmp.tmpmem[ne]; // idxb_max*ntriples*na          
              
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
             
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
    
    ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
    
    ComputeSna(snaarray, blist, alist, mask,  ncoeff, nperdim, inum, quadraticflag, backend);    
}

void ComputeMDPSnad(dstype *snad, snastruct &sna, commonstruct &common, 
        sysstruct &sys, neighborstruct &nb, tempstruct &tmp)
{    
    int inum = common.inum;
    //int anum = common.anum;
    int neighmax = common.neighmax;    
    int dim = common.dim;
    int backend = common.backend;
    
    
    int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int nperdim = sna.nperdim;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int quadraticflag = sna.quadraticflag;
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    dstype wself = sna.wself;
    dstype rmin0 = sna.rmin0;
    dstype rfac0 = sna.rfac0;
    dstype rcutfac = sna.rcutfac;
    dstype rcutmax = sna.rcutmax;        
    dstype *bzero = sna.bzero;
    dstype *rootpqarray = sna.rootpqarray;
    dstype *cglist = sna.cglist;
    dstype *rcutsq = sna.rcutsq;    
    dstype *radelem = sna.radelem;
    dstype *wjelem = sna.wjelem; 
    
    int *atomtype = nb.atomtype;       // type of each atom i      
    int *alist = nb.alist;
    int *neighnum = nb.neighnum;  // numbers of neighbors for each atom i 
    int *neighlist = nb.neighlist; // list of neighbors for each atom i    
    int *mask = nb.mask;        
  
    int na = inum;
    int ne = 0; 
    int *pairnum = &tmp.intmem[ne]; // na
    int *pairlist = &tmp.intmem[na]; // na*neighmax
    ne = na + na*neighmax;
    int *pairnumsum = &tmp.intmem[ne]; // na+1         
    
    dstype rcutsq1 = ArrayGetValueAtIndex(rcutsq, 8, backend);
    
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutsq1, alist, neighlist, 
        neighnum, inum, neighmax, dim, backend);
    
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    
    ne += na+1;
    int *ai = &tmp.intmem[ne]; // ineighmax        
    int *aj = &tmp.intmem[ne+ineighmax]; // ineighmax        
    int *ti = &tmp.intmem[ne+2*ineighmax]; // ineighmax        
    int *tj = &tmp.intmem[ne+3*ineighmax]; // ineighmax        
    int *aii = &tmp.intmem[ne+4*ineighmax]; // ineighmax    
    // 1 + 2*na + na*neighmax + 5*ineighmax
    
    ne = 0;
    dstype *rij = &tmp.tmpmem[ne]; // ineighmax*dim    
    ne += ineighmax*dim; 
    dstype *ulist_r = &tmp.tmpmem[ne]; 
    dstype *zlist_r = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na);     
    dstype *dulist_r = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;
    dstype *dulist_i = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;    
    dstype *blist; // idxb_max*ntriples*inum;
    if (quadraticflag) {
        blist = &tmp.tmpmem[ne];
        ne += idxb_max*ntriples*inum;
    }        
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    dstype *dblist = &tmp.tmpmem[ne]; // idxb_max*ntriples*dim*ineighmax                        
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;            
        
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
      ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    if (quadraticflag) {
        ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
              map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
    }
          
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
        idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);
    
}

void ComputeMDPSnav(dstype *snav, snastruct &sna, commonstruct &common, 
        sysstruct &sys, neighborstruct &nb, tempstruct &tmp)
{    
    int inum = common.inum;
    //int anum = common.anum;
    int neighmax = common.neighmax;    
    int dim = common.dim;
    int backend = common.backend;
        
    int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int nperdim = sna.nperdim;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int quadraticflag = sna.quadraticflag;
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    dstype wself = sna.wself;
    dstype rmin0 = sna.rmin0;
    dstype rfac0 = sna.rfac0;
    dstype rcutfac = sna.rcutfac;
    dstype rcutmax = sna.rcutmax;        
    dstype *bzero = sna.bzero;
    dstype *rootpqarray = sna.rootpqarray;
    dstype *cglist = sna.cglist;
    dstype *rcutsq = sna.rcutsq;    
    dstype *radelem = sna.radelem;
    dstype *wjelem = sna.wjelem; 
    
    int *atomtype = nb.atomtype;       // type of each atom i      
    int *alist = nb.alist;
    int *neighnum = nb.neighnum;  // numbers of neighbors for each atom i 
    int *neighlist = nb.neighlist; // list of neighbors for each atom i    
    int *mask = nb.mask;        
  
    int na = inum;
    int ne = 0; 
    int *pairnum = &tmp.intmem[ne]; // na
    int *pairlist = &tmp.intmem[na]; // na*neighmax
    ne = na + na*neighmax;
    int *pairnumsum = &tmp.intmem[ne]; // na+1         
    
    dstype rcutsq1 = ArrayGetValueAtIndex(rcutsq, 8, backend);
    
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutsq1, alist, neighlist, 
        neighnum, inum, neighmax, dim, backend);
    
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    
    ne += na+1;
    int *ai = &tmp.intmem[ne]; // ineighmax        
    int *aj = &tmp.intmem[ne+ineighmax]; // ineighmax        
    int *ti = &tmp.intmem[ne+2*ineighmax]; // ineighmax        
    int *tj = &tmp.intmem[ne+3*ineighmax]; // ineighmax        
    int *aii = &tmp.intmem[ne+4*ineighmax]; // ineighmax    
    // 1 + 2*na + na*neighmax + 5*ineighmax
    
    ne = 0;
    dstype *rij = &tmp.tmpmem[ne]; // ineighmax*dim    
    ne += ineighmax*dim; 
    dstype *ulist_r = &tmp.tmpmem[ne]; 
    dstype *zlist_r = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na);     
    dstype *dulist_r = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;
    dstype *dulist_i = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;    
    dstype *blist; // idxb_max*ntriples*inum;
    if (quadraticflag) {
        blist = &tmp.tmpmem[ne];
        ne += idxb_max*ntriples*inum;
    }        
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    dstype *dblist = &tmp.tmpmem[ne]; // idxb_max*ntriples*dim*ineighmax                        
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;            
        
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
      ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    if (quadraticflag) {
        ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
              map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
    }
          
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
        idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);
        
}

void ComputeMDPSnapArrays(dstype *snaarray, dstype *snad, dstype *snav, dstype *dbdx, snastruct &sna, 
        commonstruct &common, sysstruct &sys, neighborstruct &nb, tempstruct &tmp)
{    
    int inum = common.inum;
    //int anum = common.anum;
    int neighmax = common.neighmax;    
    int dim = common.dim;
    int backend = common.backend;
        
    int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int nperdim = sna.nperdim;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int quadraticflag = sna.quadraticflag;
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    dstype wself = sna.wself;
    dstype rmin0 = sna.rmin0;
    dstype rfac0 = sna.rfac0;
    dstype rcutfac = sna.rcutfac;
    dstype rcutmax = sna.rcutmax;        
    dstype *bzero = sna.bzero;
    dstype *rootpqarray = sna.rootpqarray;
    dstype *cglist = sna.cglist;
    dstype *rcutsq = sna.rcutsq;    
    dstype *radelem = sna.radelem;
    dstype *wjelem = sna.wjelem; 
    
    int *atomtype = nb.atomtype;       // type of each atom i      
    int *alist = nb.alist;
    int *neighnum = nb.neighnum;  // numbers of neighbors for each atom i 
    int *neighlist = nb.neighlist; // list of neighbors for each atom i    
    int *mask = nb.mask;        
    int *tag = nb.tag;        
  
    int na = inum;
    int ne = 0; 
    int *pairnum = &tmp.intmem[ne]; // na
    int *pairlist = &tmp.intmem[na]; // na*neighmax
    ne = na + na*neighmax;
    int *pairnumsum = &tmp.intmem[ne]; // na+1         
    
    dstype rcutsq1 = ArrayGetValueAtIndex(rcutsq, 8, backend);
    
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutsq1, alist, neighlist, 
        neighnum, inum, neighmax, dim, backend);
    
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    
    ne += na+1;
    int *ai = &tmp.intmem[ne]; // ineighmax        
    int *aj = &tmp.intmem[ne+ineighmax]; // ineighmax        
    int *ti = &tmp.intmem[ne+2*ineighmax]; // ineighmax        
    int *tj = &tmp.intmem[ne+3*ineighmax]; // ineighmax        
    int *aii = &tmp.intmem[ne+4*ineighmax]; // ineighmax    
    // 1 + 2*na + na*neighmax + 5*ineighmax
    
    ne = 0;
    dstype *rij = &tmp.tmpmem[ne]; // ineighmax*dim    
    ne += ineighmax*dim; 
    dstype *ulist_r = &tmp.tmpmem[ne]; 
    dstype *zlist_r = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    ne += max(idxu_max*ineighmax, idxz_max*ndoubles*na);     
    dstype *dulist_r = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;
    dstype *dulist_i = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;    
    dstype *blist = &tmp.tmpmem[ne]; // idxb_max*ntriples*inum;    
    ne += idxb_max*ntriples*inum;    
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    dstype *dblist = &tmp.tmpmem[ne]; // idxb_max*ntriples*dim*ineighmax                        
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;            
        
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
      ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
            nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
          
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
        idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeSna(snaarray, blist, alist, mask, atomtype, ncoeff, ntypes, nperdim, inum, quadraticflag, backend);    
    
    ComputeSnad(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);

    ComputeSnad(dbdx, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);
    
    ComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);        
}


#endif        
 

