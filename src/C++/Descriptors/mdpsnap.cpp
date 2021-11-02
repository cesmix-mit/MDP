/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDPSNAP
#define MDPSNAP

void printui(dstype *Ui, int inum, int i, int twojmax)
{
    int k = 0;
    for (int j = 0; j <= twojmax; j++) {
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                int n = i + inum*k;   
                cout << scientific << Ui[n] << "   ";
                k += 1;
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}
        
void setui(dstype *Ui, int inum, int twojmax, int idxu_max)
{    
    for (int d=0; d<3; d++)
    for (int i = 0; i < inum; i++) {
        int k = 0;
        for (int j = 0; j <= twojmax; j++) {
            for (int mb = 0; mb <= j; mb++) {
                for (int ma = 0; ma <= j; ma++) {
                    if (j%2==1 && 2*mb==j+1) {
                        int n = i + inum*k + inum*idxu_max*d;   
                        Ui[n] = 0.0;
                    }
                    k += 1;
                }
            }
        }
    }
}

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
        
        sna.idxcg_max = sna.idx_max[0];
        sna.idxu_max = sna.idx_max[1];
        sna.idxb_max = sna.idx_max[2];
        sna.idxz_max = sna.idx_max[3];        
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
        
        sna.idxcg_max = idx_max[0];
        sna.idxu_max = idx_max[1];
        sna.idxb_max = idx_max[2];
        sna.idxz_max = idx_max[3];
        
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

void InitSnap(snastruct &sna, commonstruct &common)
{
  int backend = common.backend;
  int dim = common.dim;
  int inum = common.inum;
  int nelements = (int) common.snaparam[0];  
  int ncoeffall = (int) common.snaparam[1];  
  int twojmax = (int) common.snaparam[2];  
  dstype rcutfac = common.snaparam[3];
  dstype rfac0 = common.snaparam[4];
  dstype rmin0 = common.snaparam[5];
  int bzeroflag = (int) common.snaparam[6];
  int switchflag = (int) common.snaparam[7];
  int quadraticflag = (int) common.snaparam[8];
  int chemflag = (int) common.snaparam[9];
  int bnormflag = (int) common.snaparam[10];
  int wselfallflag = (int) common.snaparam[11];    
  dstype wself=1.0;  
    
  int ncoeff, ncoeffq;
  if (!quadraticflag)
    ncoeff = ncoeffall - 1;
  else {
    // ncoeffall should be (ncoeff+2)*(ncoeff+1)/2 = 1 + ncoeff + ncoeff*(ncoeff+1)/2
    // so, ncoeff = floor(sqrt(2*ncoeffall))-1      
    ncoeff = sqrt(2*ncoeffall)-1;
    ncoeffq = (ncoeff*(ncoeff+1))/2;
    int ntmp = 1+ncoeff+ncoeffq;
    if (ntmp != ncoeffall) {
      error("Incorrect SNAP coeff file");
    }
  }

  // Calculate maximum cutoff for all elements
  dstype rcutmax = 0.0;
  for (int ielem = 0; ielem < nelements; ielem++)
    rcutmax = max(2.0*common.snapelemradius[ielem]*rcutfac,rcutmax);
      
  SetSnap(sna, twojmax, nelements, backend);  
  TemplateCopytoDevice(&sna.radelem[1], common.snapelemradius, nelements, backend);
  TemplateCopytoDevice(&sna.wjelem[1], common.snapelemweight, nelements, backend);  
  ArrayFill(&sna.map[1], (int) 0, nelements, backend);
  
  dstype cutsq[100];
  for (int i=0; i<nelements; i++)
      for (int j=0; j<nelements; j++) {
          dstype cut = (common.snapelemradius[i] + common.snapelemradius[j])*rcutfac;
          cutsq[j+1 + (i+1)*(nelements+1)] = cut*cut;
      }
  TemplateCopytoDevice(sna.rcutsq, cutsq, (nelements+1)*(nelements+1), backend);  
  
  if (bzeroflag) {
    dstype www = wself*wself*wself;
    dstype bzero[100];
    for (int j = 0; j <= twojmax; j++)
      if (bnormflag)
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
    TemplateCopytoDevice(sna.bzero, bzero, twojmax+1, backend);
  }
  TemplateMalloc(&sna.coeffelem, nelements*ncoeffall, backend);        
  TemplateCopytoDevice(sna.coeffelem, common.snapcoeff, nelements*ncoeffall, backend);        
  
  if (!chemflag)    
    nelements = 1;
  
  sna.ncoeff = ncoeff;
  sna.ncoeffall = ncoeffall;
  sna.nperdim = ncoeffall-1;
  sna.nelements = nelements;    
  sna.ndoubles = nelements*nelements;   // number of multi-element pairs
  sna.ntriples = nelements*nelements*nelements;   // number of multi-element triplets      
  sna.beta_max = inum;                 // length of beta
  sna.bnormflag = bnormflag;
  sna.chemflag = chemflag;    
  sna.quadraticflag = quadraticflag;
  sna.switchflag = switchflag;
  sna.bzeroflag = bzeroflag;
  sna.wselfallflag = wselfallflag;
  sna.wself = wself;
  sna.rmin0 = rmin0;
  sna.rfac0 = rfac0;
  sna.rcutfac = rcutfac;
  sna.rcutmax = rcutmax;      
}

void ComputePairSnap(dstype *eatom, dstype *fatom, dstype *vatom, dstype *x, snastruct &sna, commonstruct &common, 
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
        
    INIT_TIMING;
    
    START_TIMING;
    //dstype *x = &sys.x[0];
//     NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
//         neighnum, inum, neighmax, dim, backend);
    NeighborPairList(pairnum, pairlist, x, rcutsq, alist, neighlist, 
                            neighnum, atomtype, alist, inum, neighmax, dim, ntypes, backend);    
    END_TIMING(50);
    
    START_TIMING;
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    END_TIMING(51);
    
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
              
    START_TIMING;
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
    END_TIMING(52);
        
    START_TIMING;    
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
    END_TIMING(53);
    
    START_TIMING;
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);
    END_TIMING(54);
        
//     cout<<inum<<"  "<<idxu_max<<endl;
//     dstype *Stot_r, *Stot_i;
//     TemplateMalloc(&Stot_r, inum*idxu_max, backend);  
//     TemplateMalloc(&Stot_i, inum*idxu_max, backend);  
//     ZeroUarraytot(Stot_r, Stot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);        
//     ArrayAXPBY(Stot_r, Stot_r, ulisttot_r, one, minusone, inum*idxu_max, backend);    
//     dstype norme = PNORM(common.cublasHandle, inum*idxu_max, Stot_r, backend);
//     cout<<"Error : "<<norme<<endl;        
    
    START_TIMING;
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
          inum, switchflag, chemflag, backend);
//     AddUarraytot2(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
//      idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
//           inum,  ineighmax, switchflag, chemflag, backend);
    END_TIMING(55);
        
//     SnapComputeUi(Stot_r, Stot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, map, ai, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);        
//     ArrayAXPBY(Stot_r, Stot_r, ulisttot_r, one, minusone, inum*idxu_max, backend);    
//     norme = PNORM(common.cublasHandle, inum*idxu_max, Stot_r, backend);
//     cout<<"Error : "<<norme<<endl;    
//     error("here");
    
//     ComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
//         T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag, int backend)
        
//     writearray2file("x1.bin", x, 3*common.anum, backend);        
//     writearray2file("rij1.bin", rij, 3*ineighmax, backend);        
//     writearray2file("ulistr1.bin", ulist_r, idxu_max*ineighmax, backend);        
//     writearray2file("ulisttotr1.bin", ulisttot_r, idxu_max*nelem*inum, backend);  
    
    //printArray2D(ulisttot_r, idxu_max, inum, backend);
        
    START_TIMING;
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
      ai, aj, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    END_TIMING(56);
    
//    writearray2file("dulistr1.bin", dulist_r, 3*idxu_max*ineighmax, backend);        
    
    START_TIMING;
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
    END_TIMING(57);
    
//    writearray2file("zlistr1.bin", zlist_r, idxz_max*ndoubles*inum, backend);     
    
    START_TIMING;
    ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);
    END_TIMING(58);
    
//    writearray2file("blist1.bin", blist, idxb_max*ntriples*inum, backend);     
    
//     dstype *dblist;  TemplateMalloc(&dblist, ineighmax*3*(ncoeffall-1), backend);  
//     ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, ai, aj, ti, tj, 
//         twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    START_TIMING;
    SnapTallyEnergyFull(eatom, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);    
    END_TIMING(59);
//     printArray2D(map, 1, nelements+1, backend);
//     printArray2D(coeffelem, 1, ncoeff+1, backend);
//     printArray2D(blist, ncoeff, inum, backend);
//     printArray2D(eatom, 1, inum, backend);
    
//     writearray2file("eatom1.bin", eatom, inum, backend);    
//     writearray2file("coeffelem1.bin", coeffelem, ntypes*ncoeffall, backend);    
    
    START_TIMING;
    ComputeBeta2(beta, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);
    END_TIMING(60);
    
//    writearray2file("beta1.bin", beta, ncoeff*inum, backend);    
    
    START_TIMING;
    ComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, beta, idxz, idxb_block, 
          idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bnormflag, ncoeff, inum, backend);
    END_TIMING(61);
    
//    writearray2file("ylistr1.bin", ylist_r, idxu_max*inum*nelem, backend);    
    
    START_TIMING;
    ComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, map, 
          ai, aj, ti, tj, nelem, twojmax, idxu_max, chemflag, ineighmax, backend);        
    END_TIMING(62);
    
//    writearray2file("dedr1.bin", dedr, 3*ineighmax, backend);    
    
    START_TIMING;
    SnapTallyForceFull(fatom, dedr, ai, aj, alist, ineighmax, backend);
    END_TIMING(63);
    
    START_TIMING;
    SnapTallyVirialFull(vatom, dedr, rij, ai, aj, inum, ineighmax, backend);  
    END_TIMING(64);
        
//     printf("%i %i %i %i %i %i %i %i\n", common.anum, inum, ineighmax, idxu_max, idxb_max, idxz_max, bnormflag, ncoeff);
//     error("here");
    
//     dstype *bi, *bd, *bv, *e, *f, *ri1, *bi2, *blist2;
//     int M = ncoeffall-1;
//     int nperdim = M;
//     TemplateMalloc(&bi, M, backend);  
//     TemplateMalloc(&bd, ineighmax*3*M, backend);  
//     TemplateMalloc(&bv, 6*M, backend);  
//     TemplateMalloc(&e, inum, backend);  
//     TemplateMalloc(&f, inum*3, backend);  
//     TemplateMalloc(&ri1, inum*3, backend);  
//     TemplateMalloc(&bi2, M, backend);      
//     TemplateMalloc(&blist2, inum*M, backend);  
//         
//     for (int i=0; i<M; i++) {
//         bi[i] = 0.0;
//         for (int j=0; j<inum; j++)
//             bi[i] += blist[i + M*j];
//     }    
// 
//     dstype epsil = 1e-3;
//     for (int i = 0; i<dim*ineighmax; i++) {
//         cpuArrayCopy(ri1, rij, dim*ineighmax);  
//         ri1[i] = ri1[i] + epsil;
//         int d = i%3;
//         int m = (i-d)/3;
//                 
//         ComputeUij(ulist_r, ulist_i, rootpqarray, ri1, radelem, rmin0, rfac0, rcutfac, 
//               idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);   
// 
//         ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
//                 map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);
// 
//         AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, ri1, wjelem, radelem, rmin0, rcutfac, 
//          idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
//               inum, switchflag, chemflag, backend);
//         
//         ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
//               idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
// 
//         ComputeBi(blist2, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
//               map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
//               nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);
//             
//         for (int ia=0; ia<M; ia++) {
//             bi2[ia] = 0.0;
//             for (int j=0; j<inum; j++)
//                 bi2[ia] += blist2[ia + M*j];
//         }    
//         
//         printArray2D(bi, 1, M, common.backend);
//         printArray2D(bi2, 1, M, common.backend);
//         for (int j=0; j<M; j++) 
//             cout<<(bi2[j]-bi[j])/epsil<<" ";
//         cout<<endl;
//         for (int j=0; j<M; j++) 
//             cout<<(blist2[j+M*m]-blist[j+M*m])/epsil<<" ";
//         cout<<endl;
//         for (int j=0; j<M; j++) 
//             cout<<dblist[d + 3*j + 3*M*m]<<" "; // dim*M*ijnum
//         cout<<endl;
//         
//         error("here");
//     }
//     
        
// #ifdef USE_CUDA    
// //     START_TIMING;    
// //     int n = idxu_max*ineighmax;
// //     gpuComputeSij(ulist_r, ulist_i, &dulist_r[0*n], &dulist_i[0*n], &dulist_r[1*n], &dulist_i[1*n],
// //             &dulist_r[2*n], &dulist_i[2*n], rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
// //           idxu_block, ti, tj, twojmax, idxu_max, ineighmax, switchflag);        
// //     END_TIMING(65);
// //     
// //     START_TIMING;
// //     AddUarraytot2(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
// //      idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
// //           inum,  ineighmax, switchflag, chemflag, backend);        
// //     END_TIMING(66);    
// //         
// //     START_TIMING;
// //     gpuComputeZi2(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
// //           idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum);
// //     END_TIMING(67);    
// #else
//     int ijnum = ineighmax;
//     int m = idxu_max*inum;    
//     int n = idxu_max*ineighmax;
//     int z = idxz_max*inum;    
//     int b = idxb_max*inum;        
//     
//     dstype *tm, *S_r, *S_i, *dS_r, *dS_i, *Stot_r, *Stot_i, *Z_r, *Z_i, *bl, *be;
//     dstype *Y_r, *Y_i, *dE, *ea, *fa, *va;
//     TemplateMalloc(&tm, n*3, backend);    
//     TemplateMalloc(&S_r, n, backend);    
//     TemplateMalloc(&S_i, n, backend);    
//     TemplateMalloc(&dS_r, n*3, backend);    
//     TemplateMalloc(&dS_i, n*3, backend);    
//     TemplateMalloc(&Stot_r, m, backend);    
//     TemplateMalloc(&Stot_i, m, backend);    
//     TemplateMalloc(&Z_r, z*ndoubles, backend);    
//     TemplateMalloc(&Z_i, z*ndoubles, backend);    
//     TemplateMalloc(&bl, b*ntriples, backend);        
//     TemplateMalloc(&be, ncoeff*inum, backend);                    
//     TemplateMalloc(&Y_r, m*nelements, backend);    
//     TemplateMalloc(&Y_i, m*nelements, backend);    
//     TemplateMalloc(&dE, ijnum*3, backend);                    
//     TemplateMalloc(&ea, inum, backend);                    
//     TemplateMalloc(&fa, inum*3, backend);                    
//     TemplateMalloc(&va, inum*6, backend);                    
//     
//     ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
//     
//     ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);
//         
//     AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
//      idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelem, 
//           inum, switchflag, chemflag, backend);
//     
//     ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
//       rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
//       ai, aj, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);     
//             
//     //ArraySetValue(S_r, 100.0, n, backend);
//     //ArraySetValue(S_i, 100.0, n, backend);
//     ComputeSij(S_r, S_i, &dS_r[0*n], &dS_i[0*n], &dS_r[1*n], &dS_i[1*n],
//             &dS_r[2*n], &dS_i[2*n], rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//             idxu_block, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);      
//             
//     ZeroUarraytot2(Stot_r, Stot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);
// 
//     AddUarraytot(Stot_r, Stot_i, S_r, S_i, pairnum, pairnumsum, 
//             map, tj, idxu_max, nelem, inum, ijnum, chemflag, backend);
//     
//     cpuArrayTranspose(tm, Stot_r, inum, idxu_max);
//     printf("Max error = %g\n",cpuArrayMaxError(ulisttot_r, tm, m));
//     cpuArrayTranspose(tm, Stot_i, inum, idxu_max);
//     printf("Max error = %g\n",cpuArrayMaxError(ulisttot_i, tm, m));
//     
//     //  ijnum*idxu_max*3 -> 3*idxu_max*ijnum (permute 13)
//     cpuPermute13(tm, dS_r, ijnum, idxu_max, 3, 1);
//     printf("Max error = %g\n",cpuArrayMaxError(dulist_r, tm, n*3));
//     cpuPermute13(tm, dS_i, ijnum, idxu_max, 3, 1);
//     printf("Max error = %g\n",cpuArrayMaxError(dulist_i, tm, n*3));
//     
//     ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
//       idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);    
//     ComputeZi2(Z_r, Z_i, Stot_r, Stot_i, cglist, idxz, idxu_block, 
//         idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);            
//     cpuArrayTranspose(tm, Z_r, inum, idxz_max*ndoubles);
//     printf("Max error = %g\n",cpuArrayMaxError(zlist_r, tm, inum*idxz_max*ndoubles));
//     cpuArrayTranspose(tm, Z_i, inum, idxz_max*ndoubles);
//     printf("Max error = %g\n",cpuArrayMaxError(zlist_i, tm, inum*idxz_max*ndoubles));
//     
//     ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
//       map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
//       nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);
//     ComputeBi2(bl, Z_r, Z_i, Stot_r, Stot_i, bzero, alist, atomtype, 
//       map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
//       nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);    
//     cpuArrayTranspose(tm, bl, inum, idxb_max*ntriples);
//     printf("Max error = %g\n",cpuArrayMaxError(blist, tm, inum*idxb_max*ntriples));
// 
// //     ArraySetValue(eatom, 0.0, inum, backend);
// //     ArraySetValue(ea, 0.0, inum, backend);
// //     SnapTallyEnergyFull(eatom, blist, coeffelem, alist, map, atomtype, 
// //         inum, ncoeff, ncoeffall, quadraticflag, backend);        
// //     cpuSnapTallyEnergyFull2(ea, bl, coeffelem, alist, map, atomtype, 
// //         inum, ncoeff, ncoeffall, quadraticflag);    
// //     printf("Max error = %g\n",cpuArrayMaxError(eatom, ea, inum));
//         
//     ComputeBeta2(beta, blist, coeffelem, alist, map, atomtype, 
//         inum, ncoeff, ncoeffall, quadraticflag, backend);
//     ComputeBeta(be, bl, coeffelem, alist, map, atomtype, 
//         inum, ncoeff, ncoeffall, quadraticflag, backend);
//     cpuArrayTranspose(tm, be, inum, ncoeff);
//     printf("Max error = %g\n",cpuArrayMaxError(beta, tm, inum*ncoeff));    
//     
//     ComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, beta, idxz, idxb_block, 
//           idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, bnormflag, ncoeff, inum, backend);    
//     ComputeYi(Y_r, Y_i, Stot_r, Stot_i, cglist, be, idxz, idxb_block, 
//           idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, bnormflag, inum, backend);
//     cpuArrayTranspose(tm, Y_r, inum, idxu_max*nelements);
//     printf("Max error = %g\n",cpuArrayMaxError(ylist_r, tm, inum*idxu_max*nelements));
//     cpuArrayTranspose(tm, Y_i, inum, idxu_max*nelements);
//     printf("Max error = %g\n",cpuArrayMaxError(ylist_i, tm, inum*idxu_max*nelements));
// 
//     ComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, map, 
//           ai, aj, ti, tj, nelem, twojmax, idxu_max, chemflag, ineighmax, backend);            
//     ComputeDeidrj(dE, Y_r, Y_i, dS_r, dS_i, idxu_block, map, 
//           ai, tj, twojmax, idxu_max, chemflag, inum, ijnum, backend);        
//     printf("Max error = %g\n",cpuArrayMaxError(dedr, dE, ijnum*3));
//     
// //     ArraySetValue(fatom, 0.0, inum*3, backend);
// //     ArraySetValue(fa, 0.0, inum*3, backend);
// //     SnapTallyForceFull(fatom, dedr, ai, aj, alist, ineighmax, backend);
// //     SnapTallyForceFull(fa, dE, ai, aj, alist, ineighmax, backend);
// //     printf("Max error = %g\n",cpuArrayMaxError(fatom, fa, inum*3));
//     
// //     ArraySetValue(vatom, 0.0, inum*6, backend);
// //     ArraySetValue(va, 0.0, inum*6, backend);
// //     SnapTallyVirialFull(vatom, dedr, rij, ai, aj, inum, ineighmax, backend);  
// //     SnapTallyVirialFull(va, dE, rij, ai, aj, inum, ineighmax, backend);  
// //     printf("Max error = %g\n",cpuArrayMaxError(vatom, va, inum*6));    
//     
//     TemplateFree(tm, backend);    
//     TemplateFree(S_r, backend);    
//     TemplateFree(S_i, backend);    
//     TemplateFree(dS_r, backend);    
//     TemplateFree(dS_i, backend);    
//     TemplateFree(Stot_r, backend);    
//     TemplateFree(Stot_i, backend);    
//     TemplateFree(Z_r, backend);    
//     TemplateFree(Z_i, backend);    
//     TemplateFree(Y_r, backend);    
//     TemplateFree(Y_i, backend);    
//     TemplateFree(bl, backend);    
//     TemplateFree(be, backend);    
//     TemplateFree(dE, backend);    
//     TemplateFree(ea, backend);    
//     TemplateFree(fa, backend);    
//     TemplateFree(va, backend);            
// #endif                 
}

void ComputePairSnap2(dstype *eatom, dstype *fatom, dstype *vatom, dstype *x, snastruct &sna, commonstruct &common, 
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
        
    INIT_TIMING;
    
    START_TIMING;
    //dstype *x = &sys.x[0];
//     NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
//         neighnum, inum, neighmax, dim, backend);
    NeighborPairList(pairnum, pairlist, x, rcutsq, alist, neighlist, 
                            neighnum, atomtype, alist, inum, neighmax, dim, ntypes, backend);
    END_TIMING(50);
            
//    writearray2file("x.bin", x, 3*common.anum, backend);        
//     //writearray2file("pairlist.bin", x, 3*N, backend);        
//     //writearray2file("pairnum.bin", pairnum, inum, backend);        
//     printArray2D(pairlist, neighmax, inum, backend);
//     printArray2D(pairnum, 1, inum, backend);
//     printArray2D(neighnum, 1, inum, backend);
//     error("here");
    
//     printf("%g ", rcutmax);
//     printArray2D(rcutsq, 1, 9, backend);
//     error("here");            
    START_TIMING;
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    END_TIMING(51);
    
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
    dstype *dedr = &tmp.tmpmem[ne];    
    //dstype *ylist_r = &tmp.tmpmem[ne]; 
    ne += max(max(idxu_max*ineighmax, idxz_max*ndoubles*na), idxu_max*nelements*na); 
    dstype *ulist_i = &tmp.tmpmem[ne]; 
    dstype *zlist_i = &tmp.tmpmem[ne]; 
    //dstype *ylist_i = &tmp.tmpmem[ne]; 
    ne += max(max(idxu_max*ineighmax, idxz_max*ndoubles*na), idxu_max*nelements*na); 
    dstype *dulist_r = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;
    dstype *dulist_i = &tmp.tmpmem[ne];
    ne += idxu_max*dim*ineighmax;    
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    dstype *ylist_r = &tmp.tmpmem[ne]; 
    //dstype *dedr = &tmp.tmpmem[ne];    
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    dstype *ylist_i = &tmp.tmpmem[ne]; 
    ne += idxu_max*nelements*na;    
    dstype *beta = &tmp.tmpmem[ne];
    //dstype *dedr = &tmp.tmpmem[ne];    
    ne += ncoeff*na;    
    dstype *blist = &tmp.tmpmem[ne]; // dxb_max*ntriples*na          
              
    START_TIMING;
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
    END_TIMING(52);
        
    START_TIMING;    
    int n = idxu_max*ineighmax;
    ComputeSij(ulist_r, ulist_i, &dulist_r[0*n], &dulist_i[0*n], &dulist_r[1*n], &dulist_i[1*n],
            &dulist_r[2*n], &dulist_i[2*n], rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
            idxu_block, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);      
    END_TIMING(53);
        
//     writearray2file("x2.bin", x, 3*common.anum, backend);        
//     writearray2file("rij2.bin", rij, 3*ineighmax, backend);        
//     writearray2file("ulistr2.bin", ulist_r, idxu_max*ineighmax, backend);  
    //printArray2D(ulist_r, ineighmax, 10, backend);  
    
//     for (int i=0; i<pairnum[0]; i++) {
//         for (int j=0; j<3; j++)
//             printf("%g ", rij[j+3*(i+pairnum[0])]);
//         printf("\n");
//     }
//     for (int i=0; i<pairnum[0]; i++) {
//         for (int j=0; j<12; j++)
//             printf("%g ", ulist_r[(i+pairnum[0]) + ineighmax*j]);
//         printf("\n");
//     }
    
    //printArray2D(ulist_r, ineighmax, 12, backend);
//     printArray2D(ulist_r, ineighmax, 10, backend);
//     printArray2D(wjelem, 1, ntypes+1, backend);
//     printArray2D(radelem, 1, ntypes+1, backend);
//     printArray2D(rij, 3, pairnum[0], backend);
//     printArray2D(ti, 1, inum, backend);
//     printArray2D(tj, 1, inum, backend);
//     printArray2D(atomtype, 1, inum, backend);
//     error("here");
        
    START_TIMING;
    ZeroUarraytot2(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);    
    END_TIMING(54);
        
    
    START_TIMING;
//     AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, pairnum, pairnumsum, 
//             map, tj, idxu_max, nelem, inum, ineighmax, chemflag, backend);
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, aii, tj, idxu_max, 
            inum, ineighmax, chemflag, backend);    
    END_TIMING(55);
    
    dstype *Utot_r, *Utot_i;
    TemplateMalloc(&Utot_r, inum*idxu_max*nelem, backend);  
    TemplateMalloc(&Utot_i, inum*idxu_max*nelem, backend);  
    ArrayCopy(Utot_r, ulisttot_r, inum*idxu_max*nelem, backend);
    ArrayCopy(Utot_i, ulisttot_i, inum*idxu_max*nelem, backend);    
            
//     dstype *Stot_r, *Stot_i;
//     TemplateMalloc(&Stot_r, inum*idxu_max*nelem, backend);  
//     TemplateMalloc(&Stot_i, inum*idxu_max*nelem, backend);  
//     ArraySetValue(Stot_r, 0.0, inum*idxu_max*nelem, backend);
//     ArraySetValue(Stot_i, 0.0, inum*idxu_max*nelem, backend);    
//     ComputeUi(Stot_r, Stot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, map, aii, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);        
//     AddWself2Ui(Stot_r, Stot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);        
        
//     printui(ulisttot_r, inum, 10, twojmax);
//     printui(Stot_r, inum, 10, twojmax);
//     
//     printui(ulisttot_i, inum, 10, twojmax);
//     printui(Stot_i, inum, 10, twojmax);
    
    //ArrayAXPBY(Stot_r, Stot_r, ulisttot_r, one, minusone, inum*idxu_max*nelem, backend);    
    //dstype norme = PNORM(common.cublasHandle, inum*idxu_max*nelem, Stot_r, backend);
    //cout<<"Error : "<<norme<<endl;        
    //printArray2D(idxu_block, 1, twojmax+1, backend);
    //error("here");
    

//     writearray2file("ulisttotr2.bin", ulisttot_r, inum*idxu_max*nelem, backend);  
//     writearray2file("dulistr2.bin", dulist_r, 3*idxu_max*ineighmax, backend);        
    
     //printArray2D(x, 3, inum, backend);
     //printArray2D(ulisttot_r, inum, 12, backend);
     //printArray2D(pairnum, 1, inum, backend);
//     printf("%i %i %i %i",nelements, ndoubles, ntriples, nelem);
//     
    
    START_TIMING;
    ComputeZi2(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
    END_TIMING(57);
    
//    writearray2file("zlistr2.bin", zlist_r, idxz_max*inum*ndoubles, backend);    
    
//    printArray2D(ulisttot_r, inum, 12, backend);
//     display(" ");
    //printArray2D(zlist_r, inum, 12, backend);    
    //error("here");
    
    START_TIMING;
    ComputeBi2(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);           
    END_TIMING(58);
        
//    writearray2file("blist2.bin", blist, idxb_max*inum*ntriples, backend);       
//     dstype *dblist;  TemplateMalloc(&dblist, ineighmax*3*(ncoeffall-1), backend);      
//     ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, aii, tj, 
//             twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, inum, ineighmax, backend);
    
    // template <typename T> void cpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, T *cglist, T *bzero, 
//         int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
//         int idxz_max, int nelements, int bnorm_flag, int bzero_flag, int wselfall_flag, int inum)

//     dstype *blist2;
//     TemplateMalloc(&blist2, inum*idxb_max*nelem*nelem*nelem, backend);  
//     ComputeBi(blist2, zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, bzero, idxb, idxcg_block, idxu_block, idxz_block, 
//             twojmax, idxb_max, idxu_max, idxz_max, nelem, bnormflag, bzeroflag, wselfallflag, inum, backend);           
//     ArrayAXPBY(blist2, blist2, blist, one, minusone, inum*idxb_max*nelem*nelem*nelem, backend);    
//     dstype norme = PNORM(common.cublasHandle, inum*idxb_max*nelem*nelem*nelem, blist2, backend);
//     cout<<"Error : "<<norme<<endl;              
    
//     template <typename T> void cpuSnapComputeEi(T *eatom, T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, T *cglist, 
//         T *bzero, T *coeffelem, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, 
//         int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, 
//         int bzero_flag, int wselfall_flag, int quadraticflag, int inum);
// 

    dstype *eatom2, *fatom2, *vatom2;
    TemplateMalloc(&eatom2, inum, backend);   
    TemplateMalloc(&fatom2, inum*3, backend);   
    TemplateMalloc(&vatom2, inum*6, backend);   
    ArrayCopy(eatom2, eatom, inum, backend);
    ArrayCopy(fatom2, fatom, inum*3, backend);
    ArrayCopy(vatom2, vatom, inum*6, backend);

    START_TIMING;
    SnapTallyEnergyFull2(eatom, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);    
    END_TIMING(59);
                
//     dstype *ea;
//     TemplateMalloc(&ea, inum, backend);  
//     ArraySetValue(ea, 0.0, inum, backend);
//     SnapTallyEnergyFull2(ea, blist, coeffelem, alist, map, atomtype, 
//         inum, ncoeff, ncoeffall, quadraticflag, backend);    
    
//     writearray2file("eatom2.bin", eatom, inum, backend);    
//     writearray2file("coeffelem2.bin", coeffelem, ntypes*ncoeffall, backend);    
    
    START_TIMING;    
    ComputeBeta(beta, blist, coeffelem, alist, map, atomtype, 
        inum, ncoeff, ncoeffall, quadraticflag, backend);
    END_TIMING(60);
    
//    writearray2file("beta2.bin", beta, ncoeff*inum, backend);    
    
    START_TIMING;    
//     ComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, beta, idxz, idxb_block, 
//           idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, bnormflag, inum, backend);        
    ComputeYi(ylist_r, ylist_i, zlist_r, zlist_i, beta, idxz, idxb_block, 
          twojmax, idxb_max, idxu_max, idxz_max, nelem, bnormflag, inum, backend);            
    END_TIMING(61);
        
    //printf("%i %i\n",bzeroflag,  wselfallflag);
    
//     dstype *Y_r, *Y_i;
//     TemplateMalloc(&Y_r, inum*idxu_max*nelem, backend);  
//     TemplateMalloc(&Y_i, inum*idxu_max*nelem, backend);  
//     ComputeYi(Y_r, Y_i, Stot_r, Stot_i, cglist, beta, idxz, idxb_block, 
//           idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, bnormflag, inum, backend);        
//     ArrayAXPBY(Y_r, Y_r, ylist_r, one, minusone, inum*idxu_max*nelem, backend);    
//     dstype norma = PNORM(common.cublasHandle, inum*idxu_max*nelem, Y_r, backend);
//     ArrayAXPBY(Y_i, Y_i, ylist_i, one, minusone, inum*idxu_max*nelem, backend);    
//     dstype normb = PNORM(common.cublasHandle, inum*idxu_max*nelem, Y_i, backend);
//     cout<<"Error : "<<norma<<endl;  
//     cout<<"Error : "<<normb<<endl;  
    
//     dstype *ei;
//     TemplateMalloc(&ei, inum, backend);  
//     ArraySetValue(ei, 0.0, inum, backend);
//     ComputeEi(ei, ylist_r, ylist_i, Stot_r, Stot_i, idxu_max, nelem, inum, backend);
// //     SnapTallyEnergyFull2(ei, blist, coeffelem, alist, map, atomtype, 
// //         inum, ncoeff, ncoeffall, quadraticflag, backend);    
//     ArrayAXPBY(ei, ei, ea, one, minusone, inum, backend);    
//     dstype norma = PNORM(common.cublasHandle, inum, ei, backend);
//     cout<<"Error : "<<norma<<endl;        
//     error("here");
    
    
//    writearray2file("ylistr2.bin", ylist_r, idxu_max*inum*nelem, backend);    
    
    START_TIMING;
    ComputeDeidrj(dedr, ylist_r, ylist_i, dulist_r, dulist_i, idxu_block, map, 
          ai, tj, twojmax, idxu_max, chemflag, inum, ineighmax, backend);            
    END_TIMING(62);
        
    dstype *Stot_r, *Stot_i;
    TemplateMalloc(&Stot_r, inum*idxu_max*nelem, backend);  
    TemplateMalloc(&Stot_i, inum*idxu_max*nelem, backend);  
    ArraySetValue(Stot_r, 0.0, inum*idxu_max*nelem, backend);
    ArraySetValue(Stot_i, 0.0, inum*idxu_max*nelem, backend);            
    SnapComputeUi(Stot_r, Stot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, map, aii, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);        
    AddWself2Ui(Stot_r, Stot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);        
                    
    dstype *Y_r, *Y_i;
    TemplateMalloc(&Y_r, inum*idxu_max*nelem, backend);  
    TemplateMalloc(&Y_i, inum*idxu_max*nelem, backend);  
//     SnapComputeYi(Y_r, Y_i, Stot_r, Stot_i, cglist, beta, idxz, idxb_block, 
//           idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, bnormflag, inum, backend);        
    SnapComputeYi(Y_r, Y_i, Stot_r, Stot_i, cglist, coeffelem, map, atomtype, idxz, idxb_block, 
          idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, ncoeffall, bnormflag, inum, backend);        
    
//     printf("%i %i %i\n", idxu_max, idxb_max, idxz_max);
//     printArray2D(idxb, 3, idxb_max, backend);
//     printArray2D(idxz, 10, idxz_max, backend);
//     int jdim = twojmax+1;
//     for(int jjb = 0; jjb < idxb_max; jjb++) {
//         const int j1 = idxb[jjb*3 + 0];
//         const int j2 = idxb[jjb*3 + 1];
//         const int j = idxb[jjb*3 + 2];  
//         int n;
//         if (j % 2 == 0) 
//             n = (j+1)*(j/2) + 1 + j/2;
//         else
//             n = (j+1)*((j+1)/2);
//         
//         printf("%i  %i  %i   %i\n",j1, j2, j, n);  
//     }
//     for(int jjb = 0; jjb < idxz_max; jjb++) {
//         const int j1 = idxz[jjb*10 + 0];
//         const int j2 = idxz[jjb*10 + 1];
//         const int j = idxz[jjb*10 + 2];  
//         int n;
//         if (j % 2 == 0) 
//             n = (j+1)*(j/2) + 1 + j/2;
//         else
//             n = (j+1)*((j+1)/2);
//         
//         printf("%i  %i  %i   %i   %i\n",j1, j2, j, n, idxz_block[j + j2*jdim + j1*jdim*jdim]);  
//     }        
//     for(int j1 = 0; j1 <= twojmax; j1++)
//         for(int j2 = 0; j2 <= j1; j2++)
//           for(int j = j1 - j2; j <= min(twojmax, j1 + j2); j += 2) 
//             printf("%i  %i  %i   %i\n",j1, j2, j, idxz_block[j + j2*jdim + j1*jdim*jdim]);            
//     TemplateMalloc(&dU_r, ineighmax*idxu_max*3, backend);  
//     TemplateMalloc(&dU_i, ineighmax*idxu_max*3, backend);  
//     ArraySetValue(dU_r, 0.0, ineighmax*idxu_max*3, backend);
//     ArraySetValue(dU_i, 0.0, ineighmax*idxu_max*3, backend);       
    
    SnapComputeEi(eatom2, Stot_r, Stot_i, cglist, bzero, coeffelem, map, atomtype, 
            idxb, idxcg_block, idxu_block, twojmax, idxb_max, idxu_max, nelem, 
            ncoeffall, bnormflag, bzeroflag, wselfallflag, quadraticflag, inum, backend);           
    ArrayAXPBY(eatom2, eatom2, eatom, one, minusone, inum, backend);    
    dstype norme = PNORM(common.cublasHandle, inum, eatom2, backend);
    cout<<"Error : "<<norme<<endl;              
    
    SnapComputeFi(fatom2, vatom2, Y_r, Y_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
     idxu_block, map, ai, aj, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);  
        
//         ComputeDeidrj(fatom2, Y_r, Y_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//                 idxu_block, map, ai, aj, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);                    
//     ComputeEi(eatom2, ylist_r, ylist_i, Stot_r, Stot_i, twojmax, idxu_max, nelem, inum, backend);
//     printArray2D(eatom, 1, 12, backend);
//     printArray2D(eatom2, 1, 12, backend);
//     for (int n=0; n<12; n++)
//         cout<<eatom2[n]/eatom[n]<<"  ";
//     cout<<endl;
//     ArrayAXPBY(eatom2, eatom2, eatom, one, minusone, inum, backend);
//     printArray2D(eatom2,  1, 12, backend);
//     printArray2D(coeffelem, 1, 12, backend);    
    
//     AddWself2Ui(Stot_r, Stot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);            
//     ArrayAXPBY(Stot_r, Stot_r, Utot_r, one, minusone, inum*idxu_max*nelem, backend);    
//     dstype norme = PNORM(common.cublasHandle, inum*idxu_max*nelem, Stot_r, backend);
//     cout<<"Error : "<<norme<<endl;      
//     
//     ArraySetValue(Stot_r, 0.0, inum*idxu_max*nelem, backend);
//     ArraySetValue(Stot_i, 0.0, inum*idxu_max*nelem, backend);        
//     ComputeUi(Stot_r, Stot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, map, aii, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);        
//     AddWself2Ui(Stot_r, Stot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);        
//     ArrayAXPBY(Stot_r, Stot_r, Utot_r, one, minusone, inum*idxu_max*nelem, backend);    
//     dstype normb = PNORM(common.cublasHandle, inum*idxu_max*nelem, Stot_r, backend);    
//     cout<<"Error : "<<normb<<endl;  
//         
//     setui(dulist_r, ineighmax, twojmax, idxu_max);
//     setui(dulist_i, ineighmax, twojmax, idxu_max);
// //     setui(dU_r, ineighmax, twojmax, idxu_max);
// //     setui(dU_i, ineighmax, twojmax, idxu_max);
//     
// //     printui(dulist_r, ineighmax, 0, twojmax);
// //     printui(dU_r, ineighmax, 0, twojmax);
//         
//     ComputeDeidrj(dedr2, ylist_r, ylist_i, dU_r, dU_i, idxu_block, map, 
//           ai, tj, twojmax, idxu_max, chemflag, inum, ineighmax, backend);                
//     ArrayAXPBY(dedr2, dedr2, dedr, one, minusone, ineighmax*3, backend);    
//     norma = PNORM(common.cublasHandle, ineighmax*3, dedr2, backend);
//     cout<<"Error : "<<norma<<endl;         
//     
//     ArrayAXPBY(dU_r, dU_r, dulist_r, one, minusone, ineighmax*idxu_max*3, backend);    
//     normb = PNORM(common.cublasHandle, ineighmax*idxu_max*3, dU_r, backend);    
//     cout<<"Error : "<<normb<<endl;      
//     error("here");
    
//     ComputeDeidrj(dstype *dedr, dstype *ylist_r, dstype *ylist_i, dstype *rootpqarray, dstype *rij, 
//         dstype *wjelem, dstype *radelem,  dstype rmin0, dstype rfac0, dstype rcutfac, int *idxu_block, 
//         int *map, int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int inum, int ijnum, 
//         int switchflag, int chemflag, int backend)        
    
//    writearray2file("dedr2.bin", dedr, 3*ineighmax, backend);   
    
    START_TIMING;
    SnapTallyForceFull(fatom, dedr, ai, aj, alist, ineighmax, backend);
    END_TIMING(63);
        
    START_TIMING;
    SnapTallyVirialFull(vatom, dedr, rij, ai, aj, inum, ineighmax, backend);  
    END_TIMING(64);    
    
    ArrayAXPBY(fatom2, fatom2, fatom, one, minusone, inum*3, backend);    
    dstype norma = PNORM(common.cublasHandle, inum*3, fatom2, backend);
    cout<<"Error : "<<norma<<endl;        
    
    ArrayAXPBY(vatom2, vatom2, vatom, one, minusone, inum*6, backend);    
    norma = PNORM(common.cublasHandle, inum*6, vatom2, backend);
    cout<<"Error : "<<norma<<endl;        
    //error("here");
    
//     printf("%i %i %i %i %i %i %i %i %i %i %i\n", common.anum, inum, ineighmax, idxu_max, idxb_max, 
//             idxz_max, ncoeff, ncoeffall, chemflag, bnormflag,  nelem);
//     error("here");
    
//     dstype *bi, *bd, *bv, *e, *f, *ri1, *bi2, *blist2;
//     int M = ncoeffall-1;
//     int nperdim = M;
//     TemplateMalloc(&bi, M, backend);  
//     TemplateMalloc(&bd, inum*dim*M, backend);  
//     TemplateMalloc(&bv, 6*M, backend);  
//     TemplateMalloc(&e, inum, backend);  
//     TemplateMalloc(&f, inum*3, backend);  
//     TemplateMalloc(&ri1, inum*3, backend);  
//     TemplateMalloc(&bi2, M, backend);      
//     TemplateMalloc(&blist2, inum*M, backend);  
//     
//     SnapTallyBispectrum(bi, blist, alist, atomtype, inum, ncoeff, nperdim, ntypes, quadraticflag, backend);    
//     SnapTallyBispectrumDeriv(bd, blist, dblist, aii, ai, aj, 
//           ti, inum, ineighmax, ncoeff, nperdim, ntypes, quadraticflag, backend);        
//     SnapTallyBispectrumVirial(bv, blist, dblist, rij, aii, ai, aj, 
//           ti, inum, ineighmax, ncoeff, nperdim, ntypes, quadraticflag, backend);                
//     ArraySetValue(e, 0.0, inum, backend);
//     SnapTallyEnergyFull2(e, blist, coeffelem, alist, map, atomtype, inum, ncoeff, ncoeffall, quadraticflag, backend);    
//     dstype energy = ArraySum(e, tmp.tmpmem, inum, backend);    
//     dstype energy2;
//     PDOT(common.cublasHandle, M*ntypes, bi, inc1, &coeffelem[1], inc1, &energy2, backend);
//     energy2 += inum*coeffelem[0];    
//     printf("Energy: %g  %g\n", energy, energy2);        
//     
//     ArraySetValue(f, 0.0, inum*3, backend);
//     SnapTallyForceFull(f, dedr, ai, aj, alist, ineighmax, backend);
//     for (int i=0; i<inum*3; i++) {
//         dedr[i] = 0.0;
//         for (int j=0; j<M*ntypes; j++)
//             dedr[i] += bd[i + inum*3*j]*coeffelem[1+j];
//     }    
//     printArray2D(f, 3, inum, backend);
//     printArray2D(dedr, 3, inum, backend);    
//         
//     dstype epsil = 1e-6;
//     for (int i = 0; i<dim*ineighmax; i++) {
//         cpuArrayCopy(ri1, rij, dim*ineighmax);  
//         ri1[i] = ri1[i] + epsil;
//         int d = i%3;
//         int m = (i-d)/3;
//                 
//         ComputeSij(ulist_r, ulist_i, &dulist_r[0*n], &dulist_i[0*n], &dulist_r[1*n], &dulist_i[1*n],
//                     &dulist_r[2*n], &dulist_i[2*n], rootpqarray, ri1, wjelem, radelem, rmin0, rfac0, rcutfac, 
//                     idxu_block, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);      
//     
//         ZeroUarraytot2(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
//                 map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);    
//         
//         AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, ai, tj, idxu_max, 
//                 inum, ineighmax, chemflag, backend);    
// 
//         ComputeZi2(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
//               idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
// 
//         ComputeBi2(blist2, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
//               map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
//               nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);                   
// 
//         SnapTallyBispectrum(bi2, blist2, alist, atomtype, inum, ncoeff, nperdim, ntypes, quadraticflag, backend);    
// 
//         printArray2D(bi, 1, M, common.backend);
//         printArray2D(bi2, 1, M, common.backend);
//         for (int j=0; j<M; j++) 
//             cout<<(bi2[j]-bi[j])/epsil<<" ";
//         cout<<endl;
//         for (int j=0; j<M; j++) 
//             cout<<dblist[m + ineighmax*d + dim*ineighmax*j]<<" ";
//         cout<<endl;
//         
// //         for (int j=0; j<M; j++) 
// //             bi2[j] = fabs((bi2[j]-bi[j])/epsil - dblist[i + dim*ineighmax*j]);
// //         cout<<"Maximum absolute error: "<<cpuArrayMax(bi2, M)<<endl;   
//         error("here");
//     }
//     
//     
//         
    //PGEMNV(common.cublasHandle, inum*3, M*ntypes, &one, bd, M*ntypes, &coeffelem[1], inc1, &zero, dedr, inc1, backend);    
    //PGEMNV(CCal.common.cublasHandle, M, M, &one, A, M, b, inc1, &zero, c, inc1, backend);    
//     ArrayAXPBY(f, f, dedr, one, minusone, inum*3, backend);    
//     dstype norme = PNORM(common.cublasHandle, inum*3, f, backend);
//     
//     printf("Force error: %g\n", norme);    
    
//     TemplateFree(dblist, backend);  
//     TemplateFree(bi, backend);  
//     TemplateFree(bd, backend);  
//     TemplateFree(bv, backend);     
//     TemplateFree(e, backend);     
}


void ComputePairSnap3(dstype *eatom, dstype *fatom, dstype *vatom, dstype *x, snastruct &sna, commonstruct &common, 
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
        
    INIT_TIMING;
    
    START_TIMING;
    NeighborPairList(pairnum, pairlist, x, rcutsq, alist, neighlist, 
                            neighnum, atomtype, alist, inum, neighmax, dim, ntypes, backend);
    END_TIMING(50);
            
    START_TIMING;
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    END_TIMING(51);
    
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
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];    
    ne += idxu_max*nelements*na;        
    dstype *ylist_r = &tmp.tmpmem[ne]; 
    ne += idxu_max*nelements*na;
    dstype *ylist_i = &tmp.tmpmem[ne];     
              
    START_TIMING;
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
    END_TIMING(52);
                    
    ArraySetValue(ulisttot_r, 0.0, inum*idxu_max*nelem, backend);
    ArraySetValue(ulisttot_i, 0.0, inum*idxu_max*nelem, backend);            
    SnapComputeUi(ulisttot_r, ulisttot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, map, aii, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);        
    AddWself2Ui(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);        

    SnapComputeEi(eatom, ulisttot_r, ulisttot_i, cglist, bzero, coeffelem, map, atomtype, 
            idxb, idxcg_block, idxu_block, twojmax, idxb_max, idxu_max, nelem, 
            ncoeffall, bnormflag, bzeroflag, wselfallflag, quadraticflag, inum, backend);           
    
    SnapComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, coeffelem, map, atomtype, 
          idxz, idxb_block, idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, ncoeffall, bnormflag, inum, backend);        
            
    SnapComputeFi(fatom, vatom, ylist_r, ylist_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        idxu_block, map, ai, aj, ti, tj, twojmax, idxu_max, inum, ineighmax, switchflag, chemflag, backend);                      
}

void ComputeSnap(dstype *bi, dstype *bd, dstype *bv, dstype *x, snastruct &sna, commonstruct &common, 
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
    int nperdim = sna.nperdim;
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
        
    INIT_TIMING;
    
    START_TIMING;
    //dstype *x = &sys.x[0];
//     NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
//         neighnum, inum, neighmax, dim, backend);
    NeighborPairList(pairnum, pairlist, x, rcutsq, alist, neighlist, 
                            neighnum, atomtype, alist, inum, neighmax, dim, ntypes, backend);    
    END_TIMING(50);
    
    START_TIMING;
    Cumsum(pairnumsum, pairnum, &tmp.intmem[ne+na+1], &tmp.intmem[2*na+ne+2], na+1, backend);                                         
    int ineighmax = IntArrayGetValueAtIndex(pairnumsum, na, backend);        
    END_TIMING(51);
    
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
    dstype *blist = &tmp.tmpmem[ne]; // idxb_max*ntriples*na          
    ne += idxb_max*ntriples*na;    
    dstype *dblist = &tmp.tmpmem[ne]; // idxb_max*ntriples*dim*ineighmax          
    dstype *ulisttot_r = &tmp.tmpmem[ne];
    ne += idxu_max*nelements*na;
    dstype *ulisttot_i = &tmp.tmpmem[ne];
    //ne += idxu_max*nelements*na;            
                  
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
        
    int n = idxu_max*ineighmax;
    ComputeSij(ulist_r, ulist_i, &dulist_r[0*n], &dulist_i[0*n], &dulist_r[1*n], &dulist_i[1*n],
            &dulist_r[2*n], &dulist_i[2*n], rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
            idxu_block, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);      
    
    ZeroUarraytot2(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, inum, backend);    
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, ai, tj, idxu_max, 
            inum, ineighmax, chemflag, backend);    
        
    ComputeZi2(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, inum, backend);
        
    ComputeBi2(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bzeroflag,  wselfallflag, chemflag, inum, backend);           
        
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, ai, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, inum, ineighmax, backend);
    
    SnapTallyBispectrum(bi, blist, alist, atomtype, inum, ncoeff, nperdim, ntypes, quadraticflag, backend);    
    
    SnapTallyBispectrumDeriv(bd, blist, dblist, aii, ai, aj, 
          ti, inum, ineighmax, ncoeff, nperdim, ntypes, quadraticflag, backend);    
    
    SnapTallyBispectrumVirial(bv, blist, dblist, rij, aii, ai, aj, 
          ti, inum, ineighmax, ncoeff, nperdim, ntypes, quadraticflag, backend);        
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
        
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
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
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
    
//     ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
    
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
       
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
        
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
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
        
//     NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
//           pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
    
//     ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
//     ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
//       rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
//       ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
      ai, aj, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    if (quadraticflag) {
        ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
              map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
    }
          
//     ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
//         idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
//         idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);

    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, map, ai, aj, ti, tj, twojmax, 
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
        
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
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
//         
//     NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
//           pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);
       
    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
       
//     ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
//     ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
//       rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
//       ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
      ai, aj, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    if (quadraticflag) {
        ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
              map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
                nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
    }
          
//     ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
//         idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
//         idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, map, ai, aj, ti, tj, twojmax, 
        idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeSnav(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);
        
//     ComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
//           ti, mask, ncoeff, ntypes, nperdim, ineighmax, quadraticflag, backend);        
    
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
        
    dstype *x = &sys.x[0];
    NeighborPairList(pairnum, pairlist, x, rcutmax*rcutmax, alist, neighlist, 
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
        
//     NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
//           pairlist, pairnumsum, alist, atomtype, inum, neighmax, dim, backend);

    NeighborPairs(rij, x, aii, ai, aj, ti, tj, pairnum, 
          pairlist, pairnumsum, alist, atomtype, alist, inum, neighmax, dim, backend);
       
//     ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, atomtype, ai, aj, twojmax, idxu_max, ineighmax, backend);
    ComputeUij(ulist_r, ulist_i, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, 
          idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ineighmax, backend);        
       
    ZeroUarraytot(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
            map, ai, wselfallflag, chemflag, idxu_max, nelements, twojmax, inum, backend);
    
    AddUarraytot(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, rcutfac, 
     idxu_block, alist, atomtype, pairnum, pairnumsum, map, tj, twojmax, idxu_max, nelements, 
          inum, switchflag, chemflag, backend);
    
//     ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
//       rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, atomtype, 
//       ai, aj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeDuijdrj(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, 
      rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
      ai, aj, ti, tj, twojmax, idxu_max, ineighmax, switchflag, backend);     
    
    ComputeZi(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelements, bnormflag, inum, backend);
        
    ComputeBi(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
            nelements, bzeroflag,  wselfallflag, chemflag, inum, backend);
          
//     ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
//         idxb, idxu_block, idxz_block, atomtype, map, ai, aj, twojmax, 
//         idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, ineighmax, backend);
    
    ComputeDbidrj(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
        idxb, idxu_block, idxz_block, map, ai, aj, ti, tj, twojmax, 
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
 

