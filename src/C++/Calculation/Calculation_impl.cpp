#ifndef __CALCULATIONIMPL
#define __CALCULATIONIMPL

void cpuComputeSingleEnergy(dstype *ei, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int inum, int potnum)
{
}
void cpuComputeSingleForce(dstype *fi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int inum, int potnum)
{
}
void cpuComputeSingleParam(dstype *gi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int inum, int potnum)
{
}
void cpuComputeSingleEnergyForce(dstype *ei, dstype *fi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int inum, int potnum)
{
}
void cpuComputeSingle(dstype *ei, dstype *fi, dstype *gi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int inum, int potnum)
{
}

void cpuComputePairForce(dstype *fij, dstype *xij, dstype *qi, dstype *qj, int *ti, int *tj, 
        int *ai, int *aj, dstype *pairparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}
void cpuComputePairEnergyForce(dstype *eij, dstype *fij, dstype *xij, dstype *qi, dstype *qj, int *ti, int *tj, 
        int *ai, int *aj, dstype *pairparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}

void cpuElectronDensity(dstype *ei, dstype *eij, int *pairnum, int *pairnumsum, int inum) 
{
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        ei[i] = 0.0;
        for (int j=0; j<jnum; j++) 
            ei[i] += eij[start+j];        
    }
}
void cpuComputeEmbedingEnergy(dstype *ei, dstype *gi, dstype *rhoi, int inum) 
{
    
}
void cpuComputeEmbedingForce(dstype *fij, dstype *gi, int *pairnum, int *pairnumsum, int inum)
{    
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        for (int j=0; j<jnum; j++)             
            fij[start+j] = gi[i]*fij[start+j];                
    }        
}
void cpuComputeTripletForce(dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
     int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *tripletparam, int dim, int ncq, int nparam, int ijknum, int potnum)
{
}
void cpuComputeTripletEnergyForce(dstype *eijk, dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
     int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *tripletparam, int dim, int ncq, int nparam, int ijknum, int potnum)
{
}
void cpuComputeQuadrupletForce(dstype *fij, dstype *fik, dstype *fil, dstype *xij, dstype *xik, dstype *xil, 
      dstype *qi, dstype *qj, dstype *qk, dstype *ql, int *ti, int *tj, int *tk, int *tl, 
      int *ai, int *aj, int *ak, int *al, dstype *tripletparam, int dim, int ncq, int nparam, int ijklnum, int potnum)
{
}
void cpuComputeQuadrupletEnergyForce(dstype *eijkl, dstype *fij, dstype *fik, dstype *fil, dstype *xij, dstype *xik, dstype *xil, 
      dstype *qi, dstype *qj, dstype *qk, dstype *ql, int *ti, int *tj, int *tk, int *tl, 
      int *ai, int *aj, int *ak, int *al, dstype *tripletparam, int dim, int ncq, int nparam, int ijklnum, int potnum)
{
}

void implSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *ai = &tmp.intmem[na]; // na
        Int *ti = &tmp.intmem[2*na]; // na
        dstype *xi = &tmp.tmpmem[0]; // na*dim
        dstype *qi = &tmp.tmpmem[na*dim]; // na*ncq
        cpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, na, ncq, dim);

        dstype *fi = &tmp.tmpmem[na*(dim+ncq)]; // na*dim
        dstype *ei = &tmp.tmpmem[na*(2*dim+ncq)]; // na*dim
        cpuComputeSingleEnergyForce(ei, fi, xi, qi, ti, ai, param, dim, ncq, nparam, na, potnum);
        cpuSingleDecomposition(e, f, ei, fi, ai, na, dim);        
    }        
}

void implNonbondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int potnum)
{    
    implSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, atomtype, nparam, 0, potnum);            
}

void implBondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{
    implSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, atomtype, nparam, typei, potnum);            
}

void implFullNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, ntuples, potnum);
    
        if (decomp==0)
            cpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);        
        
//         if (dim==2) {
//             if (decomp==0)
//                 cpuFullForceDecomposition2D(f, fij, ai, aj, ntuples);        
//             else
//                 cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);                        
//         }
//         else {
//             if (decomp==0)
//                 cpuFullForceDecomposition3D(f, fij, ai, aj, ntuples);    
//             else
//                 cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
//         }                
    }        
}

void implFullNeighNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int decomp, Int potnum)
{    
    implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, 0, 0, decomp, potnum);            
}

void implFullNeighBondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{
    if (typei==typej) {
        implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);        
    }
    else {
        implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);   
        implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typej, typei, decomp, potnum);        
    }
}


void implHalfNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);            
        else
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, dim);
                        
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, ntuples, potnum);

        if (decomp==0)
            cpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, ntuples, dim);
        else {
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);   
            cpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
            cpuNeighborAtomPairDecomposition(e, f, eij, fij, jlist, bnumsum, index, naj, dim);                 
        }
                
//         if (dim==2) {
//             if (decomp==0)
//                 cpuHalfForceDecomposition2D(f, fij, ai, aj, ntuples);       
//             else {
//                 cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
//                 cpuArrayCopy(tmp.intmem, aj, ntuples);
//                 Int *jlist = &tmp.intmem[ntuples];   
//                 Int *bnumsum = &tmp.intmem[2*ntuples]; 
//                 Int *index = &tmp.intmem[3*ntuples]; // ntuples       
//                 Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
//                 Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
//                 Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
//                 Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//             }
//         }
//         else {
//             if (decomp==0)
//                 cpuHalfForceDecomposition3D(f, fij, ai, aj, ntuples);        
//             else {
//                 cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);      
//                 cpuArrayCopy(tmp.intmem, aj, ntuples);
//                 Int *jlist = &tmp.intmem[ntuples];   
//                 Int *bnumsum = &tmp.intmem[2*ntuples]; 
//                 Int *index = &tmp.intmem[3*ntuples]; // ntuples       
//                 Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
//                 Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
//                 Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
//                 Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);                                
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                   
//             }
//         }            
    }        
}

void implHalfNeighNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int decomp, Int potnum)
{    
    implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, 0, 0, decomp, potnum);            
}

void implHalfNeighBondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{
    if (typei==typej) {
        implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);        
    }
    else {
        implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);   
        implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typej, typei, decomp, potnum);        
    }
}

void implBO2EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq+1)]; // ntuples*dim        
        cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, ntuples, potnum); 
        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        dstype *gi = &tmp.tmpmem[2*na]; // na
        cpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        cpuComputeEmbedingEnergy(ei, gi, rhoi, na); 
        cpuComputeEmbedingForce(fij, gi, pairnum, pairnumsum, na);        
        
        cpuPutArrayAtIndex(e, ei, ilist, na);
        
        if (dim==2) {
            if (decomp==0)
                cpuHalfForceDecomposition2D(f, fij, ai, aj, ntuples);        
            else {                
                cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
                cpuArrayCopy(tmp.intmem, aj, ntuples);
                Int *jlist = &tmp.intmem[ntuples];   
                Int *bnumsum = &tmp.intmem[2*ntuples]; 
                Int *index = &tmp.intmem[3*ntuples]; // ntuples       
                Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
                Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
                Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
                Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
                cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);                
            }
        }
        else {
            if (decomp==0)
                cpuHalfForceDecomposition3D(f, fij, ai, aj, ntuples);    
            else {
                cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
                cpuArrayCopy(tmp.intmem, aj, ntuples);
                Int *jlist = &tmp.intmem[ntuples];   
                Int *bnumsum = &tmp.intmem[2*ntuples]; 
                Int *index = &tmp.intmem[3*ntuples]; // ntuples       
                Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
                Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
                Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
                Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
                cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                                   
            }
        }                
    }        
}

void implBO3EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum       
        if (typej == typei) {
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);        
        }
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1    
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int npairs = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // npairs    
        Int *aj = &tmp.intmem[1+3*na+npairs+na*jnum]; // npairs        
        Int *ti = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs        
        Int *tj = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs        
        dstype *eij = &tmp.tmpmem[0]; // npairs
        dstype *fij = &tmp.tmpmem[npairs]; // npairs*dim
        dstype *xij = &tmp.tmpmem[npairs*(1+dim)]; // npairs*dim
        dstype *qi = &tmp.tmpmem[npairs*(1+2*dim)]; // npairs*ncq
        dstype *qj = &tmp.tmpmem[npairs*(1+2*dim+ncq)]; // npairs*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                                
        cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, npairs, potnum);
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs*jnum        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, typek, dim);                
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[1+3*na+3*npairs+(na+npairs)*jnum]; // npairs+1        
        cpuCumsum(tripletnumsum, tripletnum, npairs+1);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, npairs, common.backend);     
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj, tripletnum, tripletnumsum, tripletlist)         
        Int *a3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum]; // ntuples        
        Int *a3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+ntuples]; // ntuples   
        Int *a3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+2*ntuples]; // ntuples   
        Int *t3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+3*ntuples]; // ntuples        
        Int *t3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+4*ntuples]; // ntuples        
        Int *t3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+5*ntuples]; // ntuples        
        dstype *x3ij = &tmp.tmpmem[npairs*(1+dim)]; // ntuples*dim
        dstype *x3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*dim]; // ntuples*dim
        dstype *q3i = &tmp.tmpmem[npairs*(1+dim)+2*ntuples*dim]; // ntuples*ncq
        dstype *q3j = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+ncq)]; // ntuples*ncq
        dstype *q3k = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+2*ncq)]; // ntuples*ncq
        cpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, jnum, npairs, ncq, dim);                      
        
        dstype *f3ij = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *f3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *e3ijk = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq)]; // ntuples
        cpuComputeTripletEnergyForce(e3ijk, f3ij, f3ik, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                                param, dim, ncq, nparam, ntuples, potnum);
        
        dstype *h3ij = &tmp.tmpmem[npairs*(1+dim)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2+dim)]; // npairs
        dstype *d3ij = &tmp.tmpmem[npairs*(3+dim)]; // npairs
        dstype *g3ij = &tmp.tmpmem[npairs*(4+dim)]; // dim*npairs
        cpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        cpuElectronDensity(g3ij, f3ij, tripletnum, tripletnumsum, npairs);
        cpuComputeEmbedingEnergy(c3ij, d3ij, h3ij, npairs); 
        
        for (int i=0; i<npairs; i++) 
            for (int j=0; j<dim; j++)
                fij[j+3*i] = fij[j+3*i]*c3ij[i] + eij[i]*d3ij[i]*g3ij[j+3*i];   
        
        for (int i=0; i<npairs; i++) 
            eij[i] = eij[i]*c3ij[i];
                    
        // pairnum, pairnumsum, pairlist, tripletnum, tripletlist, tripletnumsum, t3i, t3j, t3k -> ai, aj, a3i, a3j, a3k
        for (int i=0; i<npairs; i++) {
            tmp.intmem[i] = tmp.intmem[1+3*na+na*jnum+i]; // ai
            tmp.intmem[npairs+i] = tmp.intmem[1+3*na+npairs+na*jnum+i]; // aj
        }
        for (int i=0; i<ntuples; i++) {
            tmp.intmem[2*npairs+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum]; //a3i
            tmp.intmem[2*npairs+ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+ntuples]; //a3j
            tmp.intmem[2*npairs+2*ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+2*ntuples]; //a3k
        }
                
        if (decomp==0)
            cpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, npairs, dim);
        else {
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);   
            Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //npairs  
            Int *bnumsum = &tmp.intmem[3*npairs+3*ntuples]; //npairs  
            Int *index = &tmp.intmem[4*npairs+3*ntuples]; // npairs
            Int *p0 = &tmp.intmem[5*npairs+3*ntuples]; // npairs       
            Int *p1 = &tmp.intmem[6*npairs+3*ntuples]; // npairs       
            Int *p2 = &tmp.intmem[7*npairs+3*ntuples]; // npairs       
            Int *p3 = &tmp.intmem[8*npairs+3*ntuples]; // npairs                   
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
            cpuNeighborAtomPairDecomposition(e, f, eij, fij, jlist, bnumsum, index, naj, dim);                 
        }
        
//         if (dim==2) {
//             if (decomp==0)
//                 cpuHalfForceDecomposition2D(f, fij, ai, aj, npairs);        
//             else {                
//                 cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //npairs  
//                 Int *bnumsum = &tmp.intmem[3*npairs+3*ntuples]; //npairs  
//                 Int *index = &tmp.intmem[4*npairs+3*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[5*npairs+3*ntuples]; // npairs       
//                 Int *p1 = &tmp.intmem[6*npairs+3*ntuples]; // npairs       
//                 Int *p2 = &tmp.intmem[7*npairs+3*ntuples]; // npairs       
//                 Int *p3 = &tmp.intmem[8*npairs+3*ntuples]; // npairs       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);                
//             }
//         }
//         else {
//             if (decomp==0)
//                 cpuHalfForceDecomposition3D(f, fij, ai, aj, npairs);    
//             else {
//                 cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //npairs  
//                 Int *bnumsum = &tmp.intmem[3*npairs+3*ntuples]; //npairs  
//                 Int *index = &tmp.intmem[4*npairs+3*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[5*npairs+3*ntuples]; // npairs       
//                 Int *p1 = &tmp.intmem[6*npairs+3*ntuples]; // npairs       
//                 Int *p2 = &tmp.intmem[7*npairs+3*ntuples]; // npairs       
//                 Int *p3 = &tmp.intmem[8*npairs+3*ntuples]; // npairs       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                                   
//             }
//         }                
        
        for (int i=0; i<npairs; i++) {
            int m = tripletnum[i];
            int s = tripletnumsum[i];
            for (int j=0; j<m; j++) 
                for (int d=0; d<dim; d++)
                    f3ik[d+3*(s+j)] = eij[i]*d3ij[i]*f3ik[d+3*(s+j)];                                                
        }                
        
        if (dim==2) {            
            if (decomp==0)
                cpuHalfForceDecomposition2D(f, f3ik, a3i, a3k, ntuples);        
            else {
                Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
                Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
                Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
                Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
                Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
                Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
                Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
                Int nai = cpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
                cpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nai);                            
                Int nak = cpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
                cpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nak);          
            }
        }
        else {
            if (decomp==0)
                cpuHalfForceDecomposition3D(f, f3ik, a3i, a3k, ntuples);    
            else {
                Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
                Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
                Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
                Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
                Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
                Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
                Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
                Int nai = cpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
                cpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nai);                            
                Int nak = cpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
                cpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nak);         
            }
        }                        
    }        
}

void implBondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum        
        if ((typej>0) && (typek>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        for (int ii=0; ii<na; ii++)
            tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        cpuCumsum(tripletnumsum, tripletnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
                
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        cpuNeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        cpuArrayCopy(tripletlist, temp, ntuples);
        
        Int *ai = &tmp.intmem[1 + 3*na + 2*ntuples]; // ntuples        
        Int *aj = &tmp.intmem[1 + 3*na + 3*ntuples]; // ntuples   
        Int *ak = &tmp.intmem[1 + 3*na + 4*ntuples]; // ntuples   
        Int *ti = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples        
        Int *tj = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples        
        Int *tk = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *xik = &tmp.tmpmem[ntuples*dim]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[2*ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(2*dim+ncq)]; // ntuples*ncq
        dstype *qk = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples*ncq
        cpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, ncq, dim);                     
              
        dstype *fij = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *eijk = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples*dim
        cpuComputeTripletEnergyForce(eijk, fij, fik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                                param, dim, ncq, nparam, ntuples, potnum);
     
        if (decomp==0)            
            cpuTripletDecomposition(e, f, eijk, fij, fik, ai, aj, ak, ntuples, dim);
        else {
            cpuCenterAtomTripletDecomposition(e, f, eijk, fij, fik, ilist, tripletnumsum, na, dim);                
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            cpuNeighborAtomTripletDecomposition(e, f, eijk, fij, jlist, bnumsum, index, naj, dim);
            Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            cpuNeighborAtomTripletDecomposition(e, f, eijk, fik, jlist, bnumsum, index, nak, dim);
        }
        
//         if (dim==2) {
//             if (decomp==0)
//                 cpuForceDecompositionTriplet2D(f, fij, fik, ai, aj, ak, ntuples);
//             else {
//                 cpuAtomDecompositionTriplet2D(f, fij, fik, ilist, tripletnumsum, na);                
//                 Int *jlist = &tmp.intmem[0];  // ntuples       
//                 Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
//                 Int *index = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fik, jlist, bnumsum, index, nak);
//             }
//         }
//         else {
//             if (decomp==0)                  
//                 cpuForceDecompositionTriplet3D(f, fij, fik, ai, aj, ak, ntuples);
//             else {
//                 cpuAtomDecompositionTriplet3D(f, fij, fik, ilist, tripletnumsum, na);       
//                 Int *jlist = &tmp.intmem[0];  // ntuples       
//                 Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
//                 Int *index = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fik, jlist, bnumsum, index, nak);                
//             }
//         }            
    }        
}

void implBondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum    
        if ((typej>0) && (typek>0) && (typel>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, typel, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        for (int ii=0; ii<na; ii++)
            quadrupletnum[ii] = (pairnum[ii]-2)*(pairnum[ii]-1)*pairnum[ii]/6;       
                        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        cpuCumsum(quadrupletnumsum, quadrupletnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        cpuNeighTripletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        cpuArrayCopy(quadrupletlist, temp, ntuples);
        
        Int *ai = &tmp.intmem[1 + 3*na + 3*ntuples]; // ntuples        
        Int *aj = &tmp.intmem[1 + 3*na + 4*ntuples]; // ntuples   
        Int *ak = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples   
        Int *al = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples   
        Int *ti = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples        
        Int *tj = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples        
        Int *tk = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples        
        Int *tl = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples                   
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *xik = &tmp.tmpmem[ntuples*dim]; // ntuples*dim
        dstype *xil = &tmp.tmpmem[2*ntuples*dim]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[3*ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(3*dim+ncq)]; // ntuples*ncq
        dstype *qk = &tmp.tmpmem[ntuples*(3*dim+2*ncq)]; // ntuples*ncq
        dstype *ql = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*ncq        
        cpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, ncq, dim);       
        
        dstype *fij = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(5*dim+3*ncq)]; // ntuples*dim
        dstype *fil = &tmp.tmpmem[ntuples*(6*dim+3*ncq)]; // ntuples*dim
        dstype *eijkl = &tmp.tmpmem[ntuples*(7*dim+3*ncq)]; // ntuples
        cpuComputeQuadrupletEnergyForce(eijkl, fij, fik, fil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                                param, dim, ncq, nparam, ntuples, potnum);
             
        if (decomp==0)
            cpuQuadrupletDecomposition(e, f, eijkl, fij, fik, fil, ai, aj, ak, al, ntuples, dim);
        else {
            cpuCenterAtomQuadrupletDecomposition(e, f, eijkl, fij, fik, fil, ilist, quadrupletnumsum, na, dim);   
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fij, jlist, bnumsum, index, naj, dim);
            Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fik, jlist, bnumsum, index, nak, dim);                
            Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fil, jlist, bnumsum, index, nal, dim);
        }
        
//         if (dim==2) {
//             if (decomp==0)
//                 cpuForceDecompositionQuadruplet2D(f, fij, fik, fil, ai, aj, ak, al, ntuples);
//             else {
//                 cpuAtomDecompositionQuadruplet2D(f, fij, fik, fil, ilist, quadrupletnumsum, na);   
//                 Int *jlist = &tmp.intmem[0];  // ntuples       
//                 Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
//                 Int *index = &tmp.intmem[2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fik, jlist, bnumsum, index, nak);                
//                 Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fil, jlist, bnumsum, index, nal);
//             }
//         }
//         else {
//             if (decomp==0)                  
//                 cpuForceDecompositionQuadruplet3D(f, fij, fik, fil, ai, aj, ak, al, ntuples);
//             else {
//                 cpuAtomDecompositionQuadruplet3D(f, fij, fik, fil, ilist, quadrupletnumsum, na);       
//                 Int *jlist = &tmp.intmem[0];  // ntuples       
//                 Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
//                 Int *index = &tmp.intmem[2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fik, jlist, bnumsum, index, nak);      
//                 Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fil, jlist, bnumsum, index, nal);
//             }
//         }            
    }        
}


// void implBondedTripletForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
//         dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
// {        
//     for (Int b=0; b<common.nba; b++) {
//         Int e1 = common.ablks[b];
//         Int e2 = common.ablks[b+1];            
//         Int na = e2 - e1; // number of atoms in this block
//         Int jnum = common.jnum;
//         Int jknum = jnum*(jnum-1)/2;
//         Int ncq = common.ncq;
//         Int dim = common.dim;
//                 
//         Int *ilist = &tmp.intmem[0]; //na     
//         if (typei>0) {               
//             Int *olist = &tmp.intmem[na]; //na        
//             cpuArrayFill(olist, e1, na);        
// 
//             Int *t0 = &tmp.intmem[2*na]; //na        
//             Int *t1 = &tmp.intmem[3*na]; //na        
//             na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
//         }
//         else {
//             cpuArrayFill(ilist, e1, na);        
//         }
//                 
//         Int *tripletnum = &tmp.intmem[na]; // na
//         Int *tripletlist = &tmp.intmem[2*na]; // na*jknum        
//         cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, 
//                 nb.neighnum, na, jnum, jknum, typej, typek, dim);                
//                 
//         //a list contains the starting positions of the first neighbor 
//         Int *tripletnumsum = &tmp.intmem[2*na + na*jknum]; // na+1        
//         cpuCumsum(tripletnumsum, tripletnum, na+1);                                         
//         int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
//                 
//         Int *ai = &tmp.intmem[1+3*na+na*jknum]; // ntuples        
//         Int *aj = &tmp.intmem[1+3*na+ntuples+na*jknum]; // ntuples   
//         Int *ak = &tmp.intmem[1+3*na+2*ntuples+na*jknum]; // ntuples   
//         Int *ti = &tmp.intmem[1+3*na+3*ntuples+na*jknum]; // ntuples        
//         Int *tj = &tmp.intmem[1+3*na+4*ntuples+na*jknum]; // ntuples        
//         Int *tk = &tmp.intmem[1+3*na+5*ntuples+na*jknum]; // ntuples        
//         dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
//         dstype *xik = &tmp.tmpmem[ntuples*dim]; // ntuples*dim
//         dstype *qi = &tmp.tmpmem[2*ntuples*dim]; // ntuples*ncq
//         dstype *qj = &tmp.tmpmem[ntuples*(2*dim+ncq)]; // ntuples*ncq
//         dstype *qk = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples*ncq
//         cpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
//                 ilist, nb.alist, atomtype, na, jnum, jknum, ncq, dim);       
//               
//         dstype *fij = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples*dim
//         dstype *fik = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*dim
//         cpuComputeTripletForce(fij, fik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
//                                 param, dim, ncq, nparam, ntuples, potnum);
//      
//         if (dim==2) {
//             if (decomp==0)
//                 cpuForceDecompositionTriplet2D(f, fij, fik, ai, aj, ak, ntuples);
//             else {
//                 cpuAtomDecompositionTriplet2D(f, fij, fik, ilist, tripletnumsum, na);                
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1+3*na+3*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1+3*na+4*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1+3*na+5*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1+3*na+6*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fik, jlist, bnumsum, index, nak);
//             }
//         }
//         else {
//             if (decomp==0)                  
//                 cpuForceDecompositionTriplet3D(f, fij, fik, ai, aj, ak, ntuples);
//             else {
//                 cpuAtomDecompositionTriplet3D(f, fij, fik, ilist, tripletnumsum, na);       
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1+3*na+3*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1+3*na+4*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1+3*na+5*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1+3*na+6*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fik, jlist, bnumsum, index, nak);                
//             }
//         }            
//     }        
// }
// 
// void implBondedQuadrupletForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
//         dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
// {        
//     for (Int b=0; b<common.nba; b++) {
//         Int e1 = common.ablks[b];
//         Int e2 = common.ablks[b+1];            
//         Int na = e2 - e1; // number of atoms in this block
//         Int jnum = common.jnum;
//         Int jklnum = jnum*(jnum-1)*(jnum-2)/6;
//         Int ncq = common.ncq;
//         Int dim = common.dim;
//                 
//         Int *ilist = &tmp.intmem[0]; //na     
//         if (typei>0) {               
//             Int *olist = &tmp.intmem[na]; //na        
//             cpuArrayFill(olist, e1, na);        
// 
//             Int *t0 = &tmp.intmem[2*na]; //na        
//             Int *t1 = &tmp.intmem[3*na]; //na        
//             na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
//         }
//         else {
//             cpuArrayFill(ilist, e1, na);        
//         }
//                 
//         Int *quadrupletnum = &tmp.intmem[na]; // na
//         Int *quadrupletlist = &tmp.intmem[2*na]; // na*jklnum                
//         cpuNeighQuadrupletList(quadrupletnum, quadrupletlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, 
//                 nb.neighnum, na, jnum, jklnum, typej, typek, typel, dim);                           
//         
//         //a list contains the starting positions of the first neighbor 
//         Int *quadrupletnumsum = &tmp.intmem[2*na + na*jklnum]; // na+1        
//         cpuCumsum(quadrupletnumsum, quadrupletnum, na+1);                                         
//         int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
//                 
//         Int *ai = &tmp.intmem[1+3*na+na*jklnum]; // ntuples        
//         Int *aj = &tmp.intmem[1+3*na+ntuples+na*jklnum]; // ntuples   
//         Int *ak = &tmp.intmem[1+3*na+2*ntuples+na*jklnum]; // ntuples   
//         Int *al = &tmp.intmem[1+3*na+3*ntuples+na*jklnum]; // ntuples   
//         Int *ti = &tmp.intmem[1+3*na+4*ntuples+na*jklnum]; // ntuples        
//         Int *tj = &tmp.intmem[1+3*na+5*ntuples+na*jklnum]; // ntuples        
//         Int *tk = &tmp.intmem[1+3*na+6*ntuples+na*jklnum]; // ntuples    
//         Int *tl = &tmp.intmem[1+3*na+7*ntuples+na*jklnum]; // ntuples    
//         dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
//         dstype *xik = &tmp.tmpmem[ntuples*dim]; // ntuples*dim
//         dstype *xil = &tmp.tmpmem[2*ntuples*dim]; // ntuples*dim
//         dstype *qi = &tmp.tmpmem[3*ntuples*dim]; // ntuples*ncq
//         dstype *qj = &tmp.tmpmem[ntuples*(3*dim+ncq)]; // ntuples*ncq
//         dstype *qk = &tmp.tmpmem[ntuples*(3*dim+2*ncq)]; // ntuples*ncq
//         dstype *ql = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*ncq
//         cpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
//                 quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, 
//                 jnum, jklnum, ncq, dim);       
//               
//         dstype *fij = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples*dim
//         dstype *fik = &tmp.tmpmem[ntuples*(5*dim+3*ncq)]; // ntuples*dim
//         dstype *fil = &tmp.tmpmem[ntuples*(6*dim+3*ncq)]; // ntuples*dim
//         cpuComputeQuadrupletForce(fij, fik, fil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
//                                 param, dim, ncq, nparam, ntuples, potnum);
//      
//         if (dim==2) {
//             if (decomp==0)
//                 cpuForceDecompositionQuadruplet2D(f, fij, fik, fil, ai, aj, ak, al, ntuples);
//             else {
//                 cpuAtomDecompositionQuadruplet2D(f, fij, fik, fil, ilist, quadrupletnumsum, na);                
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1+3*na+3*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1+3*na+4*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1+3*na+5*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1+3*na+6*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fik, jlist, bnumsum, index, nak);
//                 Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition2D(f, fil, jlist, bnumsum, index, nal);
//             }
//         }
//         else {
//             if (decomp==0)                  
//                 cpuForceDecompositionQuadruplet3D(f, fij, fik, fil, ai, aj, ak, al, ntuples);
//             else {
//                 cpuAtomDecompositionQuadruplet3D(f, fij, fik, fil, ilist, quadrupletnumsum, na);       
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ntuples]; // ntuples       
//                 Int *t0 = &tmp.intmem[1+3*na+3*ntuples]; // ntuples       
//                 Int *t1 = &tmp.intmem[1+3*na+4*ntuples]; // ntuples       
//                 Int *t2 = &tmp.intmem[1+3*na+5*ntuples]; // ntuples       
//                 Int *t3 = &tmp.intmem[1+3*na+6*ntuples]; // ntuples       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fik, jlist, bnumsum, index, nak);      
//                 Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
//                 cpuJAtomDecomposition3D(f, fil, jlist, bnumsum, index, nal);
//             }
//         }            
//      }        
// }
// 

// void implFullNeighPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
//         dstype* x, dstype *q, Int *atomtype, Int typei, Int decomp)
// {        
//     for (Int b=0; b<common.nba; b++) {
//         Int e1 = common.ablks[b];
//         Int e2 = common.ablks[b+1];            
//         Int na = e2 - e1; // number of atoms in this block
//         Int jnum = common.jnum;
//                 
//         Int *ilist = &tmp.intmem[0]; //na     
//         if (typei>0) {               
//             Int *olist = &tmp.intmem[na]; //na        
//             cpuArrayFill(olist, e1, na);        
// 
//             Int *t0 = &tmp.intmem[2*na]; //na        
//             Int *t1 = &tmp.intmem[3*na]; //na        
//             na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
//         }
//         else {
//             cpuArrayFill(ilist, e1, na);        
//         }
//         
//         Int *pairnum = &tmp.intmem[na]; // na
//         Int *pairlist = &tmp.intmem[2*na]; // na*jnum
//         cpuFullNeighPairList(pairnum, pairlist, ilist, nb.neighlist, nb.neighnum, na, jnum);
// 
//         //a list contains the starting positions of the first neighbor 
//         Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
//         cpuCumsum(pairnumsum, pairnum, na+1);                                         
//         int ijnum = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
//                 
//         Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ijnum        
//         Int *aj = &tmp.intmem[1+3*na+ijnum+na*jnum]; // ijnum        
//         Int *ti = &tmp.intmem[1+3*na+2*ijnum+na*jnum]; // ijnum        
//         Int *tj = &tmp.intmem[1+3*na+3*ijnum+na*jnum]; // ijnum        
//         dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
//         dstype *qi = &tmp.tmpmem[ijnum*dim]; // ijnum*ncq
//         dstype *qj = &tmp.tmpmem[ijnum*(dim+common.ncq)]; // ijnum*ncq
//         cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
//                 atomtype, na, common.jnum, common.ncq, dim);       
//                         
//         dstype *fij = &tmp.tmpmem[ijnum*(dim+2*common.ncq)]; // ijnum*dim
//         cpuComputeForce(fij, xij, qi, qj, ti, tj, ai, aj, app.pairparam, common.ncq, common.npairparam, ijnum);
// 
//         if (dim==2) {
//             if (decomp==0)
//                 cpuFullForceDecomposition2D(f, fij, ai, aj, ijnum);        
//             else
//                 cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);                        
//         }
//         else {
//             if (decomp==0)
//                 cpuFullForceDecomposition3D(f, fij, ai, aj, ijnum);    
//             else
//                 cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
//         }                
//     }        
// }

// void implHalfNeighPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
//         dstype* x, dstype *q, Int *atomtype, Int typei, Int decomp)
// {        
//     for (Int b=0; b<common.nba; b++) {
//         Int e1 = common.ablks[b];
//         Int e2 = common.ablks[b+1];            
//         Int na = e2 - e1; // number of atoms in this block
//                         
//         Int *ilist = &tmp.intmem[0]; //na     
//         if (typei>0) {               
//             Int *olist = &tmp.intmem[na]; //na        
//             cpuArrayFill(olist, e1, na);        
// 
//             Int *t0 = &tmp.intmem[2*na]; //na        
//             Int *t1 = &tmp.intmem[3*na]; //na        
//             na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
//         }
//         else {
//             cpuArrayFill(ilist, e1, na);        
//         }
//         
//         Int *anum = &tmp.intmem[na]; // na
//         cpuNeighCountHalfPair(anum, ilist, nb.alist, nb.neighlist, nb.neighnum, na, common.jnum);
//         
//         // a list contains the starting positions of the first neighbor 
//         Int *anumsum = &tmp.intmem[2*na]; // na+1        
//         cpuCumsum(anumsum, anum, na+1);                                         
//         int ijnum = IntArrayGetValueAtIndex(anumsum, na, common.backend);     
//         
//         Int *ai = &tmp.intmem[1+3*na]; // ijnum        
//         Int *aj = &tmp.intmem[1+3*na+ijnum]; // ijnum        
//         Int *ti = &tmp.intmem[1+3*na+2*ijnum]; // ijnum        
//         Int *tj = &tmp.intmem[1+3*na+3*ijnum]; // ijnum        
//         dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
//         dstype *qi = &tmp.tmpmem[ijnum*dim]; // ijnum*ncq
//         dstype *qj = &tmp.tmpmem[ijnum*(dim+common.ncq)]; // ijnum*ncq
//         cpuGetHalfNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, ilist, nb.alist, nb.neighlist, nb.neighnum, 
//                 atomtype, na, common.jnum, common.ncq, dim);                   
// 
//         dstype *fij = &tmp.tmpmem[ijnum*(dim+2*common.ncq)]; // ijnum*dim
//         cpuComputeForce(fij, xij, qi, qj, ti, tj, ai, aj, app.pairparam, common.ncq, common.npairparam, ijnum);
// 
//         if (dim==2) {
//             if (decomp==0)
//                 cpuHalfForceDecomposition2D(f, fij, ai, aj, ijnum);       
//             else {
//                 cpuIAtomDecomposition2D(f, fij, ilist, anumsum, na);      
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
//                 Int *p0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
//                 Int *p1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
//                 Int *p2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
//                 Int *p3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, ijnum);
//                 cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
//             }
//         }
//         else {
//             if (decomp==0)
//                 cpuHalfForceDecomposition3D(f, fij, ai, aj, ijnum);        
//             else {
//                 cpuIAtomDecomposition3D(f, fij, ilist, anumsum, na);      
//                 Int *jlist = &tmp.intmem[0]; //2*na  
//                 Int *bnumsum = &tmp.intmem[2*na]; //2*na  
//                 Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
//                 Int *p0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
//                 Int *p1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
//                 Int *p2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
//                 Int *p3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
//                 Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, ijnum);
//                 cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                   
//             }
//         }
//     }        
// }
//       
//  
// void implFullNeighPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
//         dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int decomp)
// {        
//     for (Int b=0; b<common.nba; b++) {
//         Int e1 = common.ablks[b];
//         Int e2 = common.ablks[b+1];            
//         Int na = e2 - e1; // number of atoms in this block
//                 
//         Int *ilist = &tmp.intmem[0]; //na     
//         if (typei>0) {               
//             Int *olist = &tmp.intmem[na]; //na        
//             cpuArrayFill(olist, e1, na);        
// 
//             Int *t0 = &tmp.intmem[2*na]; //na        
//             Int *t1 = &tmp.intmem[3*na]; //na        
//             na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
//         }
//         else {
//             cpuArrayFill(ilist, e1, na);        
//         }
//         
//         Int *anum = &tmp.intmem[na]; // na
//         cpuNeighPairList(anum, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, nb.jnum, typej);        
//         
//         // a list contains the starting positions of the first neighbor 
//         Int *anumsum = &tmp.intmem[2*na]; // na+1        
//         cpuCumsum(anumsum, anum, na+1);                                         
//         int ijnum = IntArrayGetValueAtIndex(anumsum, na, common.backend);     
//         
//         Int *ai = &tmp.intmem[1+3*na]; // ijnum        
//         Int *aj = &tmp.intmem[1+3*na+ijnum]; // ijnum        
//         Int *ti = &tmp.intmem[1+3*na+2*ijnum]; // ijnum        
//         Int *tj = &tmp.intmem[1+3*na+3*ijnum]; // ijnum        
//         dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
//         dstype *qi = &tmp.tmpmem[ijnum*dim]; // ijnum*ncq
//         dstype *qj = &tmp.tmpmem[ijnum*(dim+common.ncq)]; // ijnum*ncq
//         cpuGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, ilist, nb.alist, nb.neighlist, nb.neighnum, 
//                 atomtype, na, common.jnum, common.ncq, dim);                   
// 
//         dstype *fij = &tmp.tmpmem[ijnum*(dim+2*common.ncq)]; // ijnum*dim
//         cpuComputeForce(fij, xij, qi, qj, ti, tj, ai, aj, app.pairparam, common.ncq, common.npairparam, ijnum);
// 
//         if (dim==2) {
//             if (decomp==0)
//                 cpuFullForceDecomposition2D(f, fij, ai, aj, ijnum);        
//             else
//                 cpuIAtomDecomposition2D(f, fij, ilist, anumsum, na);                        
//         }
//         else {
//             if (decomp==0)
//                 cpuFullForceDecomposition3D(f, fij, ai, aj, ijnum);    
//             else
//                 cpuIAtomDecomposition3D(f, fij, ilist, anumsum, na);     
//         }                
//     }        
// }
// 

#endif

