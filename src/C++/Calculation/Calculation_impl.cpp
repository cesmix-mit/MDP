#ifndef __CALCULATIONIMPL
#define __CALCULATIONIMPL

// rcutmax: maximum cut-off radius

// NONBONDED POTENTIALS
// Single potential     : 1 nonbonded function 
// Pair potential       : p nonbonded functions and p cut-off radius

// BONDED POTENTIALS
// Single potential     : p bonded functions
// Pair potential       : p bonded functions and p cut-off radii
// triplet potential    : p bonded functions and p cut-off radii
// quadruplet potential : p bonded functions and p cut-off radii

// BOND ORDER POTENTIALS
// eam potential        : 3*p bonded functions and p cut-off radii (two-body interactions)  
// tersoff potential    : 3*p bonded functions and p cut-off radii (three-body interactions)


void cpuComputeSingleEnergy(dstype *ei, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}
void cpuComputeSingleForce(dstype *fi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}
void cpuComputeSingleParam(dstype *gi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}
void cpuComputeSingleEnergyForce(dstype *ei, dstype *fi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}
void cpuComputeSingle(dstype *ei, dstype *fi, dstype *gi, dstype *xi, dstype *qi, int *ti, 
        int *ai, dstype *singleparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}

void cpuComputePairForce(dstype *fij, dstype *xij, dstype *qi, dstype *qj, int *ti, int *tj, 
        int *ai, int *aj, dstype *pairparam, int dim, int ncq, int nparam, int ijnum, int potnum)
{
}

void cpuComputeTripletForce(dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
     int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *tripletparam, int dim, int ncq, int nparam, int ijknum, int potnum)
{
}

void implSingleForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            Int *t1 = &tmp.intmem[2*na]; //na        
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
        cpuComputeSingleForce(fi, xi, qi, ti, ai, param, dim, ncq, nparam, na, potnum);

        if (dim==2)
            cpuSingleDecomposition2D(f, fi, ai, na);
        else
            cpuSingleDecomposition2D(f, fi, ai, na);        
    }        
}

void implNonbondedSingleForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int potnum)
{    
    implSingleForce(f, nb, common, app, tmp, x, q, param, atomtype, nparam, 0, potnum);            
}

void implBondedSingleForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{
    implSingleForce(f, nb, common, app, tmp, x, q, param, atomtype, nparam, typei, potnum);            
}

void implFullNeighPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            Int *t1 = &tmp.intmem[2*na]; //na        
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
        int ijnum = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ijnum        
        Int *aj = &tmp.intmem[1+3*na+ijnum+na*jnum]; // ijnum        
        Int *ti = &tmp.intmem[1+3*na+2*ijnum+na*jnum]; // ijnum        
        Int *tj = &tmp.intmem[1+3*na+3*ijnum+na*jnum]; // ijnum        
        dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
        dstype *qi = &tmp.tmpmem[ijnum*dim]; // ijnum*ncq
        dstype *qj = &tmp.tmpmem[ijnum*(dim+ncq)]; // ijnum*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ijnum*(dim+2*ncq)]; // ijnum*dim
        cpuComputePairForce(fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, ijnum, potnum);

        if (dim==2) {
            if (decomp==0)
                cpuFullForceDecomposition2D(f, fij, ai, aj, ijnum);        
            else
                cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);                        
        }
        else {
            if (decomp==0)
                cpuFullForceDecomposition3D(f, fij, ai, aj, ijnum);    
            else
                cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
        }                
    }        
}

void implFullNeighNonbondedPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int decomp, Int potnum)
{    
    implFullNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, 0, 0, decomp, potnum);            
}

void implFullNeighBondedPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{
    if (typei==typej) {
        implFullNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);        
    }
    else {
        implFullNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);   
        implFullNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typej, typei, decomp, potnum);        
    }
}


void implHalfNeighPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            Int *t1 = &tmp.intmem[2*na]; //na        
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
        int ijnum = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ijnum        
        Int *aj = &tmp.intmem[1+3*na+ijnum+na*jnum]; // ijnum        
        Int *ti = &tmp.intmem[1+3*na+2*ijnum+na*jnum]; // ijnum        
        Int *tj = &tmp.intmem[1+3*na+3*ijnum+na*jnum]; // ijnum        
        dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
        dstype *qi = &tmp.tmpmem[ijnum*dim]; // ijnum*ncq
        dstype *qj = &tmp.tmpmem[ijnum*(dim+ncq)]; // ijnum*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ijnum*(dim+2*ncq)]; // ijnum*dim
        cpuComputePairForce(fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, ijnum, potnum);

        if (dim==2) {
            if (decomp==0)
                cpuHalfForceDecomposition2D(f, fij, ai, aj, ijnum);       
            else {
                cpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
                Int *jlist = &tmp.intmem[0]; //2*na  
                Int *bnumsum = &tmp.intmem[2*na]; //2*na  
                Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
                Int *p0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
                Int *p1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
                Int *p2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
                Int *p3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, ijnum);
                cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
            }
        }
        else {
            if (decomp==0)
                cpuHalfForceDecomposition3D(f, fij, ai, aj, ijnum);        
            else {
                cpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);      
                Int *jlist = &tmp.intmem[0]; //2*na  
                Int *bnumsum = &tmp.intmem[2*na]; //2*na  
                Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
                Int *p0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
                Int *p1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
                Int *p2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
                Int *p3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, ijnum);
                cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                   
            }
        }            
    }        
}

void implHalfNeighNonbondedPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int decomp, Int potnum)
{    
    implHalfNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, 0, 0, decomp, potnum);            
}

void implHalfNeighBondedPairForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{
    if (typei==typej) {
        implHalfNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);        
    }
    else {
        implHalfNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typei, typej, decomp, potnum);   
        implHalfNeighPairForce(f, nb, common, app, tmp, x, q, param, rcutsq, atomtype, nparam, typej, typei, decomp, potnum);        
    }
}


void implBondedTripletForce(dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int jknum = jnum*(jnum-1)/2;
        Int ncq = common.ncq;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[2*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
                
        Int *tripletnum = &tmp.intmem[na]; // na
        Int *tripletlist = &tmp.intmem[2*na]; // na*jknum        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, jknum, typej, typek, dim);                
                
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na + na*jknum]; // na+1        
        cpuCumsum(tripletnumsum, tripletnum, na+1);                                         
        int ijnum = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jknum]; // ijnum        
        Int *aj = &tmp.intmem[1+3*na+ijnum+na*jknum]; // ijnum   
        Int *ak = &tmp.intmem[1+3*na+2*ijnum+na*jknum]; // ijnum   
        Int *ti = &tmp.intmem[1+3*na+3*ijnum+na*jknum]; // ijnum        
        Int *tj = &tmp.intmem[1+3*na+4*ijnum+na*jknum]; // ijnum        
        Int *tk = &tmp.intmem[1+3*na+5*ijnum+na*jknum]; // ijnum        
        dstype *xij = &tmp.tmpmem[0]; // ijnum*dim
        dstype *xik = &tmp.tmpmem[ijnum*dim]; // ijnum*dim
        dstype *qi = &tmp.tmpmem[2*ijnum*dim]; // ijnum*ncq
        dstype *qj = &tmp.tmpmem[ijnum*(2*dim+ncq)]; // ijnum*ncq
        dstype *qk = &tmp.tmpmem[ijnum*(2*dim+2*ncq)]; // ijnum*ncq
        cpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, jnum, jknum, ncq, dim);       
              
        dstype *fij = &tmp.tmpmem[ijnum*(2*dim+3*ncq)]; // ijnum*dim
        dstype *fik = &tmp.tmpmem[ijnum*(3*dim+3*ncq)]; // ijnum*dim
        cpuComputeTripletForce(fij, fik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                                param, dim, ncq, nparam, ijnum, potnum);
     
        if (dim==2) {
            if (decomp==0)
                cpuForceDecompositionTriplet2D(f, fij, fik, ai, aj, ak, ijnum);
            else {
                cpuAtomDecompositionTriplet2D(f, fij, fik, ilist, tripletnumsum, na);                
                Int *jlist = &tmp.intmem[0]; //2*na  
                Int *bnumsum = &tmp.intmem[2*na]; //2*na  
                Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
                Int *t0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
                Int *t1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
                Int *t2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
                Int *t3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ijnum);
                cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);
                Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ijnum);
                cpuJAtomDecomposition2D(f, fik, jlist, bnumsum, index, nak);
            }
        }
        else {
            if (decomp==0)                  
                cpuForceDecompositionTriplet3D(f, fij, fik, ai, aj, ak, ijnum);
            else {
                cpuAtomDecompositionTriplet3D(f, fij, fik, ilist, tripletnumsum, na);       
                Int *jlist = &tmp.intmem[0]; //2*na  
                Int *bnumsum = &tmp.intmem[2*na]; //2*na  
                Int *index = &tmp.intmem[1+3*na+2*ijnum]; // ijnum       
                Int *t0 = &tmp.intmem[1+3*na+3*ijnum]; // ijnum       
                Int *t1 = &tmp.intmem[1+3*na+4*ijnum]; // ijnum       
                Int *t2 = &tmp.intmem[1+3*na+5*ijnum]; // ijnum       
                Int *t3 = &tmp.intmem[1+3*na+6*ijnum]; // ijnum       
                Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ijnum);
                cpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);
                Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ijnum);
                cpuJAtomDecomposition3D(f, fik, jlist, bnumsum, index, nak);                
            }
        }            
    }        
}


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
//             Int *t1 = &tmp.intmem[2*na]; //na        
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
//             Int *t1 = &tmp.intmem[2*na]; //na        
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
//             Int *t1 = &tmp.intmem[2*na]; //na        
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

