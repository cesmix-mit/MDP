/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CPUPOTENTIALS
#define __CPUPOTENTIALS

#include "opuApp.cpp" 
       
void cpuElectronDensity(dstype *rhoi, dstype *rhoij, int *pairnum, int *pairnumsum, int inum) 
{
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        rhoi[i] = 0.0;
        for (int j=0; j<jnum; j++) 
            rhoi[i] += rhoij[start+j];        
    }
}

void cpuEmbedingForce(dstype *fij, dstype *d_rhoi, int *pairnum, int *pairnumsum, int inum)
{    
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        for (int j=0; j<jnum; j++)             
            fij[start+j] = d_rhoi[i]*fij[start+j];                
    }        
}

void cpuSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        dstype *ei = &tmp.tmpmem[na*(2*dim+ncq)]; // na
        dstype *du = &tmp.tmpmem[na*(2*dim+ncq+1)]; // na
        cpuComputeSingleEnergyForce(ei, du, fi, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        cpuSingleDecomposition(e, f, ei, fi, ai, na, dim);        
    }        
}
void cpuNonbondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        cpuSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void cpuBondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        cpuSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void cpuFullNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
        //a list contains the starting positions of the first neighbor              
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
        
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
        
        //printArray2D(f, dim, 10, common.backend);        
        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
                        
        if (decomp==0) // force decomposition
            cpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else // atom decomposition
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);     
                
#ifdef HAVE_DEBUG                      
        writearray2file("xijcpu.bin", xij, ntuples*dim, common.backend); 
        writearray2file("eijcpu.bin", eij, ntuples, common.backend); 
        writearray2file("fijcpu.bin", fij, ntuples*dim, common.backend); 
        writearray2file("ecpu.bin", e, common.inum, common.backend); 
        writearray2file("fcpu.bin", f, common.inum*dim, common.backend); 
#endif                                                        
    }        
}

void cpuHalfNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0) // force decomposition
            cpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, ntuples, dim);
        else { // atom decomposition
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
    }        
}


void cpuNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            cpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            cpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void cpuBondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                cpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                cpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                cpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                cpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void cpuBO2EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *du = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq+1)]; // ntuples        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq+2)]; // ntuples*dim                   
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        dstype *gi = &tmp.tmpmem[2*na]; // na
        dstype *du2 = &tmp.tmpmem[3*na]; // na
        cpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        cpuPaircDensityGradient(ei, du2, gi, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);
        cpuEmbedingForce(fij, gi, pairnum, pairnumsum, na);        
                
        cpuPutArrayAtIndex(e, ei, ilist, na);      
        if (decomp==0) // force decomposition
            cpuHalfForceDecomposition(f, fij, ai, aj, dim, ntuples);        
        else {   // atom decomposition             
            cpuIAtomDecomposition(f, fij, ilist, pairnumsum, dim, na);      
            cpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
            cpuJAtomDecomposition(f, fij, jlist, bnumsum, index, dim, naj);                
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
void cpuBoPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        cpuBO2EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void cpuTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax        
        if ((typej>0) && (typek>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        cpuTripletnum(tripletnum, pairnum, na);
        //for (int ii=0; ii<na; ii++)
        //    tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       

        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        //Cumsum(tripletnumsum, tripletnum, na+1, backend);                                         
        //int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, backend);                             
        
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        cpuNeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        cpuArrayCopy(tripletlist, temp, 2*ntuples);
        
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
        dstype *eijk = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(4*dim+3*ncq+1)]; // ntuples
        
        cpuComputeTripletEnergyForce(eijk, du, fij, fik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                         param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);             
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
    }        
}
void cpuNonbondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        cpuTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void cpuBondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        cpuTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void cpuBO3EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax       
        if (typej == typei) {
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        
        }
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1    
        //Cumsum(pairnumsum, pairnum, na+1, backend);                                         
        //int npairs = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int npairs = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // npairs    
        Int *aj = &tmp.intmem[1+3*na+npairs+na*neighmax]; // npairs        
        Int *ti = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs        
        Int *tj = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs        
        dstype *eij = &tmp.tmpmem[0]; // npairs
        dstype *fij = &tmp.tmpmem[npairs]; // npairs*dim
        dstype *xij = &tmp.tmpmem[npairs*(1+dim)]; // npairs*dim
        dstype *qi = &tmp.tmpmem[npairs*(1+2*dim)]; // npairs*ncq
        dstype *qj = &tmp.tmpmem[npairs*(1+2*dim+ncq)]; // npairs*ncq
        dstype *du = &tmp.tmpmem[npairs*(1+2*dim+2*ncq)]; // npairs
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                                
        //cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, npairs, potnum);
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs*neighmax        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, neighmax, typek, dim);                
        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[1+3*na+3*npairs+(na+npairs)*neighmax]; // npairs+1        
        //Cumsum(tripletnumsum, tripletnum, npairs+1);                                 
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax], &tmp.intmem[3+3*na+5*npairs+(na+npairs)*neighmax], npairs+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, npairs, common.backend);     
                
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj, tripletnum, tripletnumsum, tripletlist)         
        Int *a3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax]; // ntuples        
        Int *a3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+ntuples]; // ntuples   
        Int *a3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+2*ntuples]; // ntuples   
        Int *t3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+3*ntuples]; // ntuples        
        Int *t3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+4*ntuples]; // ntuples        
        Int *t3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+5*ntuples]; // ntuples        
        dstype *x3ij = &tmp.tmpmem[npairs*(1+dim)]; // ntuples*dim
        dstype *x3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*dim]; // ntuples*dim
        dstype *q3i = &tmp.tmpmem[npairs*(1+dim)+2*ntuples*dim]; // ntuples*ncq
        dstype *q3j = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+ncq)]; // ntuples*ncq
        dstype *q3k = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+2*ncq)]; // ntuples*ncq       
        
        cpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, neighmax, npairs, ncq, dim);                      
        
        dstype *f3ij = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *f3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *e3ijk = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du2 = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq+1)]; // ntuples
        cpuComputeTripletEnergyForce(e3ijk, du2, f3ij, f3ik, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *h3ij = &tmp.tmpmem[npairs*(1+dim)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2+dim)]; // npairs
        dstype *d3ij = &tmp.tmpmem[npairs*(3+dim)]; // npairs
        dstype *g3ij = &tmp.tmpmem[npairs*(4+dim)]; // dim*npairs
        
        cpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        cpuElectronDensity(g3ij, f3ij, tripletnum, tripletnumsum, npairs);
        cpuTripletcDensityGradient(c3ij, du, d3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);
        for (int i=0; i<npairs; i++) 
            for (int j=0; j<dim; j++)
                fij[j+3*i] = fij[j+3*i]*c3ij[i] + eij[i]*d3ij[i]*g3ij[j+3*i];   
        
        for (int i=0; i<npairs; i++) 
            eij[i] = eij[i]*c3ij[i];        
                    
        // pairnum, pairnumsum, pairlist, tripletnum, tripletlist, tripletnumsum, t3i, t3j, t3k -> ai, aj, a3i, a3j, a3k
//         for (int i=0; i<npairs; i++) {
//             tmp.intmem[i] = tmp.intmem[1+3*na+na*neighmax+i]; // ai
//             tmp.intmem[npairs+i] = tmp.intmem[1+3*na+npairs+na*neighmax+i]; // aj
//         }
//         for (int i=0; i<ntuples; i++) {
//             tmp.intmem[2*npairs+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+i]; //a3i
//             tmp.intmem[2*npairs+ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+ntuples+i]; //a3j
//             tmp.intmem[2*npairs+2*ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+2*ntuples+i]; //a3k
//         }
        cpuArrayCopy(&tmp.intmem[0], &tmp.intmem[1+3*na+na*neighmax], npairs);
        cpuArrayCopy(&tmp.intmem[npairs], &tmp.intmem[1+3*na+na*neighmax+npairs], npairs);
        cpuArrayCopy(&tmp.intmem[2*npairs], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax], ntuples);
        cpuArrayCopy(&tmp.intmem[2*npairs+ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+ntuples], ntuples);
        cpuArrayCopy(&tmp.intmem[2*npairs+2*ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+2*ntuples], ntuples);

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
        
        for (int i=0; i<npairs; i++) {
            int m = tripletnum[i];
            int s = tripletnumsum[i];
            for (int j=0; j<m; j++) 
                for (int d=0; d<dim; d++)
                    f3ik[d+3*(s+j)] = eij[i]*d3ij[i]*f3ik[d+3*(s+j)];                                                
        }                        

        if (decomp==0)
            cpuHalfForceDecomposition(f, f3ik, a3i, a3k, dim, ntuples);        
        else {
            Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
            Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
            Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
            Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
            Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
            Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
            Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
            Int nai = cpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
            cpuJAtomDecomposition(f, f3ik, jlist, bnumsum, index, dim, nai);                            
            Int nak = cpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
            cpuJAtomDecomposition(f, f3ik, jlist, bnumsum, index, dim, nak);          
        }

//         if (dim==2) {            
//             if (decomp==0)
//                 cpuHalfForceDecomposition2D(f, f3ik, a3i, a3k, ntuples);        
//             else {
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
//                 Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
//                 Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
//                 Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
//                 Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
//                 Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
//                 Int nai = cpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
//                 cpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nai);                            
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
//                 cpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nak);          
//             }
//         }
//         else {
//             if (decomp==0)
//                 cpuHalfForceDecomposition3D(f, f3ik, a3i, a3k, ntuples);    
//             else {
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
//                 Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
//                 Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
//                 Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
//                 Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
//                 Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
//                 Int nai = cpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
//                 cpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nai);                            
//                 Int nak = cpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
//                 cpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nak);         
//             }
//         }        
    }        
}

void cpuBoTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        cpuBO3EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void cpuQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax    
        if ((typej>0) && (typek>0) && (typel>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, typel, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        //for (int ii=0; ii<na; ii++)
        //    quadrupletnum[ii] = (pairnum[ii]-2)*(pairnum[ii]-1)*pairnum[ii]/6;                                       
        cpuQuadrupletnum(quadrupletnum, pairnum, na);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        //Cumsum(quadrupletnumsum, quadrupletnum, na+1, backend);       
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                        
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        cpuNeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        cpuArrayCopy(quadrupletlist, temp, 3*ntuples);
        
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
        dstype *du = &tmp.tmpmem[ntuples*(7*dim+3*ncq+1)]; // ntuples
        
        cpuComputeQuadrupletEnergyForce(eijkl, du, fij, fik, fil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                                

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
    }        
}
void cpuNonbondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        cpuQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void cpuBondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        cpuQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void cpuEmpiricalPotentialEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.npot1a > 0)
        cpuNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    if (common.npot1b > 0)
        cpuBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    if (common.npot2a > 0)
        cpuNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    if (common.npot2b > 0)
        cpuBondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    if (common.npot2c > 0)    
        cpuBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    if (common.npot3a > 0)    
        cpuNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    if (common.npot3b > 0)
        cpuBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    if (common.npot3c > 0)    
        cpuBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    if (common.npot4a > 0)
        cpuNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    if (common.npot4b > 0)
        cpuBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);              
    
    ArrayMinus(f, f, common.dim*common.inum, common.backend);
}

void cpuEmpiricalPotentialDescriptors(dstype *ei, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;
    
    int m = 0;
    int n = 0;
    if (common.npot1a > 0)
        cpuNonbondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    n += common.npot1a*dim*inum;
    if (common.npot1b > 0)
        cpuBondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    n += common.npot1b*dim*inum;
    if (common.npot2a > 0)
        cpuNonbondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    n += common.npot2a*dim*inum;
    if (common.npot2b > 0)
        cpuBondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    n += common.npot2b*dim*inum;
    if (common.npot2c > 0)    
        cpuBoPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    m += common.npot2c*inum;
    n += common.npot2c*dim*inum;
    if (common.npot3a > 0)    
        cpuNonbondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    n += common.npot3a*dim*inum;
    if (common.npot3b > 0)
        cpuBondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    n += common.npot3b*dim*inum;
    if (common.npot3c > 0)    
        cpuBoTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    n += common.npot3c*dim*inum;
    if (common.npot4a > 0)
        cpuNonbondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    n += common.npot4a*dim*inum;
    if (common.npot4b > 0)
        cpuBondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                                
}

void cpuSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {        
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        
        dstype *ei = &tmp.tmpmem[na*(dim+ncq)]; // na
        cpuSingle(ei, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        cpuSingleDecomposition(e, ei, ai, na, dim);        
    }        
}
void cpuNonbondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        cpuSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void cpuBondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        cpuSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void cpuFullNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);        
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim        
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                
        if (decomp==0)
            cpuFullNeighPairDecomposition(e, eij, ai, ntuples, dim);
        else
            cpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);                        
    }        
}

void cpuHalfNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);           
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            cpuHalfNeighPairDecomposition(e, eij, ai, aj, ntuples, dim);
        else {
            cpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);   
            cpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);            
            cpuNeighborAtomPairDecomposition(e, eij, jlist, bnumsum, index, naj, dim);        
        }                
    }        
}
void cpuNonbondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            cpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            cpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void cpuBondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                cpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                cpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                cpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                cpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void cpuBO2Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples    
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        cpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        opuPaircDensity(ei, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);                
        cpuPutArrayAtIndex(e, ei, ilist, na);        
    }        
}
void cpuBoPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        cpuBO2Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void cpuTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax        
        if ((typej>0) && (typek>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        cpuTripletnum(tripletnum, pairnum, na);
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                                              
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
                
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        cpuNeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        cpuArrayCopy(tripletlist, temp, 2*ntuples);
        
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
              
        dstype *eijk = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples
        
        cpuTriplet(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                         param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);             
        if (decomp==0)            
            cpuTripletDecomposition(e, eijk, ai, aj, ak, ntuples, dim);
        else {
            cpuCenterAtomTripletDecomposition(e,  eijk, ilist, tripletnumsum, na, dim);                
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            cpuNeighborAtomTripletDecomposition(e, eijk, jlist, bnumsum, index, naj, dim);
            Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            cpuNeighborAtomTripletDecomposition(e, eijk, jlist, bnumsum, index, nak, dim);
        }    
    }        
}
void cpuNonbondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        cpuTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void cpuBondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        cpuTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void cpuBO3Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax       
        if (typej == typei) {
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        
        }
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1            
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int npairs = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // npairs    
        Int *aj = &tmp.intmem[1+3*na+npairs+na*neighmax]; // npairs        
        Int *ti = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs        
        Int *tj = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs        
        dstype *eij = &tmp.tmpmem[0]; // npairs
        dstype *xij = &tmp.tmpmem[npairs*(1)]; // npairs*dim
        dstype *qi = &tmp.tmpmem[npairs*(1+dim)]; // npairs*ncq
        dstype *qj = &tmp.tmpmem[npairs*(1+dim+ncq)]; // npairs*ncq
        dstype *du = &tmp.tmpmem[npairs*(1+dim+2*ncq)]; // npairs
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       

         
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);        
                
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs*neighmax        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, neighmax, typek, dim);                
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[1+3*na+3*npairs+(na+npairs)*neighmax]; // npairs+1                 
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax], &tmp.intmem[3+3*na+5*npairs+(na+npairs)*neighmax], npairs+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, npairs, common.backend);     
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj, tripletnum, tripletnumsum, tripletlist)         
        Int *a3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax]; // ntuples        
        Int *a3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+ntuples]; // ntuples   
        Int *a3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+2*ntuples]; // ntuples   
        Int *t3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+3*ntuples]; // ntuples        
        Int *t3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+4*ntuples]; // ntuples        
        Int *t3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+5*ntuples]; // ntuples        
        dstype *x3ij = &tmp.tmpmem[npairs*(1)]; // ntuples*dim
        dstype *x3ik = &tmp.tmpmem[npairs*(1)+ntuples*dim]; // ntuples*dim
        dstype *q3i = &tmp.tmpmem[npairs*(1)+2*ntuples*dim]; // ntuples*ncq
        dstype *q3j = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+ncq)]; // ntuples*ncq
        dstype *q3k = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+2*ncq)]; // ntuples*ncq
        
        cpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, neighmax, npairs, ncq, dim);          
        dstype *e3ijk = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+3*ncq)]; // ntuples
        cpuTriplet(e3ijk, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                
        dstype *h3ij = &tmp.tmpmem[npairs*(1)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2)]; // npairs
        cpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        opuTripletcDensity(c3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);
                
        for (int i=0; i<npairs; i++) 
            eij[i] = eij[i]*c3ij[i];
                    
        // pairnum, pairnumsum, pairlist, tripletnum, tripletlist, tripletnumsum, t3i, t3j, t3k -> ai, aj, a3i, a3j, a3k
        for (int i=0; i<npairs; i++) {
            tmp.intmem[i] = tmp.intmem[1+3*na+na*neighmax+i]; // ai
            tmp.intmem[npairs+i] = tmp.intmem[1+3*na+npairs+na*neighmax+i]; // aj
        }
                
        if (decomp==0)
            cpuHalfNeighPairDecomposition(e, eij, ai, aj, npairs, dim);
        else {
            cpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);   
            Int *jlist = &tmp.intmem[2*npairs]; //npairs  
            Int *bnumsum = &tmp.intmem[3*npairs]; //npairs  
            Int *index = &tmp.intmem[4*npairs]; // npairs
            Int *p0 = &tmp.intmem[5*npairs]; // npairs       
            Int *p1 = &tmp.intmem[6*npairs]; // npairs       
            Int *p2 = &tmp.intmem[7*npairs]; // npairs       
            Int *p3 = &tmp.intmem[8*npairs]; // npairs                   
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
            cpuNeighborAtomPairDecomposition(e, eij, jlist, bnumsum, index, naj, dim);                 
        }                 
    }        
}
void cpuBoTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        cpuBO3Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void cpuQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax    
        if ((typej>0) && (typek>0) && (typel>0))
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, typel, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        cpuQuadrupletnum(quadrupletnum, pairnum, na);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                                
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        cpuNeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        cpuArrayCopy(quadrupletlist, temp, 3*ntuples);
        
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
        
        dstype *eijkl = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
    
        cpuQuadruplet(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                                
        if (decomp==0)
            cpuQuadrupletDecomposition(e, eijkl,  ai, aj, ak, al, ntuples, dim);
        else {
            cpuCenterAtomQuadrupletDecomposition(e, eijkl, ilist, quadrupletnumsum, na, dim);   
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, naj, dim);
            Int nak = cpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, nak, dim);                
            Int nal = cpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
            cpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, nal, dim);
        }        
    }        
}
void cpuNonbondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        cpuQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void cpuBondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        cpuQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void cpuEmpiricalPotentialDescriptors(dstype *ei, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;        
            
    int m = 0;
    if (common.npot1a > 0)
        cpuNonbondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    if (common.npot1b > 0)
        cpuBondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    if (common.npot2a > 0)
        cpuNonbondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    if (common.npot2b > 0)
        cpuBondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    if (common.npot2c > 0)    
        cpuBoPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
 
    m += common.npot2c*inum;
    if (common.npot3a > 0)    
        cpuNonbondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    if (common.npot3b > 0)
        cpuBondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    if (common.npot3c > 0)    
        cpuBoTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    if (common.npot4a > 0)
        cpuNonbondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    if (common.npot4b > 0)
        cpuBondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                            
}

void cpuFullNeighPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
        //a list contains the starting positions of the first neighbor              
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
        
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
        
        //printArray2D(f, dim, 10, common.backend);        
        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
                        
        if (decomp==0) {// force decomposition
            cpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);            
        } else {// atom decomposition
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);     
        }
                
        // virial tally
        cpuVirialPairTally(v, fij, xij, -0.5, ai, dim, common.inum, ntuples);        
                
#ifdef HAVE_DEBUG                      
        writearray2file("xijcpu.bin", xij, ntuples*dim, common.backend); 
        writearray2file("eijcpu.bin", eij, ntuples, common.backend); 
        writearray2file("fijcpu.bin", fij, ntuples*dim, common.backend); 
        writearray2file("ecpu.bin", e, common.inum, common.backend); 
        writearray2file("fcpu.bin", f, common.inum*dim, common.backend); 
#endif                                                        
    }        
}

void cpuHalfNeighPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
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
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            cpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0) // force decomposition
            cpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, ntuples, dim);
        else { // atom decomposition
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
        
        // virial tally
        cpuVirialPairTally(v, fij, xij, -0.5, ai, aj, dim, common.inum, ntuples);                
    }        
}

void cpuNonbondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            cpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            cpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void cpuBondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                cpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                cpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                cpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                cpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                cpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void cpuEmpiricalPotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.npot1a > 0)
        cpuNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    if (common.npot1b > 0)
        cpuBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    if (common.npot2a > 0)
        cpuNonbondedPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    if (common.npot2b > 0)
        cpuBondedPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              
    
    
//     if (common.npot2c > 0)    
//         cpuBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
// 
//     if (common.npot3a > 0)    
//         cpuNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
//     
//     if (common.npot3b > 0)
//         cpuBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              
// 
//     if (common.npot3c > 0)    
//         cpuBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              
// 
//     if (common.npot4a > 0)
//         cpuNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              
// 
//     if (common.npot4b > 0)
//         cpuBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);              
    
    ArrayMinus(f, f, common.dim*common.inum, common.backend);
}

void cpuEmpiricalPotentialDescriptors(dstype *ei, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;
    
    int m = 0;
    int n = 0;
    int k = 0;
    if (common.npot1a > 0)
        cpuNonbondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    n += common.npot1a*dim*inum;
    if (common.npot1b > 0)
        cpuBondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    n += common.npot1b*dim*inum;    
    if (common.npot2a > 0)
        cpuNonbondedPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    n += common.npot2a*dim*inum;
    k += common.npot2a*6*inum;
    if (common.npot2b > 0)
        cpuBondedPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

//     m += common.npot2b*inum;
//     n += common.npot2b*dim*inum;
//    k += common.npot2b*6*inum;
//     if (common.npot2c > 0)    
//         cpuBoPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
// 
//     m += common.npot2c*inum;
//     n += common.npot2c*dim*inum;
//    k += common.npot2c*6*inum;
//     if (common.npot3a > 0)    
//         cpuNonbondedTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
//     
//     m += common.npot3a*inum;
//     n += common.npot3a*dim*inum;
//    k += common.npot3a*6*inum;
//     if (common.npot3b > 0)
//         cpuBondedTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              
// 
//     m += common.npot3b*inum;
//     n += common.npot3b*dim*inum;
//    k += common.npot3b*6*inum;
//     if (common.npot3c > 0)    
//         cpuBoTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              
// 
//     m += common.npot3c*inum;
//     n += common.npot3c*dim*inum;
//    k += common.npot3c*6*inum;
//     if (common.npot4a > 0)
//         cpuNonbondedQuadrupletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              
// 
//     m += common.npot4a*inum;
//     n += common.npot4a*dim*inum;
//    k+= common.npot4a*6*inum;
//     if (common.npot4b > 0)
//         cpuBondedQuadrupletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                                
}

#endif

