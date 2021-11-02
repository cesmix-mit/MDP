/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __GPUPOTENTIALS
#define __GPUPOTENTIALS

#include "gpuApp.cpp" 
       
void gpuSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *ai = &tmp.intmem[na]; // na
        Int *ti = &tmp.intmem[2*na]; // na
        dstype *xi = &tmp.tmpmem[0]; // na*dim
        dstype *qi = &tmp.tmpmem[na*dim]; // na*ncq
        gpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, na, ncq, dim);

        dstype *fi = &tmp.tmpmem[na*(dim+ncq)]; // na*dim
        dstype *ei = &tmp.tmpmem[na*(2*dim+ncq)]; // na
        dstype *du = &tmp.tmpmem[na*(2*dim+ncq+1)]; // na
        gpuComputeSingleEnergyForce(ei, du, fi, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        gpuSingleDecomposition(e, f, ei, fi, ai, na, dim);        
    }        
}
void gpuNonbondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        gpuSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void gpuBondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        gpuSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void gpuFullNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);                       
        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            gpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);    
        
#ifdef HAVE_DEBUG                      
        writearray2file("xijgpu.bin", xij, ntuples*dim, common.backend); 
        writearray2file("eijgpu.bin", eij, ntuples, common.backend); 
        writearray2file("fijgpu.bin", fij, ntuples*dim, common.backend); 
        writearray2file("egpu.bin", e, common.inum, common.backend); 
        writearray2file("fgpu.bin", f, common.inum*dim, common.backend); 
    //error("here");    
#endif                                        
    }        
}

void gpuHalfNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            gpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, ntuples, dim);
        else {
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);   
            gpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);            
            gpuNeighborAtomPairDecomposition(e, f, eij, fij, jlist, bnumsum, index, naj, dim);        
        }                        
    }        
}

void gpuNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            gpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            gpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}

void gpuBondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    //error("here");
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                gpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                gpuFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                gpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                gpuHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void gpuBO2EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);

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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *du = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq+1)]; // ntuples        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq+2)]; // ntuples*dim                   
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        dstype *gi = &tmp.tmpmem[2*na]; // na
        dstype *du2 = &tmp.tmpmem[3*na]; // na
        gpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        gpuPairDensityGradient(ei, du2, gi, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);
        gpuEmbedingForce(fij, gi, pairnum, pairnumsum, na);        
                
        gpuPutArrayAtIndex(e, ei, ilist, na);        
        if (decomp==0)
            gpuHalfForceDecomposition(f, fij, ai, aj, dim, ntuples);        
        else {                
            gpuIAtomDecomposition(f, fij, ilist, pairnumsum, dim, na);      
            gpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
            gpuJAtomDecomposition(f, fij, jlist, bnumsum, index, dim, naj);                
        }
        
//         if (dim==2) {
//             if (decomp==0)
//                 gpuHalfForceDecomposition2D(f, fij, ai, aj, ntuples);        
//             else {                
//                 gpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
//                 gpuArrayCopy(tmp.intmem, aj, ntuples);
//                 Int *jlist = &tmp.intmem[ntuples];   
//                 Int *bnumsum = &tmp.intmem[2*ntuples]; 
//                 Int *index = &tmp.intmem[3*ntuples]; // ntuples       
//                 Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
//                 Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
//                 Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
//                 Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
//                 Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);                
//             }
//         }
//         else {
//             if (decomp==0)
//                 gpuHalfForceDecomposition3D(f, fij, ai, aj, ntuples);    
//             else {
//                 gpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
//                 gpuArrayCopy(tmp.intmem, aj, ntuples);
//                 Int *jlist = &tmp.intmem[ntuples];   
//                 Int *bnumsum = &tmp.intmem[2*ntuples]; 
//                 Int *index = &tmp.intmem[3*ntuples]; // ntuples       
//                 Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
//                 Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
//                 Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
//                 Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
//                 Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                                   
//             }
//         }                                 
    }        
}
void gpuBoPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        gpuBO2EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void gpuTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax        
        if ((typej>0) && (typek>0))
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        gpuTripletnum(tripletnum, pairnum, na);
        //for (int ii=0; ii<na; ii++)
        //    tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       

        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        //Cumsum(tripletnumsum, tripletnum, na+1, backend);                                         
        //int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, backend);                             
        
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        gpuNeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        gpuArrayCopy(tripletlist, temp, 2*ntuples);
        
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
        gpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, ncq, dim);                     
              
        dstype *fij = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *eijk = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(4*dim+3*ncq+1)]; // ntuples
        
        gpuComputeTripletEnergyForce(eijk, du, fij, fik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                         param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);             
        if (decomp==0)            
            gpuTripletDecomposition(e, f, eijk, fij, fik, ai, aj, ak, ntuples, dim);
        else {
            gpuCenterAtomTripletDecomposition(e, f, eijk, fij, fik, ilist, tripletnumsum, na, dim);                
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            gpuNeighborAtomTripletDecomposition(e, f, eijk, fij, jlist, bnumsum, index, naj, dim);
            Int nak = gpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            gpuNeighborAtomTripletDecomposition(e, f, eijk, fik, jlist, bnumsum, index, nak, dim);
        }        
    }        
}
void gpuNonbondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        gpuTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void gpuBondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        gpuTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void gpuBO3EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax       
        if (typej == typei) {
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        
        }
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        

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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                                
        //gpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, npairs, potnum);
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs*neighmax        
        gpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
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
        
        gpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, neighmax, npairs, ncq, dim);                      
        
        dstype *f3ij = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *f3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *e3ijk = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du2 = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq+1)]; // ntuples
        gpuComputeTripletEnergyForce(e3ijk, du2, f3ij, f3ik, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *h3ij = &tmp.tmpmem[npairs*(1+dim)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2+dim)]; // npairs
        dstype *d3ij = &tmp.tmpmem[npairs*(3+dim)]; // npairs
        dstype *g3ij = &tmp.tmpmem[npairs*(4+dim)]; // dim*npairs
        
        gpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        gpuElectronDensity(g3ij, f3ij, tripletnum, tripletnumsum, npairs);
        gpuTripletDensityGradient(c3ij, du, d3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);
        gpuArrayABPXYZ(fij, c3ij, eij, d3ij, g3ij, dim, npairs);
        gpuArrayAXY(eij, eij, c3ij, 1.0, npairs);
                            
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
        gpuArrayCopy(&tmp.intmem[0], &tmp.intmem[1+3*na+na*neighmax], npairs);
        gpuArrayCopy(&tmp.intmem[npairs], &tmp.intmem[1+3*na+na*neighmax+npairs], npairs);
        gpuArrayCopy(&tmp.intmem[2*npairs], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax], ntuples);
        gpuArrayCopy(&tmp.intmem[2*npairs+ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+ntuples], ntuples);
        gpuArrayCopy(&tmp.intmem[2*npairs+2*ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*neighmax+2*ntuples], ntuples);

        if (decomp==0)
            gpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, npairs, dim);
        else {
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);   
            Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //npairs  
            Int *bnumsum = &tmp.intmem[3*npairs+3*ntuples]; //npairs  
            Int *index = &tmp.intmem[4*npairs+3*ntuples]; // npairs
            Int *p0 = &tmp.intmem[5*npairs+3*ntuples]; // npairs       
            Int *p1 = &tmp.intmem[6*npairs+3*ntuples]; // npairs       
            Int *p2 = &tmp.intmem[7*npairs+3*ntuples]; // npairs       
            Int *p3 = &tmp.intmem[8*npairs+3*ntuples]; // npairs                   
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
            gpuNeighborAtomPairDecomposition(e, f, eij, fij, jlist, bnumsum, index, naj, dim);                 
        }
        
        gpuTripletForceDecomposition(f3ik, eij, d3ij, tripletnum, tripletnumsum, npairs, dim);                

        if (decomp==0)
            gpuHalfForceDecomposition(f, f3ik, a3i, a3k, dim, ntuples);        
        else {
            Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
            Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
            Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
            Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
            Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
            Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
            Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
            Int nai = gpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
            gpuJAtomDecomposition(f, f3ik, jlist, bnumsum, index, dim, nai);                            
            Int nak = gpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
            gpuJAtomDecomposition(f, f3ik, jlist, bnumsum, index, dim, nak);          
        }        
//         if (dim==2) {            
//             if (decomp==0)
//                 gpuHalfForceDecomposition2D(f, f3ik, a3i, a3k, ntuples);        
//             else {
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
//                 Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
//                 Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
//                 Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
//                 Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
//                 Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
//                 Int nai = gpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nai);                            
//                 Int nak = gpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nak);          
//             }
//         }
//         else {
//             if (decomp==0)
//                 gpuHalfForceDecomposition3D(f, f3ik, a3i, a3k, ntuples);    
//             else {
//                 Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
//                 Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
//                 Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
//                 Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
//                 Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
//                 Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
//                 Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
//                 Int nai = gpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nai);                            
//                 Int nak = gpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
//                 gpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nak);         
//             }
//         }        
    }        
}

void gpuBoTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        gpuBO3EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void gpuQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax    
        if ((typej>0) && (typek>0) && (typel>0))
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, typel, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        //for (int ii=0; ii<na; ii++)
        //    quadrupletnum[ii] = (pairnum[ii]-2)*(pairnum[ii]-1)*pairnum[ii]/6;                                       
        gpuQuadrupletnum(quadrupletnum, pairnum, na);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        //Cumsum(quadrupletnumsum, quadrupletnum, na+1, backend);       
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                        
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        gpuNeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        gpuArrayCopy(quadrupletlist, temp, 3*ntuples);
        
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
        gpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, ncq, dim);       
        
        dstype *fij = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(5*dim+3*ncq)]; // ntuples*dim
        dstype *fil = &tmp.tmpmem[ntuples*(6*dim+3*ncq)]; // ntuples*dim
        dstype *eijkl = &tmp.tmpmem[ntuples*(7*dim+3*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(7*dim+3*ncq+1)]; // ntuples
        
        gpuComputeQuadrupletEnergyForce(eijkl, du, fij, fik, fil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                                

        if (decomp==0)
            gpuQuadrupletDecomposition(e, f, eijkl, fij, fik, fil, ai, aj, ak, al, ntuples, dim);
        else {
            gpuCenterAtomQuadrupletDecomposition(e, f, eijkl, fij, fik, fil, ilist, quadrupletnumsum, na, dim);   
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fij, jlist, bnumsum, index, naj, dim);
            Int nak = gpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fik, jlist, bnumsum, index, nak, dim);                
            Int nal = gpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, f, eijkl, fil, jlist, bnumsum, index, nal, dim);
        }         
    }        
}
void gpuNonbondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        gpuQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void gpuBondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        gpuQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void gpuEmpiricalPotentialEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    //error("here");
    if (common.npot1a > 0)
        gpuNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    if (common.npot1b > 0)
        gpuBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    if (common.npot2a > 0)
        gpuNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    if (common.npot2b > 0)
        gpuBondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    if (common.npot2c > 0)    
        gpuBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    if (common.npot3a > 0)    
        gpuNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    if (common.npot3b > 0)
        gpuBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    if (common.npot3c > 0)    
        gpuBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    if (common.npot4a > 0)
        gpuNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    if (common.npot4b > 0)
        gpuBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);              
    
    ArrayMinus(f, f, common.dim*common.inum, common.backend);
}

void gpuEmpiricalPotentialDescriptors(dstype *ei, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;
    
    int m = 0;
    int n = 0;
    if (common.npot1a > 0)
        gpuNonbondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    n += common.npot1a*dim*inum;
    if (common.npot1b > 0)
        gpuBondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    n += common.npot1b*dim*inum;
    if (common.npot2a > 0)
        gpuNonbondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    n += common.npot2a*dim*inum;
    if (common.npot2b > 0)
        gpuBondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    n += common.npot2b*dim*inum;
    if (common.npot2c > 0)    
        gpuBoPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    m += common.npot2c*inum;
    n += common.npot2c*dim*inum;
    if (common.npot3a > 0)    
        gpuNonbondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    n += common.npot3a*dim*inum;
    if (common.npot3b > 0)
        gpuBondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    n += common.npot3b*dim*inum;
    if (common.npot3c > 0)    
        gpuBoTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    n += common.npot3c*dim*inum;
    if (common.npot4a > 0)
        gpuNonbondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    n += common.npot4a*dim*inum;
    if (common.npot4b > 0)
        gpuBondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                                
}

void gpuSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *ai = &tmp.intmem[na]; // na
        Int *ti = &tmp.intmem[2*na]; // na
        dstype *xi = &tmp.tmpmem[0]; // na*dim
        dstype *qi = &tmp.tmpmem[na*dim]; // na*ncq
        gpuNeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, na, ncq, dim);
        
        dstype *ei = &tmp.tmpmem[na*(dim+ncq)]; // na
        gpuSingle(ei, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        gpuSingleDecomposition(e, ei, ai, na, dim);        
    }        
}
void gpuNonbondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        gpuSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void gpuBondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        gpuSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void gpuFullNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim        
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                
        if (decomp==0)
            gpuFullNeighPairDecomposition(e, eij, ai, ntuples, dim);
        else
            gpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);                        
    }        
}

void gpuHalfNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            gpuHalfNeighPairDecomposition(e, eij, ai, aj, ntuples, dim);
        else {
            gpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);   
            gpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);            
            gpuNeighborAtomPairDecomposition(e, eij, jlist, bnumsum, index, naj, dim);        
        }                
    }        
}
void gpuNonbondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            gpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            gpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void gpuBondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                gpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                gpuFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                gpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                gpuHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void gpuBO2Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);

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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples    
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        gpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        gpuPaircDensity(ei, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);                

        gpuPutArrayAtIndex(e, ei, ilist, na);        
    }        
}
void gpuBoPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        gpuBO2Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void gpuTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax        
        if ((typej>0) && (typek>0))
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        gpuTripletnum(tripletnum, pairnum, na);
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                                              
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
                
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        gpuNeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        gpuArrayCopy(tripletlist, temp, 2*ntuples);
        
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
        gpuNeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, ncq, dim);                     
              
        dstype *eijk = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples
        
        gpuTriplet(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak,
                         param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);             
        if (decomp==0)            
            gpuTripletDecomposition(e, eijk, ai, aj, ak, ntuples, dim);
        else {
            gpuCenterAtomTripletDecomposition(e,  eijk, ilist, tripletnumsum, na, dim);                
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 5*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 6*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            gpuNeighborAtomTripletDecomposition(e, eijk, jlist, bnumsum, index, naj, dim);
            Int nak = gpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            gpuNeighborAtomTripletDecomposition(e, eijk, jlist, bnumsum, index, nak, dim);
        }    
    }        
}
void gpuNonbondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        gpuTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void gpuBondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        gpuTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void gpuBO3Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax       
        if (typej == typei) {
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        
        }
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);        

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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       

         
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);        
                
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*neighmax]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*neighmax]; // npairs*neighmax        
        gpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
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
        
        gpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, neighmax, npairs, ncq, dim);          
        dstype *e3ijk = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+3*ncq)]; // ntuples
        gpuTriplet(e3ijk, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                
        dstype *h3ij = &tmp.tmpmem[npairs*(1)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2)]; // npairs
        gpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        gpuTripletcDensity(c3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);
            
        gpuArrayAXY(eij, eij, c3ij, 1.0, npairs);                        
        ArrayCopy(&tmp.intmem[0], &tmp.intmem[1+3*na+na*neighmax], npairs, backend);
        ArrayCopy(&tmp.intmem[npairs], &tmp.intmem[1+3*na+na*neighmax+npairs], npairs, backend);
        
        // pairnum, pairnumsum, pairlist, tripletnum, tripletlist, tripletnumsum, t3i, t3j, t3k -> ai, aj, a3i, a3j, a3k        
                
        if (decomp==0)
            gpuHalfNeighPairDecomposition(e, eij, ai, aj, npairs, dim);
        else {
            gpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);   
            Int *jlist = &tmp.intmem[2*npairs]; //npairs  
            Int *bnumsum = &tmp.intmem[3*npairs]; //npairs  
            Int *index = &tmp.intmem[4*npairs]; // npairs
            Int *p0 = &tmp.intmem[5*npairs]; // npairs       
            Int *p1 = &tmp.intmem[6*npairs]; // npairs       
            Int *p2 = &tmp.intmem[7*npairs]; // npairs       
            Int *p3 = &tmp.intmem[8*npairs]; // npairs                   
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, aj, p1, p2, p3, npairs);
            gpuNeighborAtomPairDecomposition(e, eij, jlist, bnumsum, index, naj, dim);                 
        }                 
    }        
}
void gpuBoTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        gpuBO3Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void gpuQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*neighmax    
        if ((typej>0) && (typek>0) && (typel>0))
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, typek, typel, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        gpuQuadrupletnum(quadrupletnum, pairnum, na);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                                
        Int *temp = &tmp.intmem[1 + 4*na + na*neighmax]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        gpuNeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, neighmax);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        gpuArrayCopy(quadrupletlist, temp, 3*ntuples);
        
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
        gpuNeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, ncq, dim);       
        
        dstype *eijkl = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
    
        gpuQuadruplet(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                                
        if (decomp==0)
            gpuQuadrupletDecomposition(e, eijkl,  ai, aj, ak, al, ntuples, dim);
        else {
            gpuCenterAtomQuadrupletDecomposition(e, eijkl, ilist, quadrupletnumsum, na, dim);   
            Int *jlist = &tmp.intmem[0];  // ntuples       
            Int *bnumsum = &tmp.intmem[ntuples]; // ntuples       
            Int *index = &tmp.intmem[2*ntuples]; // ntuples       
            Int *t0 = &tmp.intmem[1 + 3*na + 7*ntuples]; // ntuples       
            Int *t1 = &tmp.intmem[1 + 3*na + 8*ntuples]; // ntuples       
            Int *t2 = &tmp.intmem[1 + 3*na + 9*ntuples]; // ntuples       
            Int *t3 = &tmp.intmem[1 + 3*na + 10*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, t0, aj, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, naj, dim);
            Int nak = gpuUniqueSort(jlist, bnumsum, index, t0, ak, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, nak, dim);                
            Int nal = gpuUniqueSort(jlist, bnumsum, index, t0, al, t1, t2, t3, ntuples);
            gpuNeighborAtomQuadrupletDecomposition(e, eijkl, jlist, bnumsum, index, nal, dim);
        }        
    }        
}
void gpuNonbondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        gpuQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void gpuBondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        gpuQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void gpuEmpiricalPotentialDescriptors(dstype *ei, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;        
            
    int m = 0;
    if (common.npot1a > 0)
        gpuNonbondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    if (common.npot1b > 0)
        gpuBondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    if (common.npot2a > 0)
        gpuNonbondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    if (common.npot2b > 0)
        gpuBondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    if (common.npot2c > 0)    
        gpuBoPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
 
    m += common.npot2c*inum;
    if (common.npot3a > 0)    
        gpuNonbondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    if (common.npot3b > 0)
        gpuBondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    if (common.npot3c > 0)    
        gpuBoTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    if (common.npot4a > 0)
        gpuNonbondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    if (common.npot4b > 0)
        gpuBondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                            
}

void gpuFullNeighPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);
        else
            gpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim);        
                        
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);                       
        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            gpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);    
        
        // virial tally
        gpuVirialPairTally(v, fij, xij, -0.5, ai, dim, common.inum, ntuples);                
        
#ifdef HAVE_DEBUG                      
        writearray2file("xijgpu.bin", xij, ntuples*dim, common.backend); 
        writearray2file("eijgpu.bin", eij, ntuples, common.backend); 
        writearray2file("fijgpu.bin", fij, ntuples*dim, common.backend); 
        writearray2file("egpu.bin", e, common.inum, common.backend); 
        writearray2file("fgpu.bin", f, common.inum*dim, common.backend); 
    //error("here");    
#endif                                        
    }        
}

void gpuHalfNeighPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
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
            gpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = gpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            gpuArrayFill(ilist, e1, na);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim);            
        else
            gpuHalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, dim);
                                
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
        gpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            gpuHalfNeighPairDecomposition(e, f, eij, fij, ai, aj, ntuples, dim);
        else {
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);   
            gpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);            
            gpuNeighborAtomPairDecomposition(e, f, eij, fij, jlist, bnumsum, index, naj, dim);        
        }                        
        // virial tally
        gpuVirialPairTally(v, fij, xij, -0.5, ai, aj, dim, common.inum, ntuples);                
    }        
}

void gpuNonbondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            gpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            gpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}

void gpuBondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    //error("here");
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                gpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                gpuFullNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                gpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                gpuHalfNeighPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void gpuEmpiricalPotentialEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    //error("here");
    if (common.npot1a > 0)
        gpuNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    if (common.npot1b > 0)
        gpuBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    if (common.npot2a > 0)
        gpuNonbondedPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    if (common.npot2b > 0)
        gpuBondedPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

//     if (common.npot2c > 0)    
//         gpuBoPairEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
// 
//     if (common.npot3a > 0)    
//         gpuNonbondedTripletEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
//     
//     if (common.npot3b > 0)
//         gpuBondedTripletEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              
// 
//     if (common.npot3c > 0)    
//         gpuBoTripletEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              
// 
//     if (common.npot4a > 0)
//         gpuNonbondedQuadrupletEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              
// 
//     if (common.npot4b > 0)
//         gpuBondedQuadrupletEnergyForceVirial(e, f, v, nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);              
    
    ArrayMinus(f, f, common.dim*common.inum, common.backend);
}

void gpuNonbondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, int dim, int inum, int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        gpuSingleEnergyForce(&e[i*inum], &f[i*dim*inum], nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void gpuBondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, int dim, int inum, int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        gpuSingleEnergyForce(&e[i*inum], &f[i*dim*inum], nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void gpuNonbondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, int dim, int inum, int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            gpuFullNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            gpuHalfNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}

void gpuBondedPairEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, int dim, int inum, int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuFullNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                gpuFullNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                gpuFullNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                gpuHalfNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                gpuHalfNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                gpuHalfNeighPairEnergyForceVirial(&e[i*inum], &f[i*dim*inum], &v[i*6*inum], nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void gpuEmpiricalPotentialDescriptors(dstype *ei, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;
    
    int m = 0;
    int n = 0;
    int k = 0;
    if (common.npot1a > 0)
        gpuNonbondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[0]], dim, inum, nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    n += common.npot1a*dim*inum;
    if (common.npot1b > 0)
        gpuBondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[1]], dim, inum, nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    n += common.npot1b*dim*inum;    
    if (common.npot2a > 0)
        gpuNonbondedPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[2]], dim, inum, nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    n += common.npot2a*dim*inum;
    k += common.npot2a*6*inum;
    if (common.npot2b > 0)
        gpuBondedPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[3]], dim, inum, nparam[4]-nparam[3]);              

//     m += common.npot2b*inum;
//     n += common.npot2b*dim*inum;
//     k += common.npot2b*6*inum;
//     if (common.npot2c > 0)    
//         gpuBoPairEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
// 
//     m += common.npot2c*inum;
//     n += common.npot2c*dim*inum;
//     k += common.npot2c*6*inum;
//     if (common.npot3a > 0)    
//         gpuNonbondedTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
//     
//     m += common.npot3a*inum;
//     n += common.npot3a*dim*inum;
//     k += common.npot3a*6*inum;
//     if (common.npot3b > 0)
//         gpuBondedTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              
// 
//     m += common.npot3b*inum;
//     n += common.npot3b*dim*inum;
//     k += common.npot3b*6*inum;
//     if (common.npot3c > 0)    
//         gpuBoTripletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              
// 
//     m += common.npot3c*inum;
//     n += common.npot3c*dim*inum;
//     k += common.npot3c*6*inum;
//     if (common.npot4a > 0)
//         gpuNonbondedQuadrupletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              
// 
//     m += common.npot4a*inum;
//     n += common.npot4a*dim*inum;
//     k+= common.npot4a*6*inum;
//     if (common.npot4b > 0)
//         gpuBondedQuadrupletEnergyForceVirial(&ei[m], &f[n], &v[k], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
//     
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                                
}

#endif

