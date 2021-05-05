#ifndef __POTENTIALS
#define __POTENTIALS

#include "opuApp.cpp" 

#ifdef HAVE_CUDA    
#include "gpuApp.cpp" 
#endif
        
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
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *ai = &tmp.intmem[na]; // na
        Int *ti = &tmp.intmem[2*na]; // na
        dstype *xi = &tmp.tmpmem[0]; // na*dim
        dstype *qi = &tmp.tmpmem[na*dim]; // na*ncq
        NeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, na, ncq, dim, backend);

        dstype *fi = &tmp.tmpmem[na*(dim+ncq)]; // na*dim
        dstype *ei = &tmp.tmpmem[na*(2*dim+ncq)]; // na
        dstype *du = &tmp.tmpmem[na*(2*dim+ncq+1)]; // na
#ifdef HAVE_CUDA    
        gpuComputeSingleEnergyForce(ei, du, fi, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);
        gpuSingleDecomposition(e, f, ei, fi, ai, na, dim); 
#else        
        cpuComputeSingleEnergyForce(ei, du, fi, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        cpuSingleDecomposition(e, f, ei, fi, ai, na, dim);        
#endif        
    }        
}
void implNonbondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        implSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void implBondedSingleEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        implSingleEnergyForce(e, f, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void implFullNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
                
//         printArray2D(pairnum, 1, na, common.backend);  
//         printArray2D(nb.neighnum, 1, na, common.backend);  
//         printArray2D(pairlist, jnum, na, common.backend);  
//         printArray2D(nb.neighlist, jnum, na, common.backend);          
        
        //a list contains the starting positions of the first neighbor              
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
        
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
#ifdef HAVE_CUDA            
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
        if (decomp==0)
            gpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else
            gpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);                        
#else                        
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        if (decomp==0)
            cpuFullNeighPairDecomposition(e, f, eij, fij, ai, ntuples, dim);
        else
            cpuCenterAtomPairDecomposition(e, f, eij, fij, ilist, pairnumsum, na, dim);                        
#endif        
    }        
}
void implHalfNeighPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);            
        else
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);
                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
        dstype *eij = &tmp.tmpmem[ntuples*(2*dim+2*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(2*dim+2*ncq+1)]; // ntuples
#ifdef HAVE_CUDA            
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
#else        
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
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
#endif        
    }        
}


void implNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void implBondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                implFullNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                implHalfNeighPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void implBO2EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                                 
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                        
        dstype *du = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq+1)]; // ntuples        
        dstype *fij = &tmp.tmpmem[ntuples*(dim+2*ncq+2)]; // ntuples*dim                   
#ifdef HAVE_CUDA    
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        dstype *gi = &tmp.tmpmem[2*na]; // na
        dstype *du2 = &tmp.tmpmem[3*na]; // na
        gpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        gpuPaircDensity(ei, du2, gi, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);
        gpuEmbedingForce(fij, gi, pairnum, pairnumsum, na);        
                
        gpuPutArrayAtIndex(e, ei, ilist, na);        
        if (dim==2) {
            if (decomp==0)
                gpuHalfForceDecomposition2D(f, fij, ai, aj, ntuples);        
            else {                
                gpuIAtomDecomposition2D(f, fij, ilist, pairnumsum, na);      
                gpuArrayCopy(tmp.intmem, aj, ntuples);
                Int *jlist = &tmp.intmem[ntuples];   
                Int *bnumsum = &tmp.intmem[2*ntuples]; 
                Int *index = &tmp.intmem[3*ntuples]; // ntuples       
                Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
                Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
                Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
                Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
                Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
                gpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, naj);                
            }
        }
        else {
            if (decomp==0)
                gpuHalfForceDecomposition3D(f, fij, ai, aj, ntuples);    
            else {
                gpuIAtomDecomposition3D(f, fij, ilist, pairnumsum, na);     
                gpuArrayCopy(tmp.intmem, aj, ntuples);
                Int *jlist = &tmp.intmem[ntuples];   
                Int *bnumsum = &tmp.intmem[2*ntuples]; 
                Int *index = &tmp.intmem[3*ntuples]; // ntuples       
                Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
                Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
                Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
                Int *p3 = &tmp.intmem[7*ntuples]; // ntuples                       
                Int naj = gpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);
                gpuJAtomDecomposition3D(f, fij, jlist, bnumsum, index, naj);                                   
            }
        }                        
#else        
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        dstype *gi = &tmp.tmpmem[2*na]; // na
        dstype *du2 = &tmp.tmpmem[3*na]; // na
        cpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        cpuPaircDensity(ei, du2, gi, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);
        cpuEmbedingForce(fij, gi, pairnum, pairnumsum, na);        
                
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
#endif        
    }        
}
void implBoPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        implBO2EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void implTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum        
        if ((typej>0) && (typek>0))
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        ArrayTripletnum(tripletnum, pairnum, na, backend);
        //for (int ii=0; ii<na; ii++)
        //    tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       

        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        //cpuCumsum(tripletnumsum, tripletnum, na+1);                                         
        //int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, backend);                             
        
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        NeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum, backend);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        ArrayCopy(tripletlist, temp, 2*ntuples, backend);
        
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
        NeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, ncq, dim, backend);                     
              
        dstype *fij = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *eijk = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(4*dim+3*ncq+1)]; // ntuples
        
#ifdef HAVE_CUDA      
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
#else        
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
#endif    
}
void implNonbondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        implTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void implBondedTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        implTripletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void implBO3EnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum       
        if (typej == typei) {
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);        
        }
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1    
        //cpuCumsum(pairnumsum, pairnum, na+1);                                         
        //int npairs = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int npairs = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // npairs    
        Int *aj = &tmp.intmem[1+3*na+npairs+na*jnum]; // npairs        
        Int *ti = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs        
        Int *tj = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs        
        dstype *eij = &tmp.tmpmem[0]; // npairs
        dstype *fij = &tmp.tmpmem[npairs]; // npairs*dim
        dstype *xij = &tmp.tmpmem[npairs*(1+dim)]; // npairs*dim
        dstype *qi = &tmp.tmpmem[npairs*(1+2*dim)]; // npairs*ncq
        dstype *qj = &tmp.tmpmem[npairs*(1+2*dim+ncq)]; // npairs*ncq
        dstype *du = &tmp.tmpmem[npairs*(1+2*dim+2*ncq)]; // npairs
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                                
        //cpuComputePairEnergyForce(eij, fij, xij, qi, qj, ti, tj, ai, aj, param, dim, ncq, nparam, npairs, potnum);
#ifdef HAVE_CUDA      
        gpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);        
#else
        cpuComputePairEnergyForce(eij, du, fij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);
#endif        
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs*jnum        
#ifdef HAVE_CUDA              
        gpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, typek, dim);                
#else        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, typek, dim);                
#endif        
        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[1+3*na+3*npairs+(na+npairs)*jnum]; // npairs+1        
        //cpuCumsum(tripletnumsum, tripletnum, npairs+1);                                 
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum], &tmp.intmem[3+3*na+5*npairs+(na+npairs)*jnum], npairs+1, backend);                                         
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
#ifdef HAVE_CUDA                      
        gpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, jnum, npairs, ncq, dim);                              
#else        
        cpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, jnum, npairs, ncq, dim);                      
#endif
        
        dstype *f3ij = &tmp.tmpmem[npairs*(1+dim)+ntuples*(2*dim+3*ncq)]; // ntuples*dim
        dstype *f3ik = &tmp.tmpmem[npairs*(1+dim)+ntuples*(3*dim+3*ncq)]; // ntuples*dim
        dstype *e3ijk = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq)]; // ntuples
        dstype *du2 = &tmp.tmpmem[npairs*(1+dim)+ntuples*(4*dim+3*ncq+1)]; // ntuples
#ifdef HAVE_CUDA              
        gpuComputeTripletEnergyForce(e3ijk, du2, f3ij, f3ik, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
#else        
        cpuComputeTripletEnergyForce(e3ijk, du2, f3ij, f3ik, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
#endif  
        
        dstype *h3ij = &tmp.tmpmem[npairs*(1+dim)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2+dim)]; // npairs
        dstype *d3ij = &tmp.tmpmem[npairs*(3+dim)]; // npairs
        dstype *g3ij = &tmp.tmpmem[npairs*(4+dim)]; // dim*npairs
#ifdef HAVE_CUDA             
        gpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        gpuElectronDensity(g3ij, f3ij, tripletnum, tripletnumsum, npairs);
        gpuTripletcDensity(c3ij, du, d3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);        
        gpuArrayABPXYZ(fij, c3ij, eij, d3ij, g3ij, dim, npairs);
        gpuArrayAXY(eij, eij, c3ij, 1.0, npairs);
#else                
        cpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        cpuElectronDensity(g3ij, f3ij, tripletnum, tripletnumsum, npairs);
        cpuTripletcDensity(c3ij, du, d3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);
        for (int i=0; i<npairs; i++) 
            for (int j=0; j<dim; j++)
                fij[j+3*i] = fij[j+3*i]*c3ij[i] + eij[i]*d3ij[i]*g3ij[j+3*i];   
        
        for (int i=0; i<npairs; i++) 
            eij[i] = eij[i]*c3ij[i];        
#endif                  
                    
        // pairnum, pairnumsum, pairlist, tripletnum, tripletlist, tripletnumsum, t3i, t3j, t3k -> ai, aj, a3i, a3j, a3k
//         for (int i=0; i<npairs; i++) {
//             tmp.intmem[i] = tmp.intmem[1+3*na+na*jnum+i]; // ai
//             tmp.intmem[npairs+i] = tmp.intmem[1+3*na+npairs+na*jnum+i]; // aj
//         }
//         for (int i=0; i<ntuples; i++) {
//             tmp.intmem[2*npairs+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+i]; //a3i
//             tmp.intmem[2*npairs+ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+ntuples+i]; //a3j
//             tmp.intmem[2*npairs+2*ntuples+i] = tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+2*ntuples+i]; //a3k
//         }
        ArrayCopy(&tmp.intmem[0], &tmp.intmem[1+3*na+na*jnum], npairs, backend);
        ArrayCopy(&tmp.intmem[npairs], &tmp.intmem[1+3*na+na*jnum+npairs], npairs, backend);
        ArrayCopy(&tmp.intmem[2*npairs], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum], ntuples, backend);
        ArrayCopy(&tmp.intmem[2*npairs+ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+ntuples], ntuples, backend);
        ArrayCopy(&tmp.intmem[2*npairs+2*ntuples], &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+2*ntuples], ntuples, backend);

#ifdef HAVE_CUDA           
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
        
        if (dim==2) {            
            if (decomp==0)
                gpuHalfForceDecomposition2D(f, f3ik, a3i, a3k, ntuples);        
            else {
                Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
                Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
                Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
                Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
                Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
                Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
                Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
                Int nai = gpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
                gpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nai);                            
                Int nak = gpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
                gpuJAtomDecomposition2D(f, f3ik, jlist, bnumsum, index, nak);          
            }
        }
        else {
            if (decomp==0)
                gpuHalfForceDecomposition3D(f, f3ik, a3i, a3k, ntuples);    
            else {
                Int *jlist = &tmp.intmem[2*npairs+3*ntuples]; //ntuples  
                Int *bnumsum = &tmp.intmem[2*npairs+4*ntuples]; //ntuples  
                Int *index = &tmp.intmem[2*npairs+5*ntuples]; // npairs
                Int *p0 = &tmp.intmem[2*npairs+6*ntuples]; //ntuples       
                Int *p1 = &tmp.intmem[2*npairs+7*ntuples]; //ntuples       
                Int *p2 = &tmp.intmem[2*npairs+8*ntuples]; //ntuples       
                Int *p3 = &tmp.intmem[2*npairs+9*ntuples]; //ntuples       
                Int nai = gpuUniqueSort(jlist, bnumsum, index, p0, a3i, p1, p2, p3, ntuples);
                gpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nai);                            
                Int nak = gpuUniqueSort(jlist, bnumsum, index, p0, a3k, p1, p2, p3, ntuples);
                gpuJAtomDecomposition3D(f, f3ik, jlist, bnumsum, index, nak);         
            }
        }                
#else        
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
#endif                
    }        
}

void implBoTripletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        implBO3EnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void implQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum    
        if ((typej>0) && (typek>0) && (typel>0))
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, typel, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        //for (int ii=0; ii<na; ii++)
        //    quadrupletnum[ii] = (pairnum[ii]-2)*(pairnum[ii]-1)*pairnum[ii]/6;                                       
        ArrayQuadrupletnum(quadrupletnum, pairnum, na, backend);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        //cpuCumsum(quadrupletnumsum, quadrupletnum, na+1);       
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                        
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        NeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum, backend);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        ArrayCopy(quadrupletlist, temp, 3*ntuples, backend);
        
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
        NeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, ncq, dim, backend);       
        
        dstype *fij = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples*dim
        dstype *fik = &tmp.tmpmem[ntuples*(5*dim+3*ncq)]; // ntuples*dim
        dstype *fil = &tmp.tmpmem[ntuples*(6*dim+3*ncq)]; // ntuples*dim
        dstype *eijkl = &tmp.tmpmem[ntuples*(7*dim+3*ncq)]; // ntuples
        dstype *du = &tmp.tmpmem[ntuples*(7*dim+3*ncq+1)]; // ntuples
        
#ifdef HAVE_CUDA           
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
#else        
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
#endif        
    }        
}
void implNonbondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        implQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void implBondedQuadrupletEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        implQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void implEmpiricalPotentialEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    if (common.npot1a > 0)
        implNonbondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    if (common.npot1b > 0)
        implBondedSingleEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    if (common.npot2a > 0)
        implNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    if (common.npot2b > 0)
        implBondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    if (common.npot2c > 0)    
        implBoPairEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    if (common.npot3a > 0)    
        implNonbondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    if (common.npot3b > 0)
        implBondedTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    if (common.npot3c > 0)    
        implBoTripletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    if (common.npot4a > 0)
        implNonbondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    if (common.npot4b > 0)
        implBondedQuadrupletEnergyForce(e, f, nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);              
    
    ArrayMinus(f, f, common.dim*common.inum, common.backend);
}

void implEmpiricalPotentialDescriptors(dstype *ei, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;
    
    int m = 0;
    int n = 0;
    if (common.npot1a > 0)
        implNonbondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    n += common.npot1a*dim*inum;
    if (common.npot1b > 0)
        implBondedSingleEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    n += common.npot1b*dim*inum;
    if (common.npot2a > 0)
        implNonbondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    n += common.npot2a*dim*inum;
    if (common.npot2b > 0)
        implBondedPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    n += common.npot2b*dim*inum;
    if (common.npot2c > 0)    
        implBoPairEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             

    m += common.npot2c*inum;
    n += common.npot2c*dim*inum;
    if (common.npot3a > 0)    
        implNonbondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    n += common.npot3a*dim*inum;
    if (common.npot3b > 0)
        implBondedTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    n += common.npot3b*dim*inum;
    if (common.npot3c > 0)    
        implBoTripletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    n += common.npot3c*dim*inum;
    if (common.npot4a > 0)
        implNonbondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    n += common.npot4a*dim*inum;
    if (common.npot4b > 0)
        implBondedQuadrupletEnergyForce(&ei[m], &f[n], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                                
}

void implSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int *atomtype, Int nparam, Int typei, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {        
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *ai = &tmp.intmem[na]; // na
        Int *ti = &tmp.intmem[2*na]; // na
        dstype *xi = &tmp.tmpmem[0]; // na*dim
        dstype *qi = &tmp.tmpmem[na*dim]; // na*ncq
        NeighSingles(xi, qi, x, q, ai, ti, ilist, atomtype, na, ncq, dim, backend);
        
        dstype *ei = &tmp.tmpmem[na*(dim+ncq)]; // na
#ifdef HAVE_CUDA          
        gpuSingle(ei, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        gpuSingleDecomposition(e, ei, ai, na, dim);                
#else                
        cpuSingle(ei, xi, qi, ti, ai, param, app.eta, app.kappa, 
                dim, ncq, nparam, common.neta, common.nkappa, na, potnum, common.bondtype);        
        cpuSingleDecomposition(e, ei, ai, na, dim);        
#endif                
    }        
}
void implNonbondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot1a; i++)
        implSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, 0, common.pot1a[i]);            
    
}
void implBondedSingleEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction    
    for (int i = 0; i < common.npot1b; i++)
        implSingleEnergy(e, nb, common, app, tmp, x, q, param, nb.atomtype, nparam, common.atom1b[i], common.pot1b[i]);                
}

void implFullNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
                        
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);        
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim        
#ifdef HAVE_CUDA            
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
                
        if (decomp==0)
            gpuFullNeighPairDecomposition(e, eij, ai, ntuples, dim);
        else
            gpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);                
#else                        
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);
                
        if (decomp==0)
            cpuFullNeighPairDecomposition(e, eij, ai, ntuples, dim);
        else
            cpuCenterAtomPairDecomposition(e, eij, ilist, pairnumsum, na, dim);                        
#endif                        
    }        
}

void implHalfNeighPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);            
        else
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);
                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);           
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples*dim
#ifdef HAVE_CUDA            
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
#else        
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
#endif        
    }        
}
void implNonbondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2a; i++)
            implFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);            
    }
    else {
        for (int i = 0; i < common.npot2a; i++)
            implHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2a[i], 
                nb.atomtype, nparam, 0, 0, common.decomposition, common.pot2a[i]);                    
    }        
}
void implBondedPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // bonded interaction
    if (common.neighpair == 0) {
        for (int i = 0; i < common.npot2b; i++) {
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                implFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
            }
            else {
                implFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);            
                implFullNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);            
            }
        }
    }
    else {
        for (int i = 0; i < common.npot2b; i++)
            if (common.atom2b[2*i] == common.atom2b[1+2*i]) {
                implHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
            }
            else {
                implHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[2*i], common.atom2b[1+2*i], common.decomposition, common.pot2b[i]);                    
                implHalfNeighPairEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2b[i], 
                    nb.atomtype, nparam, common.atom2b[1+2*i], common.atom2b[2*i], common.decomposition, common.pot2b[i]);                    
            }
    }            
}

void implBO2Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1                
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        dstype *qi = &tmp.tmpmem[ntuples*dim]; // ntuples*ncq
        dstype *qj = &tmp.tmpmem[ntuples*(dim+ncq)]; // ntuples*ncq
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       
                        
        dstype *eij = &tmp.tmpmem[ntuples*(dim+2*ncq)]; // ntuples    
#ifdef HAVE_CUDA            
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        gpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        gpuPaircDensity(ei, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);                
        gpuPutArrayAtIndex(e, ei, ilist, na);                
#else        
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);        
        dstype *rhoi = &tmp.tmpmem[0]; // na
        dstype *ei = &tmp.tmpmem[na]; // na
        cpuElectronDensity(rhoi, eij, pairnum, pairnumsum, na); 
        opuPaircDensity(ei, rhoi, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, na, potnum);                
        cpuPutArrayAtIndex(e, ei, ilist, na);        
#endif        
    }        
}
void implBoPairEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot2c; i++)
        implBO2Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq2c[i], 
                    nb.atomtype, nparam, common.atom2c[2*i], common.atom2c[1+2*i], common.decomposition, common.pot2c[i]);                    
}

void implTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum        
        if ((typej>0) && (typek>0))
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
                
        Int *tripletnum = &tmp.intmem[na]; // na        
        ArrayTripletnum(tripletnum, pairnum, na, backend);
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                                              
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, na, common.backend);     
                
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, tripletnum, tripletnumsum, pairnum, pairlist, temp)            
        NeighTripletList(temp, tripletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum, backend);                                
        Int *tripletlist = &tmp.intmem[1 + 3*na]; // 2*ntuples    (ilist, tripletnum, tripletnumsum, tripletlist)            
        ArrayCopy(tripletlist, temp, 2*ntuples, backend);
        
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
        NeighTriplets(xij, xik, qi, qj, qk, x, q, ai, aj, ak, ti, tj, tk, tripletnum, tripletlist, tripletnumsum, 
                ilist, nb.alist, atomtype, na, ncq, dim, backend);                     
              
        dstype *eijk = &tmp.tmpmem[ntuples*(2*dim+3*ncq)]; // ntuples
#ifdef HAVE_CUDA      
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
#else        
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
#endif        
    }        
}
void implNonbondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot3a; i++)
        implTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3a[i], 
            nb.atomtype, nparam, 0, 0, 0, common.decomposition, common.pot3a[i]);            

}
void implBondedTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot3b; i++)
        implTripletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3b[i], 
            nb.atomtype, nparam, common.atom3b[3*i], common.atom3b[1+3*i], common.atom3b[2+3*i], common.decomposition, common.pot3b[i]);                
}

void implBO3Energy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum       
        if (typej == typei) {
            HalfNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);        
        }
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim, backend);        

        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1            
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int npairs = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // npairs    
        Int *aj = &tmp.intmem[1+3*na+npairs+na*jnum]; // npairs        
        Int *ti = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs        
        Int *tj = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs        
        dstype *eij = &tmp.tmpmem[0]; // npairs
        dstype *xij = &tmp.tmpmem[npairs*(1)]; // npairs*dim
        dstype *qi = &tmp.tmpmem[npairs*(1+dim)]; // npairs*ncq
        dstype *qj = &tmp.tmpmem[npairs*(1+dim+ncq)]; // npairs*ncq
        dstype *du = &tmp.tmpmem[npairs*(1+dim+2*ncq)]; // npairs
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim, backend);       

#ifdef HAVE_CUDA              
        gpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);
#else        
        cpuPair(eij, xij, qi, qj, ti, tj, ai, aj, param, app.eta, app.kappa, dim, ncq, 
                nparam, common.neta, common.nkappa, npairs, potnum, 3);        
#endif        
                
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj)                 
        Int *tripletnum = &tmp.intmem[1+3*na+2*npairs+na*jnum]; // npairs
        Int *tripletlist = &tmp.intmem[1+3*na+3*npairs+na*jnum]; // npairs*jnum        
#ifdef HAVE_CUDA              
        gpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, typek, dim);                
#else        
        cpuNeighTripletList(tripletnum, tripletlist, x, rcutsq, pairnum, pairnumsum, pairlist, atomtype, ilist, nb.alist, nb.neighlist, 
                nb.neighnum, na, jnum, typek, dim);                
#endif        
                        
        //a list contains the starting positions of the first neighbor 
        Int *tripletnumsum = &tmp.intmem[1+3*na+3*npairs+(na+npairs)*jnum]; // npairs+1                 
        Cumsum(tripletnumsum, tripletnum, &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum], &tmp.intmem[3+3*na+5*npairs+(na+npairs)*jnum], npairs+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(tripletnumsum, npairs, common.backend);     
        
        // (ilist, pairnum, pairlist, pairnumsum, ai, aj, tripletnum, tripletnumsum, tripletlist)         
        Int *a3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum]; // ntuples        
        Int *a3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+ntuples]; // ntuples   
        Int *a3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+2*ntuples]; // ntuples   
        Int *t3i = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+3*ntuples]; // ntuples        
        Int *t3j = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+4*ntuples]; // ntuples        
        Int *t3k = &tmp.intmem[2+3*na+4*npairs+(na+npairs)*jnum+5*ntuples]; // ntuples        
        dstype *x3ij = &tmp.tmpmem[npairs*(1)]; // ntuples*dim
        dstype *x3ik = &tmp.tmpmem[npairs*(1)+ntuples*dim]; // ntuples*dim
        dstype *q3i = &tmp.tmpmem[npairs*(1)+2*ntuples*dim]; // ntuples*ncq
        dstype *q3j = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+ncq)]; // ntuples*ncq
        dstype *q3k = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+2*ncq)]; // ntuples*ncq
        
#ifdef HAVE_CUDA                      
        gpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, jnum, npairs, ncq, dim);               
        dstype *e3ijk = &tmp.tmpmem[npairs*(1)+ntuples*(2*dim+3*ncq)]; // ntuples
        gpuTriplet(e3ijk, x3ij, x3ik, q3i, q3j, q3k, t3i, t3j, t3k, a3i, a3j, a3k,
                       param, app.eta, app.kappa, dim, ncq, nparam, common.neta, common.nkappa, ntuples, potnum, common.bondtype);                
        dstype *h3ij = &tmp.tmpmem[npairs*(1)]; // npairs        
        dstype *c3ij = &tmp.tmpmem[npairs*(2)]; // npairs
        gpuElectronDensity(h3ij, e3ijk, tripletnum, tripletnumsum, npairs);
        gpuTripletcDensity(c3ij, h3ij, param, app.eta, app.kappa, 1, nparam, common.neta, common.nkappa, npairs, potnum);        
        gpuArrayAXY(eij, eij, c3ij, 1.0, npairs);        
        ArrayCopy(&tmp.intmem[0], &tmp.intmem[1+3*na+na*jnum], npairs, backend);
        ArrayCopy(&tmp.intmem[npairs], &tmp.intmem[1+3*na+na*jnum+npairs], npairs, backend);
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
#else        
        cpuNeighTriplets(x3ij, x3ik, q3i, q3j, q3k, x, q, a3i, a3j, a3k, t3i, t3j, t3k, tripletnum, tripletlist, tripletnumsum, 
                nb.alist, atomtype, jnum, npairs, ncq, dim);          
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
            tmp.intmem[i] = tmp.intmem[1+3*na+na*jnum+i]; // ai
            tmp.intmem[npairs+i] = tmp.intmem[1+3*na+npairs+na*jnum+i]; // aj
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
#endif                
    }        
}
void implBoTripletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 2; // non-bonded interaction
    for (int i = 0; i < common.npot3c; i++)
        implBO3Energy(e, nb, common, app, tmp, x, q, param, &app.rcutsq3c[i], nb.atomtype, nparam, 
                common.atom3c[3*i], common.atom3c[1+3*i], common.atom3c[2+3*i], common.decomposition, common.pot3c[i]);                    
}

void implQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, 
        dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int typek, Int typel, Int decomp, Int potnum)
{        
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int ncq = common.ncq;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                                
        Int *pairnum  = &tmp.intmem[1 + 3*na]; // na
        Int *pairlist = &tmp.intmem[1 + 4*na]; // na*jnum    
        if ((typej>0) && (typek>0) && (typel>0))
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, typek, typel, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim, backend);        
        
        Int *quadrupletnum = &tmp.intmem[na]; // na        
        ArrayQuadrupletnum(quadrupletnum, pairnum, na, backend);
        
        //a list contains the starting positions of the first neighbor 
        Int *quadrupletnumsum = &tmp.intmem[2*na]; // na+1        
        Cumsum(quadrupletnumsum, quadrupletnum, &tmp.intmem[3*na+na*jnum+1], &tmp.intmem[4*na+na*jnum+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(quadrupletnumsum, na, common.backend);     
                                
        Int *temp = &tmp.intmem[1 + 4*na + na*jnum]; //  (ilist, quadrupletnum, quadrupletnumsum, pairnum, pairlist, tmp)            
        NeighQuadrupletList(temp, quadrupletnumsum, pairnum, pairlist, ilist, nb.alist, na, jnum, backend);                                
        Int *quadrupletlist = &tmp.intmem[1 + 3*na]; // 3*ntuples    (ilist, quadrupletnum, quadrupletnumsum, quadrupletlist)            
        ArrayCopy(quadrupletlist, temp, 3*ntuples, backend);
        
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
        NeighQuadruplets(xij, xik, xil, qi, qj, qk, ql, x, q, ai, aj, ak, al, ti, tj, tk, tl,
                quadrupletnum, quadrupletlist, quadrupletnumsum, ilist, nb.alist, atomtype, na, ncq, dim, backend);       
        
        dstype *eijkl = &tmp.tmpmem[ntuples*(4*dim+3*ncq)]; // ntuples
#ifdef HAVE_CUDA           
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
#else        
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
#endif        
    }        
}
void implNonbondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 0; // non-bonded interaction    
    for (int i = 0; i < common.npot4a; i++)
        implQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4a[i], 
            nb.atomtype, nparam, 0, 0, 0, 0, common.decomposition, common.pot4a[i]);            

}
void implBondedQuadrupletEnergy(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, dstype *q, dstype* param, Int nparam)
{    
    common.bondtype = 1; // non-bonded interaction
    for (int i = 0; i < common.npot4b; i++)
        implQuadrupletEnergy(e, nb, common, app, tmp, x, q, param, &app.rcutsq4b[i], nb.atomtype, nparam, 
           common.atom4b[4*i], common.atom4b[1+4*i], common.atom4b[2+4*i], common.atom4b[3+4*i], common.decomposition, common.pot4b[i]);                
}

void implEmpiricalPotentialDescriptors(dstype *ei, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, dstype* x, dstype *q, dstype *param, Int *nparam) 
{    
    int inum = common.inum;
    int dim = common.dim;        
            
    int m = 0;
    if (common.npot1a > 0)
        implNonbondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[0]], nparam[1]-nparam[0]);              

    m += common.npot1a*inum;
    if (common.npot1b > 0)
        implBondedSingleEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[1]], nparam[2]-nparam[1]);              

    m += common.npot1b*inum;
    if (common.npot2a > 0)
        implNonbondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[2]], nparam[3]-nparam[2]);              

    m += common.npot2a*inum;
    if (common.npot2b > 0)
        implBondedPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[3]], nparam[4]-nparam[3]);              

    m += common.npot2b*inum;
    if (common.npot2c > 0)    
        implBoPairEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[4]], nparam[5]-nparam[4]);             
 
    m += common.npot2c*inum;
    if (common.npot3a > 0)    
        implNonbondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[5]], nparam[6]-nparam[5]);              
    
    m += common.npot3a*inum;
    if (common.npot3b > 0)
        implBondedTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[6]], nparam[7]-nparam[6]);              

    m += common.npot3b*inum;
    if (common.npot3c > 0)    
        implBoTripletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[7]], nparam[8]-nparam[7]);              

    m += common.npot3c*inum;
    if (common.npot4a > 0)
        implNonbondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[8]], nparam[9]-nparam[8]);              

    m += common.npot4a*inum;
    if (common.npot4b > 0)
        implBondedQuadrupletEnergy(&ei[m], nb, common, app, tmp, x, q, &param[nparam[9]], nparam[10]-nparam[9]);                  
    
//     dstype *onevec =  &tmp.tmpmem[0];  
//     ArraySetValue(onevec, 1.0, inum, common.backend);
//     PGEMTV(common.cublasHandle, inum, common.Nempot, &one, ei, inum, onevec, inc1, &one, e, inc1, common.backend);                            
}


#endif

