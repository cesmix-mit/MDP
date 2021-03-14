#ifndef __CALCULATION
#define __CALCULATION

#include "Calculation.h"
#include "Calculation_impl.cpp"

void CCalculation::NonbondedPairEnergyForce(dstype *e, dstype *f, dstype* x, dstype *q, Int *param, Int nparam) 
{    
    implNonbondedPairEnergyForce(e, f, nb, common, app, tmp, x, q, param, nparam);              
}

// void implNonbondedPairEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
//         dstype* x, dstype *q, dstype* param, Int nparam)


// void CCalculation::NeighPairs(dstype* x, dstype *q, Int *atomtype, Int istart, Int iend) 
// {
// 
// }
// 
// void CCalculation::NeighPairs(dstype* x, dstype *q, Int *atomtype, Int typei, Int istart, Int iend)
// {
// 
// }
// 
// void CCalculation::NeighPairs(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int istart, Int iend)
// {
// 
// }
// 
// void CCalculation::NeighTriplets(dstype* x, dstype *q, Int *atomtype, Int istart, Int iend)
// {
// 
// }
// 
// void NeighTriplets(dstype* x, dstype *q, Int *atomtype, Int typei, Int typej, Int typek, Int istart, Int iend)            
// {
// 
// }

#endif        
