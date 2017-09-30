#include "mex.h" /* Always include this */
#include <math.h>
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
#define ccr(ai,bi,ci) (ccr[ai+Nss[0]*(bi+Nss[1]*ci)])
#define cci(ai,bi,ci) (cci[ai+Nss[0]*(bi+Nss[1]*ci)])
#define kk1(ai,bi,ci) (kk1[ai+Nss[0]*(bi+Nss[1]*ci)])
#define kk2(ai,bi,ci) (kk2[ai+Nss[0]*(bi+Nss[1]*ci)])
#define kb(loc1,loc2,ai,bi) kb[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
#define avgdx(loc1,loc2,ai,bi) avgdx[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
#define avgdy(loc1,loc2,ai,bi) avgdy[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
    
    size_t ai, bi, ci;
    int NB1, NB2, loc1, loc2, di = 0;
    double *kk1, *kk2, *ccr, *cci, *kb, *avgdx, *avgdy;
    double EXT, aglMod, da, dr, r, agl, R_low;
    double temp_energy;
    ccr = mxGetPr(prhs[0]);
    cci = mxGetPi(prhs[0]);
    kk1 = mxGetPr(prhs[1]);
    kk2 = mxGetPr(prhs[2]);
    EXT = mxGetScalar(prhs[3]);
    aglMod = mxGetScalar(prhs[4]);
    da = mxGetScalar(prhs[5]);
    dr = mxGetScalar(prhs[6]);
    NB1 = mxGetScalar(prhs[7]);
    NB2 = mxGetScalar(prhs[8]);
    R_low = mxGetScalar(prhs[9]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    nrhs = 10;
    
    nlhs = 3;
    int ndim = 4, dims[4] = {NB1,NB2,Nss[0],Nss[1]};
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    kb = mxGetPr(plhs[0]);
    
    for (ai=0;ai<Nss[0];ai++) {
        for (bi=0;bi<Nss[1];bi++) {
            for (ci=0;ci<Nss[2];ci++) {
                if (kk1(ai,bi,ci)<EXT) {
                    r = sqrt(kk1(ai,bi,ci)*kk1(ai,bi,ci)+kk2(ai,bi,ci)*kk2(ai,bi,ci));
                    if (kk1(ai,bi,ci)>=0) {
                        agl = fmod(acos(kk2(ai,bi,ci)/r),aglMod);
                    }
                    else
                        agl = fmod(3.1415926-acos(kk2(ai,bi,ci)/r),aglMod);
                    loc1 = round((r-R_low)/dr)-1;  /* round(r)-1;*/
                    if (loc1<=-1) {
                        loc1 = loc1 + 1;
                    }
                    loc2 = ceil(agl/da)-1;
                    if (loc2<=-1) {
                        loc2 = loc2 + 1;
                    }
                    temp_energy = ccr(ai,bi,ci)*ccr(ai,bi,ci) + cci(ai,bi,ci)*cci(ai,bi,ci);
                    kb(loc1,loc2,ai,bi) = kb(loc1,loc2,ai,bi) + temp_energy;
                }
            }
        }
    }
    
    return;
}
