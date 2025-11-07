#include <R.h> 
#include <Rinternals.h> 
#include <Rmath.h>
#include <R_ext/Rdynload.h>

SEXP monotone(SEXP lfdr){ 
    int i,n;
    double *vec, *out;
    vec = REAL(lfdr);
    n = length(lfdr);
    
    SEXP Rlfdr; 
    PROTECT(Rlfdr = allocVector(REALSXP,n));
    out = REAL(Rlfdr);
    out[0] = vec[0];
    for (i = 1; i < n; i++){
        if(vec[i] < out[(i-1)]){
            out[i] = out[(i-1)];
        }else{
            out[i] = vec[i];
        }
    }  
    UNPROTECT(1);
    return(Rlfdr);
}
static const
R_CallMethodDef callMethods[]  = {
    {"monotone", (DL_FUNC) &monotone, 1},
    {NULL, NULL, 0}
};

void R_init_ComBatSuite(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}

void R_unload_ComBatSuite(DllInfo *info)
{
    (void) info;
}
