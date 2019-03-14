#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BNPmix_cDDP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BNPmix_cDDP",    (DL_FUNC) &_BNPmix_cDDP,    18},
    {"_BNPmix_cICS",    (DL_FUNC) &_BNPmix_cICS,    15},
    {"_BNPmix_cICS_mv", (DL_FUNC) &_BNPmix_cICS_mv, 16},
    {"_BNPmix_cSLI",    (DL_FUNC) &_BNPmix_cSLI,    14},
    {"_BNPmix_cSLI_mv", (DL_FUNC) &_BNPmix_cSLI_mv, 16},
    {"_BNPmix_MAR",     (DL_FUNC) &_BNPmix_MAR,     14},
    {"_BNPmix_MAR_mv",  (DL_FUNC) &_BNPmix_MAR_mv,  16},
    {NULL, NULL, 0}
};

void R_init_BNPmix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
