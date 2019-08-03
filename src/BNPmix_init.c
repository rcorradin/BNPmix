#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _BNPmix_cICS_mv_MKR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_mv_MKR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_mv_MKR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cDDP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS_mv_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cICS_mv_P(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_mv_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_cSLI_mv_P(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_mv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_mv_L(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_MAR_mv_P(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BNPmix_BNPmix_psm(SEXP);
extern SEXP _BNPmix_clean_partition(SEXP);
extern SEXP _BNPmix_BNPmix_VI_LB(SEXP, SEXP);
extern SEXP _BNPmix_BNPmix_BIN(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    // REG
    {"_BNPmix_cICS_mv_MKR",     (DL_FUNC) &_BNPmix_cICS_mv_MKR,     25},
    {"_BNPmix_cSLI_mv_MKR",     (DL_FUNC) &_BNPmix_cSLI_mv_MKR,     24},
    {"_BNPmix_MAR_mv_MKR",      (DL_FUNC) &_BNPmix_MAR_mv_MKR,      25},
    // DDP
    {"_BNPmix_cDDP",      (DL_FUNC) &_BNPmix_cDDP,      18},
    // ICS
    {"_BNPmix_cICS",      (DL_FUNC) &_BNPmix_cICS,      22},
    {"_BNPmix_cICS_L",    (DL_FUNC) &_BNPmix_cICS_L,    20},
    {"_BNPmix_cICS_mv",   (DL_FUNC) &_BNPmix_cICS_mv,   23},
    {"_BNPmix_cICS_mv_L", (DL_FUNC) &_BNPmix_cICS_mv_L, 21},
    {"_BNPmix_cICS_mv_P", (DL_FUNC) &_BNPmix_cICS_mv_P, 23},
    // SLI
    {"_BNPmix_cSLI",      (DL_FUNC) &_BNPmix_cSLI,      21},
    {"_BNPmix_cSLI_L",    (DL_FUNC) &_BNPmix_cSLI_L,    19},
    {"_BNPmix_cSLI_mv",   (DL_FUNC) &_BNPmix_cSLI_mv,   22},
    {"_BNPmix_cSLI_mv_L", (DL_FUNC) &_BNPmix_cSLI_mv_L, 20},
    {"_BNPmix_cSLI_mv_P", (DL_FUNC) &_BNPmix_cSLI_mv_P, 22},
    // MAR
    {"_BNPmix_MAR",       (DL_FUNC) &_BNPmix_MAR,       21},
    {"_BNPmix_MAR_L",     (DL_FUNC) &_BNPmix_MAR_L,     19},
    {"_BNPmix_MAR_mv",    (DL_FUNC) &_BNPmix_MAR_mv,    22},
    {"_BNPmix_MAR_mv_L",  (DL_FUNC) &_BNPmix_MAR_mv_L,  20},
    {"_BNPmix_MAR_mv_P",  (DL_FUNC) &_BNPmix_MAR_mv_P,  22},
    // others
    {"_BNPmix_BNPmix_psm",        (DL_FUNC) &_BNPmix_BNPmix_psm,        1},
    {"_BNPmix_clean_partition",   (DL_FUNC) &_BNPmix_clean_partition,   1},
    {"_BNPmix_BNPmix_VI_LB",      (DL_FUNC) &_BNPmix_BNPmix_VI_LB,      2},
    {"_BNPmix_BNPmix_BIN",        (DL_FUNC) &_BNPmix_BNPmix_BIN,        2},
    {NULL, NULL, 0}
};

void R_init_BNPmix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
