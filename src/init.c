#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "supclust.h"

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

// pelora.c --------------------------------------------

static R_NativePrimitiveArgType R_clusterer_typ[] = {
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    /* probini */ REALSXP, REALSXP, REALSXP, INTSXP,
    /* g */       INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    /*blockflag*/ INTSXP, INTSXP, INTSXP, INTSXP,
    /*kriterium*/ REALSXP
};

// wilma.c --------------------------------------------

static R_NativePrimitiveArgType R_multicluster_typ[] = {
    REALSXP, INTSXP,
    INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP,
    REALSXP,
    INTSXP, INTSXP, INTSXP,
    REALSXP,
    INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_margin_typ[] = {
    REALSXP, INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType R_score_typ[] = {
    REALSXP, INTSXP, INTSXP, REALSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(R_clusterer),
    CDEF(R_multicluster),
    CDEF(R_margin),
    CDEF(R_score),
    {NULL, NULL, 0}
};


void R_init_supclust(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, /* CallEntries */NULL, /* FortEntries */NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
