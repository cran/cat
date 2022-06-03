#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bipf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(estepc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(impc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ipf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(istepc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(llc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mgstepc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pstep1c)(void *, void *, void *, void *);
extern void F77_NAME(rngs)(void *);
extern void F77_NAME(sjn)(void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"bipf",    (DL_FUNC) &F77_NAME(bipf),    13},
    {"estepc",  (DL_FUNC) &F77_NAME(estepc),  15},
    {"impc",    (DL_FUNC) &F77_NAME(impc),    17},
    {"ipf",     (DL_FUNC) &F77_NAME(ipf),     12},
    {"istepc",  (DL_FUNC) &F77_NAME(istepc),  15},
    {"llc",     (DL_FUNC) &F77_NAME(llc),     15},
    {"mgstepc", (DL_FUNC) &F77_NAME(mgstepc), 18},
    {"pstep1c", (DL_FUNC) &F77_NAME(pstep1c),  4},
    {"rngs",    (DL_FUNC) &F77_NAME(rngs),     1},
    {"sjn",     (DL_FUNC) &F77_NAME(sjn),      4},
    {NULL, NULL, 0}
};

void R_init_cat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
