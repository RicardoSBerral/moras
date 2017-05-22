//---------------------------------------------------------------------------

#pragma hdrstop

#include "FuncOptim.h"
#include "nrr.h"
#include <stdlib.h>
extern "C" {
#include "f2c.h"
//#include "routines.c"
}

//---------------------------------------------------------------------------

#pragma package(smart_init)

doublereal *ADoblePrec(double f[], int tam) {
  doublereal *d = (doublereal*)malloc(tam * sizeof(doublereal));
  for (int i=0;i<tam;i++) d[i] = f[i];
  return d;
}

void Optimizacion::LBFGS(double p[], int n, double gtol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
    nr::dfpmin(p, n, gtol, iter, fret, func, dfunc);
}
    extern/* extern Subroutine*/ int s_copy(char *, char *, ftnlen, ftnlen);
  extern  /*extern*/ integer s_cmp(char *, char *, ftnlen, ftnlen);
void Optimizacion::LBFGSB(double p[], int nn, double gtol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []),
  double ll[], double uu[], int nnbd[])
{
    extern /* Subroutine */ int setulb_(integer *, integer *, doublereal *,
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, char *,
	    integer *, char *, logical *, integer *, doublereal *, ftnlen,
	    ftnlen);
    bool salir;
    doublereal f;
    doublereal *g = (doublereal *)malloc(sizeof(doublereal)*nn);
    doublereal *l;
    integer m, n, *nbd;
    doublereal *u, *x, t1, t2, wa[42584];
    integer iwa[3072];
    char task[60];
    doublereal factr;
    char csave[60];
    doublereal dsave[29];
    integer isave[44];
    logical lsave[4];
    doublereal pgtol;
    integer iprint;

    iprint = 1;
    factr = 1e7;
    pgtol = gtol;
    m = 5;
    n = nn;
    x = p+1;//ADoblePrec(p+1, n);
    u = ADoblePrec(uu+1, n);
    l = ADoblePrec(ll+1, n);
    nbd = (integer*)nnbd;

    s_copy(task, "START", (ftnlen)60, (ftnlen)5);
    do {
      salir = false;
      setulb_(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
	                  iprint, csave, lsave, isave, dsave, (ftnlen)60, (ftnlen)60);
      if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) {
/*        the minimization routine has returned to request the */
/*        function f and gradient g values at the current x. */
/*        Compute function value f for the sample problem. */
/* Computing 2nd power */
        f=func(x);
        dfunc(x, g);
    }
    else if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) {
    }
    else {
      salir = true;
    }
  } while (!salir);
  *iter = isave[29];
  *fret = f;
}


