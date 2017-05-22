//---------------------------------------------------------------------------

#ifndef FuncOptimH
#define FuncOptimH
//---------------------------------------------------------------------------
class Optimizacion {
  public:
    static void LBFGS(double p[], int n, double gtol, int *iter, double *fret,
    	double (*func)(double []), void (*dfunc)(double [], double []));
    static void LBFGSB(double p[], int n, double gtol, int *iter, double *fret,
    	double (*func)(double []), void (*dfunc)(double [], double []),
      double l[], double u[], int nbd[]);
};
#endif
