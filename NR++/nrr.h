//---------------------------------------------------------------------------

#ifndef NRRH
#define NRRH
//---------------------------------------------------------------------------
//RUTINAS
namespace nr {
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	 double *f, double stpmax, int *check, double (*func)(double []));
  void lubksb(double **a, int n,int *indx,double b[]);
  void ludcmp(double **a, int n,int *indx, double *d);
	void  simp1(double **a, int mm, int *ll, int nll, int iabf, int *kp,
		double *bmax);
	void  simp2(double **a, int n, int *l2, int nl2, int *ip, int kp,
		double *q1);
	void  simp3(double **a, int i1, int k1, int ip, int kp);
	void  simplx(double **a, int m, int n, int m1, int m2, int m3,
		int *icase, int *izrov, int *iposv);
//---------------------------------------------------------------------------
//NRUTILS
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//---------------------------------------------------------------------------
void free_ivector(int *v, long nl, long nh);//void free_ivector();
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_vector(double *v,int nl,int nh);
int *ivector(long nl, long nh);//int *ivector();
double **matrix(long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
double *vector(int nl,int nh);

//---------------------------------------------------------------------------
}
#endif
