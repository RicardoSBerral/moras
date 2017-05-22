//---------------------------------------------------------------------------



#include "nrr.h"

//---------------------------------------------------------------------------
#include <math.h>
//---------------------------------------------------------------------------
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <exception>
//using namespace std;
//---------------------------------------------------------------------------

namespace nr
{

static float sqrarg;
//static double dsqrarg;
//static double dmaxarg1,dmaxarg2;
//static double dminarg1,dminarg2;
static float maxarg1,maxarg2;
//static float minarg1,minarg2;
//static long lmaxarg1,lmaxarg2;
//static long lminarg1,lminarg2;
//static int imaxarg1,imaxarg2;
//static int iminarg1,iminarg2;
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//RUTINAS
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []))
{
	void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	int check,i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double *dg,*g,*hdg,**hessin,*pnew,*xi;

	dg=vector(1,n);
	g=vector(1,n);
	hdg=vector(1,n);
	hessin=matrix(1,n,1,n);
	pnew=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g);
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac > sqrt(EPS*sumdg*sumxi)) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=i;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
					hessin[j][i]=hessin[i][j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("too many iterations in dfpmin");
	FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#define EPS 1.0e-6
#define TINY 1.0e-20;
#define FREEALL free_ivector(l3,1,m);free_ivector(l2,1,m);\
	free_ivector(l1,1,n+1);
#define NR_END 1
#define FREE_ARG char*
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(double []))
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam > 0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		alam=FMAX(tmplam,0.1*alam);
	}
}
#undef ALF
#undef TOLX
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void lubksb(double **a, int n,int *indx,double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
//---------------------------------------------------------------------------
void ludcmp(double **a, int n,int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}
//---------------------------------------------------------------------------
//void simplx(a,m,n,m1,m2,m3,icase,izrov,iposv)
//int m,n,m1,m2,m3,*icase,izrov[],iposv[];
//double **a;
void  simplx(double **a, int m, int n, int m1, int m2, int m3,
		int *icase, int *izrov, int *iposv)
{
	int i,ip,ir,is,k,kh,kp,m12,nl1,nl2;
	int *l1,*l2,*l3;//,*ivector();
	double q1,bmax;
//	void simp1(),simp2(),simp3(),free_ivector();

	if (m != (m1+m2+m3)) return;//throw new exception("Bad input constraint counts in SIMPLX");
	l1=ivector(1,n+1);
	l2=ivector(1,m);
	l3=ivector(1,m);
	nl1=n;
	for (k=1;k<=n;k++) l1[k]=izrov[k]=k;
	nl2=m;
	for (i=1;i<=m;i++) {
		if (a[i+1][1] < 0.0) return;//throw new exception("Bad input tableau in SIMPLX");
		l2[i]=i;
		iposv[i]=n+i;
	}
	for (i=1;i<=m2;i++) l3[i]=1;
	ir=0;
	if (m2+m3) {
		ir=1;
		for (k=1;k<=(n+1);k++) {
			q1=0.0;
			for (i=m1+1;i<=m;i++) q1 += a[i+1][k];
			a[m+2][k] = -q1;
		}
		do {
			simp1(a,m+1,l1,nl1,0,&kp,&bmax);
			if (bmax <= EPS && a[m+2][1] < -EPS) {
				*icase = -1;
				FREEALL return;
			} else if (bmax <= EPS && a[m+2][1] <= EPS) {
				m12=m1+m2+1;
				if (m12 <= m) {
					for (ip=m12;ip<=m;ip++) {
						if (iposv[ip] == (ip+n)) {
							simp1(a,ip,l1,
								nl1,1,&kp,&bmax);
							if (bmax > 0.0)
								goto one;
						}
					}
				}
				ir=0;
				--m12;
				if (m1+1 <= m12)
					for (i=m1+1;i<=m12;i++)
						if (l3[i-m1] == 1)
							for (k=1;k<=n+1;k++)
								a[i+1][k] = -a[i+1][k];
				break;
			}
			simp2(a,n,l2,nl2,&ip,kp,&q1);
			if (ip == 0) {
				*icase = -1;
				FREEALL return;
			}
one:		simp3(a,m+1,n,ip,kp);
			if (iposv[ip] >= (n+m1+m2+1)) {
				for (k=1;k<=nl1;k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k;is<=nl1;is++) l1[is]=l1[is+1];
				a[m+2][kp+1] += 1.0;
				for (i=1;i<=m+2;i++) a[i][kp+1] = -a[i][kp+1];
			} else {
				if (iposv[ip] >= (n+m1+1)) {
					kh=iposv[ip]-m1-n;
					if (l3[kh]) {
						l3[kh]=0;
						a[m+2][kp+1] += 1.0;
						for (i=1;i<=m+2;i++)
							a[i][kp+1] = -a[i][kp+1];
					}
				}
			}
			is=izrov[kp];
			izrov[kp]=iposv[ip];
			iposv[ip]=is;
		} while (ir);
	}
	for (;;) {
		simp1(a,0,l1,nl1,0,&kp,&bmax);
		if (bmax <= 0.0) {
			*icase=0;
			FREEALL return;
		}
		simp2(a,n,l2,nl2,&ip,kp,&q1);
		if (ip == 0) {
			*icase=1;
			FREEALL return;
		}
		simp3(a,m,n,ip,kp);
		is=izrov[kp];
		izrov[kp]=iposv[ip];
		iposv[ip]=is;
	}
}

//void simp1(a,mm,ll,nll,iabf,kp,bmax)
//double **a,*bmax;
//int mm,ll[],nll,iabf,*kp;
void  simp1(double **a, int mm, int *ll, int nll, int iabf, int *kp,
		double *bmax)
{
	int k;
	double test;

	*kp=ll[1];
	*bmax=a[mm+1][*kp+1];
	for (k=2;k<=nll;k++) {
		if (iabf == 0)
			test=a[mm+1][ll[k]+1]-(*bmax);
		else
			test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
		if (test > 0.0) {
			*bmax=a[mm+1][ll[k]+1];
			*kp=ll[k];
		}
	}
}

//void simp2(a,n,l2,nl2,ip,kp,q1)
//int n,l2[],nl2,*ip,kp;
//double **a,*q1;
void  simp2(double **a, int n, int *l2, int nl2, int *ip, int kp,
		double *q1)
{
	int k,ii,i;
	double qp,q0,q;

	*ip=0;
	for (i=1;i<=nl2;i++) {
		if (a[l2[i]+1][kp+1] < -EPS) {
			*q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
			*ip=l2[i];
			for (i=i+1;i<=nl2;i++) {
				ii=l2[i];
				if (a[ii+1][kp+1] < -EPS) {
					q = -a[ii+1][1]/a[ii+1][kp+1];
					if (q < *q1) {
						*ip=ii;
						*q1=q;
					} else if (q == *q1) {
						for (k=1;k<=n;k++) {
							qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
							q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
							if (q0 != qp) break;
						}
						if (q0 < qp) *ip=ii;
					}
				}
			}
		}
	}
}

//void simp3(a,i1,k1,ip,kp)
//int i1,k1,ip,kp;
//double **a;
void  simp3(double **a, int i1, int k1, int ip, int kp)
{
	int kk,ii;
	double piv;

	piv=1.0/a[ip+1][kp+1];
	for (ii=1;ii<=i1+1;ii++)
		if (ii-1 != ip) {
			a[ii][kp+1] *= piv;
			for (kk=1;kk<=k1+1;kk++)
				if (kk-1 != kp)
					a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
		}
	for (kk=1;kk<=k1+1;kk++)
		if (kk-1 != kp) a[ip+1][kk] *= -piv;
	a[ip+1][kp+1]=piv;
}
#undef EPS
#undef TINY
#undef FREEALL


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//NRUTILS
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
void free_vector (double *v,int nl,int nh)
{
	free((char*) (v+nl));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) return 0;//nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
 //	throw exception();//runtime_error(error_text);
}
double *vector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}


}//FIN namespace nr
