#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<time.h>
#define TDIV 100  /* Minimal resolution of computation time is 1/TDIV */

/*-----------------------------------------------------------------------
     Some Fundamental Functions for RHC 
                     T. Ohtsuka  '97/10/30~'97/10/31 (rhfunc.c)
								 '00/01/27 (rhfuncu.c)
-----------------------------------------------------------------------*/

/*---- a[m][n] -> b[m][n] ----*/
void mmov( int m, int n, double *a, double *b )
{
     int	i, j, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
               tmp = i*n + j ;
               b[tmp] = a[tmp] ;
          }
}

/*---- a[m][n] + b[m][n] -> c[m][n] ----*/
void madd( int m, int n, double *a, double *b, double *c )
{
     int	i, j, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++) {
               tmp = i*n + j ;
               c[tmp] = a[tmp] + b[tmp] ;
          }
}

/*---- a[m][n] - b[m][n] -> c[m][n] ----*/
void msub( int m, int n, double *a, double *b, double *c )
{
     int	i, j, tmp ;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++) {
               tmp = i*n + j ;
               c[tmp] = a[tmp] - b[tmp] ;
          }
}

/*---- k * a[m][n] -> b[m][n] ----*/
void mmulsc( int m, int n, double *a, double k, double *b )
{
     int	i, j, tmp;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
               tmp = i*n + j ;
               b[tmp] = k * a[tmp];
          }
}

/*---- a[m][n] / k -> b[m][n] ----*/
void mdivsc( int m, int n, double *a, double k, double *b )
{
     int	i, j, tmp;
     for(i = 0; i < m; i++)
          for(j = 0; j < n; j++){
               tmp = i*n + j ;
               b[tmp] = a[tmp] / k;
          }
}

/*---- Inner Product of a[m] and b[m] ----*/
double	mvinner( int m, double *a, double *b )
{
	int	i;
	double tmp;
	tmp =0.;
	for(i=0; i<m; i++)
		tmp += a[i] * b[i];
	return tmp;
}

/*---- Define an n-Dimensional Vector ----*/
double *defvector(int n)
{
     double *v;
     v = (double *)malloc( (size_t)( n * sizeof(double) ) );
     if(!v) {
          printf("defvector : allocation failure. \n");
          exit(1);
     }
     return v;
}

void freevector(double *v)
{
	free(v);
}

/*---- Define an n by m Matrix ----*/
double **defmatrix(int n, int m)
{
     int i;
     double **a;

     a = (double **)malloc( (size_t)( n * sizeof(double*) ) );
     if(!a) {
          printf("defmatrix() : allocation failure. \n");
          exit(1);
     }
     a[0] = (double *)malloc( (size_t)( n * m * sizeof(double) ) );
     if(!a[0]) {
          printf("defmatrix() : allocation failure. \n");
          exit(1);
     }
     for(i=0; i<n-1; i++) a[i+1] = a[i] + m;
     return a;
}

void freematrix(double **a)
{
	free(a[0]);
	free(a);
}



/*---- Search an Ordered Table ----*/
/* xx: n vector.  Find j such that xx[j] <= x < xx[j+1]  */

#define EPSLOCATE 1e-6

int locate(int n, double xx[], double x)
{
     int j, ju, jm, jl, ascnd;

     jl = -1;
     ju = n;
     ascnd = (xx[n-1] >= xx[0]);
     if( (ascnd&&(x < xx[0]-EPSLOCATE)) ||(!ascnd&&(x > xx[0]+EPSLOCATE)) ){
        /*  printf(" locate: out of range of table \n"); */
         return 0;
     }
     if( (ascnd&&(x > xx[n-1]+EPSLOCATE)) ||(!ascnd&&(x < xx[n-1]-EPSLOCATE)) ){
        /*  printf(" locate: out of range of table \n"); */
         return n-2;
     }
	if( (ascnd&&(x < xx[0])) || (!ascnd&&(x > xx[0])) ) return 0;
	if( (ascnd&&(x > xx[n-1])) || (!ascnd&&(x < xx[n-1])) ) return n-2;
/*     if( (ascnd&&( (x < xx[0])||(x > xx[n-1]) ))
	|| ( !ascnd&&( (x > xx[0])||(x < xx[n-1]) )) ){
          printf(" locate: out of range of table \n");
         return 0;
     } */
     
     while( ju-jl > 1){
          jm = (ju + jl)/2 ;
          if( (x >= xx[jm]) == ascnd)
               jl = jm;
          else
               ju = jm;
     }
     if( x == xx[0]) j = 0;
     if( x == xx[n-1]) j = n-2; 
     else j = jl;
    return j;
}


/*---- Linear Interpolation of x(t) ----*/
/* t: n vector, x: n by m matrix,  x1 = x(t1) to be determined */

#define EPSTABLE1 1e-6

void table1(int n, int m, double t[], double **x, double t1, double x1[])
{
     int j;
     double thigh, tlow, alpha, beta;
     double *vtemp1, *vtemp2;
     vtemp1 = defvector(m);
     vtemp2 = defvector(m);

     j = locate(n, t, t1);
     tlow = t[j];
     thigh = t[j+1];
     if( fabs(thigh - tlow) < EPSTABLE1 ){
       madd(1,m, x[j], x[j+1], vtemp1);
       mmulsc(1,m, vtemp1, 0.5, x1);
     }
     else {
       alpha = (thigh - t1) / (thigh - tlow);
       beta = 1. - alpha;
       mmulsc(1,m, x[j], alpha, vtemp1);
       mmulsc(1,m, x[j+1], beta, vtemp2);
       madd(1,m, vtemp1, vtemp2, x1);
     }
     free(vtemp1);
     free(vtemp2);
}



/*---- k Givens Roations used in GMRES ----*/
void givrot(int k, double *c, double *s, double *v)
{
	int i;
	double w1, w2;
	for(i=0; i<k; i++){
		w1 = c[i] * v[i] - s[i] * v[i+1];
		w2 = s[i] * v[i] + c[i] * v[i+1];
		v[i] = w1;
		v[i+1] = w2;
	}
}

/*--------------------------------------------------------
 Linear Equation Subroutine
	GMRES (Generalized Minimum Residual) 
	Cf. C.T.Kelley: Iterative Methods for Linear and Nonlinear Equations
						T.Ohtsuka  '00/01/24
	axfunc(int n, double *x, double *ax): a function that gives A*x
	n: dim x
	b: the right-hand-side of the equation, A*x = b
	x: initial guess, to be replaced with the solution
	kmax: number of iteration
	err: residual norms, |b - A*x_i| (i=1,...,n+1)
----------------------------------------------------------*/
void nfgmres(void (*axfunc)(), int n, double *b, double *x, int kmax, double *err)
{
	int i,j,k;
	double rho, nu; 
	double *cvec, *svec, *gvec, *tmpvec, **hmat, **vmat; 

	cvec = defvector(kmax+1); 
	svec = defvector(kmax+1); 
	gvec = defvector(kmax+1);
	tmpvec = defvector(n); 
	hmat = defmatrix(kmax+1,kmax+1); 
	vmat = defmatrix(kmax+1, n); 

	axfunc(n, x, tmpvec); 
	msub(1,n, b, tmpvec, tmpvec); 
	rho = sqrt(mvinner(n, tmpvec, tmpvec));
	gvec[0] = rho; 
	for(i=1; i<kmax+1; i++){
		gvec[i] = 0;
	}
	err[0] = rho;

	mdivsc(1,n, tmpvec, rho, vmat[0]);
	for(k=0; k<kmax; k++){
		axfunc(n, vmat[k], vmat[k+1]); 

		/* Modified Gram-Schmidt */
		for(j=0; j<=k; j++){
			hmat[k][j] = mvinner(n, vmat[j], vmat[k+1]);
			mmulsc(1,n, vmat[j], hmat[k][j], tmpvec);
			msub(1,n, vmat[k+1], tmpvec, vmat[k+1]); 
		}
		hmat[k][k+1] = sqrt(mvinner(n, vmat[k+1], vmat[k+1]));
		
		/* No Breakdown? */
		if( hmat[k][k+1] != 0){
			mdivsc(1,n, vmat[k+1], hmat[k][k+1], vmat[k+1]);
		}
		else{
			printf("gmress() : breakdown \n");
		}
		
		/* Givens Rotation */
		givrot(k, cvec, svec, hmat[k]); 
		nu = sqrt(hmat[k][k] * hmat[k][k] + hmat[k][k+1] * hmat[k][k+1]);
		if( nu != 0 ){
			cvec[k] = hmat[k][k] / nu;
			svec[k] = - hmat[k][k+1] / nu;
			hmat[k][k] = cvec[k] * hmat[k][k] - svec[k] * hmat[k][k+1];
			hmat[k][k+1] = 0;
			givrot(1, cvec+k, svec+k, gvec+k);
		}
		else printf("nu is zero!\n");
		
		/* Residual Update */
		rho = fabs(gvec[k+1]);
		err[k+1] = rho;
	}
	
	/* Solve hmat * y = gvec (hmat: upper triangular) */
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j][i] * cvec[j];
		}
		cvec[i] = nu / hmat[i][i] ; 
/*		for(nu=0, j=i+1; j<k; j++){
			nu += hmat[j][i] * cvec[j];
		}
		cvec[i] = (gvec[i] - nu) / hmat[i][i] ;*/
	}
	/* Ans. */
	for(i=0; i<n; i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}
	
	freevector(cvec);
	freevector(svec);
	freevector(gvec);
	freevector(tmpvec);
	freematrix(hmat);
	freematrix(vmat);
}



/*------------------------------------------------------------------------
			ODE.C: nfrkginpex() and nfadamsinp() for NLSF/SCM
                  	T.Ohtsuka  	 '90/09/26
								~'91/07/06 (for UNIX (Sun))
								 '92/10/29 ( adams() )
								 '93/06/25 ode.c
								 '99/12/19 with input u[]
-------------------------------------------------------------------------*/

#define	DIMRK	50	/* Maximum Dimension of the state vector (for the Runge-Kutta_Gill Method) */
#define	DIMAD	50	/* Maximum Dimension of the state vector (for the Adams Method) */
#define	DIMEU	50	/* Maximum Dimension of the state vector (for the Euler Method) */

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Runge-Kutta_Gill Method	(Original by T.Murayama)

---- Variation of rkg() : dx/dt is also returned. ----
----  (for Start of the Adams method)  ----
----------------------------------------------------------*/
void nfrkginpex(func,x,y,u,h,dim,ans,fxy)
void (*func)();
double x,y[],u[],h,ans[],fxy[]; 
int dim;
{
	int i;
	double fval[DIMRK],k1[DIMRK],k2[DIMRK],k3[DIMRK],k4[DIMRK],
			yp1[DIMRK],yp2[DIMRK],yp3[DIMRK],
			q1[DIMRK],q2[DIMRK],q3[DIMRK],
			c1 = 0.2928932188134528,
			c2 = 0.1213203435596426,
			c3 = 0.5857864376269054,
			c4 = 1.707106781186548,
			c5 = 4.121320343559643,
			c6 = 3.414213562373097;

/*	fval = (double *)malloc( (size_t)( dim * sizeof(double) ) );
	k1   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
	k2   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
	k3   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
	k4   = (double *)malloc( (size_t)( dim * sizeof(double) ) );
*/

	func(x,y,u,fxy);
	for(i = 0; i < dim; i++)
	{
		 k1[i] = h * fxy[i];
		yp1[i] = y[i] + 0.5 * k1[i];
		 q1[i] = k1[i];
	}
	func(x + 0.5 * h,yp1,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k2[i] = h * fval[i];
		yp2[i] = yp1[i] + c1 * (k2[i] - q1[i]);;
		 q2[i] = c2 * q1[i] + c3 * k2[i];
	}
	func(x + 0.5 * h,yp2,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k3[i] = h * fval[i];
		yp3[i] = yp2[i] + c4 * (k3[i] - q2[i]);
		 q3[i] = -c5 * q2[i] + c6 * k3[i];
	}
	func(x + h,yp3,u,fval);
	for(i = 0; i < dim; i++)
	{
		 k4[i] = h * fval[i];
		ans[i] = yp3[i] + k4[i] / 6.0 - q3[i] / 3.0;
	}
}

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Adams Method (Predictor-Corrector Method)
		Predictor: Adams-Bashforth
		Corrector: Adams-Moulton
				T.Ohtsuka  '92/10/24
----------------------------------------------------------*/
void nfadamsinp(func,x,y,u,f1,f2,f3,h,dim,ans)
void (*func)();
double x,y[],u[],f1[],f2[],f3[],h,ans[];  /* fi: func(x,y) at i-steps ago. */
int dim;						   /* *Notice! fi are updated. */
{
	int i;
	double y1[DIMAD], fval[DIMAD], cd;
/*	double *y1, *fval, cd;

	fval = (double *)malloc( (size_t)( dim * sizeof(double) ) );
	y1   = (double *)malloc( (size_t)( dim * sizeof(double) ) ); */
	cd = h / 24.;

/*---- Adams-Bashforth Predictor ----*/
	func(x,y,u,fval);
	for(i=0; i<dim; i++)
	{
		y1[i] = 55.* fval[i] -59.* f1[i] +37.* f2[i] -9.* f3[i];
		y1[i] *= cd;
		y1[i] += y[i];
		f3[i] = f2[i]; f2[i] = f1[i]; f1[i] = fval[i]; /* Shift */
	}

/*---- Adams-Moulton Corrector ----*/
	func(x+h,y1,u,fval);
	for(i=0; i<dim; i++)
	{
		y1[i] = 9.* fval[i] +19.* f1[i] -5.* f2[i] + f3[i];
		y1[i] *= cd;
		ans[i] = y1[i] + y[i];
	}

/*	free(fval);
	free(y1); */
}

/*--------------------------------------------------------
 Simultaneous Ordinaly Diferential Equation Subroutine
	Euler Method (Forward Differentce)
				T.Ohtsuka  '99/12/19
----------------------------------------------------------*/
void nfeulerinp(void (*func)(), double x, double y[], double u[], double h, int dim, double ans[])
{
	int i;
	double fval[DIMEU];

	func(x,y,u,fval);
	for(i=0; i<dim; i++)
		ans[i] = y[i] + h * fval[i];
}


/*-------------- Global Variagles Defined in Program -------------- */

int	   isim;
int		dimeq;
double onemhdir, hdirbht, onemzetahdir;
double htau;
double ts;
double tauf;
double *tau;
double **xtau, **xtau1, **ltau;
double **utau, **utau1, **hutau, **hutau1, **hutau2;
double *x1s, *bvec, *duvec, *dutmp, *errvec;
double *lmd0, *hu0;
