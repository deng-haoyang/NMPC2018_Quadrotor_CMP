/*--------------------------------------------------------
  RHMAINU.C 
                du(tau(t),t)/dt  (Cf. #29-45~)
                Continuation in terminal time
                ODE: ADAMS/RKG 
                T.Ohtsuka        '96/12/29~ (rhsf_rv.c)
                                 '97/08/21~ (rvtr.c)
                                 '97/10/30~'97/10/31  (rhmain.c)
                                 '99/10/31~ (rhmainu.c)
--------------------------------------------------------*/

/*-------------- Initial Conditions -------------- */
void dhu0func(int dimu, double *du0, double *dhu)
{
	double du[DIMUC], u[DIMUC], hu[DIMUC];
	mmulsc(1,DIMUC, du0, hdir, du);
	madd(1,DIMUC, u0, du, u);
	hufunc(tsim0, x0, lmd0, u, hu);
	msub(1,DIMUC, hu, hu0, dhu);
	mdivsc(1,DIMUC, dhu, hdir, dhu);
}

void init(double x[], double u[])
{
    int i, j;
    double r, b[DIMUC], du0[DIMUC], erru0[DIMUC+1];
	lmd0= defvector(DIMX);
	hu0 = defvector(DIMUC);
    mmov(1,DIMX, x0, x);
    phix(tsim0, x0, lmd0);
	hufunc(tsim0, x0, lmd0, u0, hu0);
	for(i=0; i<DIMUC; i++) du0[i] = 0;
	r = sqrt(mvinner(DIMUC, hu0, hu0));
	i=0;
	while( (r > rtol) && (i < 100) ){
		for(j=0; j<DIMUC; j++) {
			b[j] = -hu0[j];
		}
		nfgmres(dhu0func, DIMUC, b, du0, DIMUC, erru0);
#ifdef TRACE_ON
		printf("i=%d ,  r=%g \n", i, r);
		printf("u0 = {%g", u0[0]);
		for(j=1; j<DIMUC; j++){
			printf(", %g", u0[j]);
		}
		printf("} \n");
		printf("erru0 = ");
		for(j=0; j<=DIMUC; j++){ 
			printf("%g ", (float)erru0[j]);
		}
		printf("\n");
#endif
		madd(1,DIMUC, u0, du0, u0);
		hufunc(tsim0, x0, lmd0, u0, hu0);
		r = sqrt(mvinner(DIMUC, hu0, hu0));
		i++;
	}
#ifdef TRACE_ON
	printf("i=%d ,  r=%g \n", i, r);
	printf("u0 = {%g", u0[0]);
	for(j=1; j<DIMUC; j++){
		printf(", %g", u0[j]);
	}
	printf("} \n\n");
#endif
	mmov(1,DIMUC, u0, u);
    for(i=0; i<dv; i++){
        mmov(1,DIMUC, u0, utau[i]);
		mmov(1,DIMUC, hu0, hutau[i]);
    }
	for(i=0; i<dimeq; i++){
			duvec[i] = 0;
	}
	freevector(lmd0);
	freevector(hu0);
}


/*-------------- Control Update -------------- */
void errfunc(double t, double *x, double **u, double **hu)
{
     int i;
     double tauf, taut, linp[DIMX+DIMUC];
#ifdef ADAMS
	 double f3[3][DIMX];
#endif
     tauf = tf * (1. - exp(-alpha * t) ); 
     htau = tauf / dv;
	 mmov(1,DIMX, x, xtau[0]);  
#ifdef ADAMS
	for(taut = t, i=0; (i<3)&&(i < dv); taut += htau, i++){
		nfrkginpex(xpfunc,taut,xtau[i],u[i],htau,DIMX,xtau[i+1],f3[i]);
		tau[i] = taut; 
	}
	for( ; i < dv; taut += htau, i++){
		nfadamsinp(xpfunc,taut,xtau[i],u[i],f3[2],f3[1],f3[0],htau,DIMX,xtau[i+1]);
		tau[i] = taut; 
	}
	tau[i] = taut; 
	phix(taut, xtau[dv], ltau[dv]);
	for(i = dv-1; (i > dv-4)&&(i >= 0); i--){
		mmov(1,DIMX, xtau[i], linp);
		mmov(1,DIMUC, u[i], linp+DIMX);
		nfrkginpex(lpfunc,taut,ltau[i+1],linp,-htau,DIMX,ltau[i],f3[dv-1-i]);
		taut -= htau; 
		hufunc(taut, xtau[i], ltau[i+1], u[i], hu[i]);
	}
	for( ; i >= 0; i--){
		mmov(1,DIMX, xtau[i], linp);
		mmov(1,DIMUC, u[i], linp+DIMX);
		nfadamsinp(lpfunc,taut,ltau[i+1],linp,f3[2],f3[1],f3[0],-htau,DIMX,ltau[i]);
		taut -= htau; 
		hufunc(taut, xtau[i], ltau[i+1], u[i], hu[i]);
	}
#else
	for(taut = t, i=0; i < dv; taut += htau, i++){
		nfeulerinp(xpfunc,taut,xtau[i],u[i],htau,DIMX,xtau[i+1]);
		tau[i] = taut; 
	}
	tau[i] = taut; 
	phix(taut, xtau[dv], ltau[dv]);
	for(i = dv-1; i >= 0; i--){
		mmov(1,DIMX, xtau[i], linp);
		mmov(1,DIMUC, u[i], linp+DIMX);
		nfeulerinp(lpfunc,taut,ltau[i+1],linp,-htau,DIMX,ltau[i]);
		taut -= htau; 
		hufunc(taut, xtau[i], ltau[i+1], u[i], hu[i]);
	}
#endif
}

void adufunc(int n, double *du, double *adu)
{
	mmulsc(1,n, du, hdir, dutmp);
	madd(1,n, utau[0], dutmp, utau1[0]);
	errfunc(ts, x1s, utau1, hutau2);
	msub(1,n, hutau2[0], hutau1[0], adu);
	mdivsc(1,n, adu, hdir, adu);
}

void unew(double t, double x[], double x1[], double u[])
{
	int i;
	double x2[DIMX];
	ts = t + hdir;
#ifdef HDIR_EQ_HT
 mmov(1,DIMX, x1, x1s);
#else
	mmulsc(1,DIMX, x, onemhdir, x2);
	mmulsc(1,DIMX, x1, hdirbht , x1s);
	madd(1,DIMX, x2, x1s, x1s);
#endif
	errfunc(t, x, utau, hutau);
#ifdef TRACE_ON
	printf("isim=%d, t=%g \n", isim, (float)t);
	for(i=0; i<dv; i++){
		printf("taut=%10g, x0=%10g, l0=%10g, u0=%10g, hu=%10g \n", 
			(float)tau[i], (float)xtau[i][0], (float)ltau[i][0], (float)utau[i][0], (float)hutau[i][0]);
	}
	printf("taut=%10g, x0=%10g, l0=%10g \n", 
			(float)tau[i], (float)xtau[i][0], (float)ltau[i][0]);
#endif
	errfunc(ts, x1s, utau, hutau1);
	mmulsc(1,dimeq, hutau[0], onemzetahdir, bvec);
	msub(1,dimeq, bvec, hutau1[0], bvec);
	mdivsc(1,dimeq, bvec, hdir, bvec);
	nfgmres(adufunc, dimeq, bvec, duvec, kmax, errvec);
	for(i=0; i<dimeq; i++){
		utau[0][i] += ht * duvec[i];
#ifdef RESET_DU
		duvec[i] = 0;   /**** Reset dU ****/
#endif
	}
	mmov(1,DIMUC, utau[0], u);
#ifdef TRACE_ON
	printf("errvec = ");
	for(i=0; i<=kmax; i++) printf("%g ", (float)errvec[i]);
	printf("\n\n");
#endif
}


/*-------------- Save Data -------------- */
void savedata(FILE *fpx, FILE *fpu, FILE *fpe, FILE *fpp, double t, double x[], double u[])
{
     int i, j;
     double hu2;
#ifdef P_IS_USED
	 double p[DIMX];
#endif
     
     for(i=0; i<DIMX; i++) fprintf(fpx, "%g ", (float)x[i]);
     fprintf(fpx, "\n");
     for(i=0; i<DIMUC; i++) fprintf(fpu, "%g ", (float)u[i]);
     fprintf(fpu, "\n");
     htau = tf * (1. - exp(-alpha * t) ) / dv; 
	 for(i=0; i<DIMUC; i++) {
		 for(hu2=0, j=0; j<dv; j++){ 
			 hu2 += hutau[j][i] * hutau[j][i];
		 }
		 fprintf(fpe, "%g ", (float)hu2);
	 }
	 fprintf(fpe, "\n");
#ifdef P_IS_USED
     pfunc(t, p);
     for(i=0; i<DIMP; i++) fprintf(fpp, "%g ", (float)p[i]);
     fprintf(fpp, "\n");
#endif
}


/*-------------- Main -------------- */
int main(void)
{
     int nsim; /* isim is a global var. */ 
     double t, x[DIMX], x1[DIMX], u[DIMUC], f3[3][DIMX];
     FILE *fpx, *fpu, *fpe, *fpp, *fpc;  
     clock_t tstart, tend, tbefore, tafter, ttotal=0;
     float t_cpu, t_s2e; 

     nsim  = (int)((tsim -tsim0) / ht);
#ifdef HDIR_EQ_HT
     if( hdir != ht ){
          printf("*** Warning: hdir is set equal to ht. \n");
     }
     hdir = ht;
#endif
	 /* Global Var. */
	 dimeq = dv * DIMUC;
	 onemhdir = 1 - hdir / ht;
	 hdirbht = hdir / ht; 
     onemzetahdir = 1 - zeta * hdir; 
	 /* Global Arrays */
	tau   = defvector(dv+1); 
    xtau  = defmatrix(dv+1, (int)DIMX);
	xtau1 = defmatrix(dv+1, (int)DIMX);
	ltau  = defmatrix(dv+1, (int)DIMX);
	utau  = defmatrix(dv, (int)DIMUC);
	utau1 = defmatrix(dv, (int)DIMUC); 
	hutau = defmatrix(dv, (int)DIMUC);
	hutau1= defmatrix(dv, (int)DIMUC);
	hutau2= defmatrix(dv, (int)DIMUC);
	x1s   = defvector(DIMX);
	bvec  = defvector(dimeq);
	duvec = defvector(dimeq);
	dutmp = defvector(dimeq);
	errvec= defvector(kmax+1);
     
	 /* Open Files */
#define MATLAB_M_FILE

#ifdef MATLAB_M_FILE
     fpx = fopen(FNMHD "x.m", "w");
     fpu = fopen(FNMHD "u.m", "w");
     fpe = fopen(FNMHD "e.m", "w");
#ifdef P_IS_USED
     fpp = fopen(FNMHD "p.m", "w");
#else
     fpp = fpx;
#endif
     fpc = fopen(FNMHD "c.m", "w");
#else
     fpx = fopen(FNMHD "x.dat", "w");
     fpu = fopen(FNMHD "u.dat", "w");
     fpe = fopen(FNMHD "e.dat", "w");
#ifdef P_IS_USED
     fpp = fopen(FNMHD "p.dat", "w");
#else
     fpp = fpx;
#endif
     fpc = fopen(FNMHD "c.dat", "w");
#endif
         
     if( (fpx==NULL)||(fpu==NULL)||(fpe==NULL)||(fpp==NULL)||(fpc==NULL) ){
          printf(" ****** can't open data files \n");
          exit(1);
     }
	 /* Start Simulation */
     printf(" Start  \n");
     tstart = clock();
     init(x,u);
     for(t = tsim0, isim=0; (isim < 3)&&(isim < nsim); t+=ht, isim++){
          if(isim%dstep == 0 ) savedata(fpx, fpu, fpe, fpp, t, x, u);
          nfrkginpex(xpfunc,t,x,u,ht,DIMX, x1,f3[isim]);
          tbefore = clock();
          unew(t,x,x1,u);
          tafter = clock();
          ttotal += tafter - tbefore;
		  mmov(1,DIMX, x1, x);
     }
     for( ; isim < nsim; t+=ht, isim++) {
          if(isim%dstep == 0 ) savedata(fpx, fpu, fpe, fpp, t, x, u);
          nfadamsinp(xpfunc,t,x,u,f3[2],f3[1],f3[0],ht,DIMX, x1);
          tbefore = clock();
          unew(t,x,x1,u);
          tafter = clock();
          ttotal += tafter - tbefore;
		  mmov(1,DIMX, x1, x);
    }
     if(isim%dstep == 0 ) savedata(fpx, fpu, fpe, fpp, t, x, u);
     tend = clock();
     printf(" End \n");
     printf("CPU Time (Only for Simulation): %g sec\n", t_cpu = ((float)(TDIV*ttotal/(CLOCKS_PER_SEC)))/TDIV );
     printf("Start ~ End (Simulation + Data Save) : %g sec\n", t_s2e = ((float)(TDIV*(tend-tstart)/(CLOCKS_PER_SEC)))/TDIV );
     printf("CLOCKS_PER_SEC = %ld\n", (long)CLOCKS_PER_SEC);
	 /* Close Files */
     fclose(fpx);
     fclose(fpu);
     fclose(fpe);
#ifdef P_IS_USED
     fclose(fpp);
#endif
     final(fpc, t_cpu, t_s2e);
     fclose(fpc);
	 /* Del. Global Arrays */
	 freevector(tau);
	 freematrix(xtau);
	 freematrix(xtau1);
	 freematrix(ltau);
	 freematrix(utau);
	 freematrix(utau1);
	 freematrix(hutau);
	 freematrix(hutau1);
	 freematrix(hutau2);
	 freevector(x1s);
	 freevector(bvec);
	 freevector(duvec);
	 freevector(dutmp);
	 freevector(errvec);
}
