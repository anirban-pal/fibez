struct cBezier2 {
	//double PQB[48],PQB1[48];
  
  gsl_matrix * P = gsl_matrix_alloc (4, 3);
  gsl_matrix * Q = gsl_matrix_alloc (4, 3);

	gsl_matrix * PP = gsl_matrix_alloc (4, 4);
	gsl_matrix * QQ = gsl_matrix_alloc (4, 4);
  gsl_matrix * PQ = gsl_matrix_alloc (4, 4);
	
};

static int Integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void * params) {
  
  cBezier2 grparams = *(cBezier2 *) params; 
  
  double t1 = xx[0], t2 = xx[1];
  double t1c = 1.0-t1, t2c = 1.0-t2;

  double Xv[4] = {t1c*t1c*t1c, 3*t1c*t1c*t1, 3*t1c*t1*t1, t1*t1*t1};
  double Yv[4] = {t2c*t2c*t2c, 3*t2c*t2c*t2, 3*t2c*t2*t2, t2*t2*t2};

  double DXv[4] = {-3*t1c*t1c, 3*t1c*(t1c-2*t1), 3*t1*(2*t1c-t1), 3*t1*t1};
  double DYv[4] = {-3*t2c*t2c, 3*t2c*(t2c-2*t2), 3*t2*(2*t2c-t2), 3*t2*t2};
  
  double J1 = 0, J2 = 0, C = 0;

  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++){
      J1 += DXv[i]*gsl_matrix_get (grparams.PP, i, j)*DXv[j];
      J2 += DYv[i]*gsl_matrix_get (grparams.QQ, i, j)*DYv[j];
    
      C += Xv[i]*gsl_matrix_get (grparams.PP, i, j)*Xv[j] 
          + Yv[i]*gsl_matrix_get (grparams.QQ, i, j)*Yv[j] 
          -2*Xv[i]*gsl_matrix_get (grparams.PQ, i, j)*Yv[j];
      
    }
  
  C = sqrt(C);
  ff[0] = V(C)*sqrt(J1)*sqrt(J2);	

  return 0;
}

static int I1234(const int *ndim, const double xx[],
  const int *ncomp, double F[], void * params) {
  
  cBezier2 grparams = *(cBezier2 *) params; 
  
  double t1 = xx[0], t2 = xx[1];
  double t1c = 1.0-t1, t2c = 1.0-t2;

  double Xv[4] = {t1c*t1c*t1c, 3*t1c*t1c*t1, 3*t1c*t1*t1, t1*t1*t1};
  double Yv[4] = {t2c*t2c*t2c, 3*t2c*t2c*t2, 3*t2c*t2*t2, t2*t2*t2};

  double DXv[4] = {-3*t1c*t1c, 3*t1c*(t1c-2*t1), 3*t1*(2*t1c-t1), 3*t1*t1};
  double DYv[4] = {-3*t2c*t2c, 3*t2c*(t2c-2*t2), 3*t2*(2*t2c-t2), 3*t2*t2};
  
  double J1 = 0, J2 = 0, C = 0;

  //double F[48] = {0};
  double cP[12] = {0},cQ[12] = {0},JP[12] = {0},JQ[12] = {0};

  for(int i=0;i<4;i++)
    for(int j=0;j<3;j++){
      //cP[i*3+j] = 0.0; cQ[i*3+j] = 0.0; JP[i*3+j] = 0.0; JQ[i*3+j] = 0.0;
      
      for(int k=0;k<4;k++) {
        cP[i*3+j] += Xv[i] * ( Xv[k] * gsl_matrix_get (grparams.P, k, j) - Yv[k] * gsl_matrix_get (grparams.Q, k, j) ) ;
        cQ[i*3+j] += Yv[i] * ( Yv[k] * gsl_matrix_get (grparams.Q, k, j) - Xv[k] * gsl_matrix_get (grparams.P, k, j) ) ;
        
        JP[i*3+j] += DXv[i] * DXv[k] * gsl_matrix_get (grparams.P, k, j) ;
        JQ[i*3+j] += DYv[i] * DYv[k] * gsl_matrix_get (grparams.Q, k, j) ;
        
      }
    }
  
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++){
      J1 += DXv[i]*gsl_matrix_get (grparams.PP, i, j)*DXv[j];
      J2 += DYv[i]*gsl_matrix_get (grparams.QQ, i, j)*DYv[j];
    
      C += Xv[i]*gsl_matrix_get (grparams.PP, i, j)*Xv[j] 
          + Yv[i]*gsl_matrix_get (grparams.QQ, i, j)*Yv[j] 
          -2*Xv[i]*gsl_matrix_get (grparams.PQ, i, j)*Yv[j];
      
    }
  
  C = sqrt(C);  
  for(int i=0;i<12;i++)
  {
    F[i] = sqrt(J1*J2)*cP[i]*f(C);
    F[i+12] = sqrt(J2/J1)*JP[i]*V(C);
    F[i+24] = sqrt(J1*J2)*cQ[i]*f(C);
    F[i+36] = sqrt(J1/J2)*JQ[i]*V(C);
  }

  return 0;
}

static int dJP(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  
  loop(i,4) loop(j,3) {
    double sum = 0.0;
    loop(k,4) sum += DX[i]*DX[k]*gsl_matrix_get(params.P, k, j);
    ff[i*3+j] = sum/sqrt(Jsq);
  }
  
  return 0;
}

double inter_energy(cBezier *b1, cBezier *b2)
{
  cBezier2 bpair;
  
  gsl_matrix_memcpy(bpair.P, b1->P);
  gsl_matrix_memcpy(bpair.Q, b2->P);
  
  gsl_matrix_memcpy(bpair.PP, b1->PP);
  gsl_matrix_memcpy(bpair.QQ, b2->PP);
  
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, bpair.P, bpair.Q, 0.0, bpair.PQ); 
  
  //gsl_matrix_print(bpair.PP,"PP");
  //gsl_matrix_print(bpair.QQ,"QQ");
  //gsl_matrix_print(bpair.PQ,"PQ");
  
  
  int comp, nregions, neval, fail;
  double integral[1], error[1], prob[1];
  
  Cuhre(2, 1, Integrand, &bpair, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integral, error, prob);
  
  double integralc[48], errorc[48], probc[48];    
  Cuhre(2, 48, I1234, &bpair, NVEC,
      EPSREL, EPSABS, VERBOSE | LAST,
      MINEVAL, MAXEVAL, KEY,
      STATEFILE, SPIN,
      &nregions, &neval, &fail, integralc, errorc, probc);
  
  double rho1 = rho0[b1->type-1]*b1->length0/b1->length ;
  double rho2 = rho0[b2->type-1]*b2->length0/b2->length ;

  //double L12 = rho1*rho2;
  //L12 = 1.0;  

  double integralLP[12], errorLP[12], probLP[12];
  if(1) Cuhre(2, 12, dJP, b1, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integralLP, errorLP, probLP);
    
  double integralLQ[12], errorLQ[12], probLQ[12];
  if(1) Cuhre(2, 12, dJP, b2, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integralLQ, errorLQ, probLQ);
  
  double U12 = integral[0]*rho1*rho2;
  
  b1->cohE += 0.5*U12;
  b2->cohE += 0.5*U12;
  
  loop(i,12) {
    b1->f[i] += -(integralc[i]+integralc[i+12])*rho1*rho2 - 1.0*(U12/b1->length)*integralLP[i] ;
    b2->f[i] += -(integralc[i+24]+integralc[i+36])*rho1*rho2 - 1.0*(U12/b2->length)*integralLQ[i];
  }
  
  if(0) {
    double errorI = (integral[0]-b1->length0 * b2->length0)/(b1->length0 * b2->length0);
    printf("EPSREL=%le Length_product = %lf Integral = %lf Error = %le\n",EPSREL,b1->length0 * b2->length0, integral[0], errorI);
  } 
  
  return U12;
}
