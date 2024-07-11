void dH(const double q[], const double p[], void *fn_data, int cflag)
{
  fibnetwork *FN = (fibnetwork *) fn_data;

  int np = FN->np;	int nb = FN->nb;  int nf = FN->nf;
  
  loop(i,np) loop(j,3) if(!isnan(q[i*3+j])) FN->CPx[i].comp[j] = q[i*3+j];
  loop(i,np) loop(j,3) if(!isnan(p[i*3+j])) FN->CPw[i].comp[j] = p[i*3+j]; //define momenta from momenta vector, define velocity later

  loop(i,nb) {
    loop(j,4) {
      int cpid = FN->Bz_list[i].CP[j]-1;
      loop(k,3) {
        FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
        FN->Bz_list[i].w[j*3+k] = FN->CPw[cpid].comp[k];
	  }
    }
    FN->Bz_list[i].update_bezier(1);
  }
  
  loop(i,np) loop(j,3) {
    gsl_matrix_set(FN->Psys, i,j, FN->CPx[i].comp[j]);
    gsl_matrix_set(FN->Wsys, i,j, FN->CPw[i].comp[j]);
  }
  
  FN->compute_dH();
  
  gsl_vector *wsys = gsl_vector_alloc (3*np);
  gsl_vector *vsys = gsl_vector_alloc (3*np);
  
  loop(i,3*np) {
    gsl_vector_set(wsys, i, p[i]);
  }
  
  gsl_matrix *GMtmp1 = gsl_matrix_alloc (3*np, 3*np);
  gsl_matrix *GMtmp2 = gsl_matrix_alloc (3*np, 3*np);
  
  gsl_matrix_memcpy(GMtmp1, FN->GM3);
   
  gsl_permutation * p1 = gsl_permutation_alloc (3*np); int s;
  gsl_linalg_LU_decomp (GMtmp1, p1, &s);
  gsl_linalg_LU_solve (GMtmp1, p1, wsys, vsys);
  
  //determine quantites for continuity constraints
  int ncc_active = 0;
  std::vector<iVector5> active_cps;
  if(1) loop(i,FN->C2_cps.size()) {
      int bi = FN->C2_cps[i].comp[0];
      int ci = FN->C2_cps[i].comp[1];
      int di = FN->C2_cps[i].comp[2];
      int ei = FN->C2_cps[i].comp[3];
      int fi = FN->C2_cps[i].comp[4];
      
      active_cps.push_back(iVector5(bi,ci,di,ei,fi));
      ncc_active++;
  }
  
  int C1 = C1_flag, C2 = C2_flag;
  
  int ncon_active = FN->nfdof;
  
  if(C1==1) ncon_active += 3*ncc_active;
  if(C2==1) { ncon_active += 3*ncc_active; ncc_active = 2*ncc_active; }
  
  gsl_matrix *Lambda = gsl_matrix_alloc (3*ncc_active, 3*np);
  gsl_matrix *Gamma = gsl_matrix_alloc (ncon_active, 1);
  
  gsl_matrix *Lsys = gsl_matrix_alloc (ncon_active, 3*np);
  gsl_matrix *bsys = gsl_matrix_alloc (ncon_active, 1);
  
  gsl_matrix *AMA = gsl_matrix_alloc (ncon_active, ncon_active);
  gsl_matrix *AHq = gsl_matrix_alloc (ncon_active, 1);
  
  gsl_matrix_set_zero(Lambda);
  gsl_matrix_set_zero(Gamma);
  
  if(C1) loop(c,active_cps.size()) {
      int bi = FN->C2_cps[c].comp[0];
      int ci = FN->C2_cps[c].comp[1];
      int di = FN->C2_cps[c].comp[2];
      int ei = FN->C2_cps[c].comp[3];
      int fi = FN->C2_cps[c].comp[4];
            
      double q[15] = {FN->CPx[bi-1].comp[0], FN->CPx[bi-1].comp[1], FN->CPx[bi-1].comp[2],
                                  FN->CPx[ci-1].comp[0], FN->CPx[ci-1].comp[1], FN->CPx[ci-1].comp[2],
                                  FN->CPx[di-1].comp[0], FN->CPx[di-1].comp[1], FN->CPx[di-1].comp[2],
                                  FN->CPx[ei-1].comp[0], FN->CPx[ei-1].comp[1], FN->CPx[ei-1].comp[2],
                                  FN->CPx[fi-1].comp[0], FN->CPx[fi-1].comp[1], FN->CPx[fi-1].comp[2]  };
      
      double u[3] = { 3*(q[6]-q[3]), 3*(q[7]-q[4]), 3*(q[8]-q[5]) };
      double v[3] = { 3*(q[9]-q[6]), 3*(q[10]-q[7]), 3*(q[11]-q[8]) };
      
      double uu = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
      double vv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
      
      double squu = sqrt(uu), sqvv = sqrt(vv);
      
      double du[3][15] = {0}, dv[3][15] = {0};
      double udu[15] = {0}, vdv[15] = {0};
      
      du[0][3] = -3; du[0][6] = +3;
      du[1][4] = -3; du[1][7] = +3;
      du[2][5] = -3; du[2][8] = +3;
      
      dv[0][6] = -3; dv[0][9] = +3;
      dv[1][7] = -3; dv[1][10] = +3;
      dv[2][8] = -3; dv[2][11] = +3;
      
      double dudu[15][15] = {0}, dvdv[15][15] = {0};
      
      loop(j,15) loop(m,3) {
        udu[j] += u[m]*du[m][j];  
        vdv[j] += v[m]*dv[m][j];
      }
      
      loop(j,15) loop(k,15)
      loop(m,3) {
          dudu[j][k] += du[m][j]*du[m][k];
          dvdv[j][k] += dv[m][j]*dv[m][k];
      }
      
      double phi[3] = {0};
      double dphidq[3][15] = {0};
      double d2phidqq[3][15][15] = {0};
      
      loop(i,3) phi[i] = u[i]/squu - v[i]/sqvv;
      
      loop(i,3) loop(j,15) {
        dphidq[i][j] = (du[i][j]*uu - u[i]*udu[j])/(uu*squu) - (dv[i][j]*vv - v[i]*vdv[j])/(vv*sqvv); 
      }
      
      loop(i,3) loop(j,15) loop(k,15) {
        d2phidqq[i][j][k] += ( uu*( 2*du[i][k]*udu[j] - u[i]*dudu[j][k] -du[i][j]*udu[k] ) - 3*udu[j]*( du[i][k]*uu - u[i]*udu[k] ) )/(uu*uu*squu)
                                        - ( vv*( 2*dv[i][k]*vdv[j] - v[i]*dvdv[j][k] -dv[i][j]*vdv[k] ) - 3*vdv[j]*( dv[i][k]*vv - v[i]*vdv[k] ) )/(vv*vv*sqvv);
      }
      
      double vsys0[15] = {0};
      loop(m,5) {
          int loc = FN->C2_cps[c].comp[m]-1;
          loop(j,3) vsys0[3*m+j] = gsl_vector_get(vsys, 3*loc+j);
      }
      
      double Gammai[3] = {0};
      
      loop(j,15) loop(k,15) loop(m,3){
        Gammai[m] += -d2phidqq[m][j][k]*vsys0[j]*vsys0[k];
      }
    
      loop(m,5) {
        int loc = FN->C2_cps[c].comp[m]-1;
      
        loop(j,3) {
          gsl_matrix_set(Lambda, 3*c, 3*loc+j, dphidq[0][3*m+j]);
          gsl_matrix_set(Lambda, 3*c+1, 3*loc+j, dphidq[1][3*m+j]);
          gsl_matrix_set(Lambda, 3*c+2, 3*loc+j, dphidq[2][3*m+j]);
        }
      }
      
      loop(j,3) gsl_matrix_set(Gamma, FN->nfdof + 3*c+j, 0, Gammai[j]);
      
      if(C2==1) {
        
        double Phi[3] = {0};
        double dPhidq[3][15] = {0}, U[3][15] = {0}, V[3][15] = {0};
        double d2Phidqq[3][15][15] = {0}, U2[3][15][15] = {0}, V2[3][15][15] = {0};
      
        double x[3] = { 6*(q[0]-2*q[3]+q[6]), 6*(q[1]-2*q[4]+q[7]), 6*(q[2]-2*q[5]+q[8]) };
        double y[3] = { 6*(q[6]-2*q[9]+q[12]), 6*(q[7]-2*q[10]+q[13]), 6*(q[8]-2*q[11]+q[14]) };
    
        double dx[3][15] = {0};
        double dy[3][15] = {0};
        
        dx[0][0] = +6;  dx[0][3] = -12; dx[0][6] = +6;
        dx[1][1] = +6;  dx[1][4] = -12; dx[1][7] = +6;
        dx[2][2] = +6;  dx[2][5] = -12; dx[2][8] = +6;
        
        dy[0][6] = +6;  dy[0][9] = -12; dy[0][12] = +6;
        dy[1][7] = +6;  dy[1][10] = -12; dy[1][13] = +6;
        dy[2][8] = +6;  dy[2][11] = -12; dy[2][14] = +6;
        
        double ux = u[0]*x[0] + u[1]*x[1] + u[2]*x[2];
        double vy = v[0]*y[0] + v[1]*y[1] + v[2]*y[2];
      
        loop(i,3) Phi[i] = (uu*x[i] - ux*u[i])/(uu*uu) - (vv*y[i] - vy*v[i])/(vv*vv);
        
        double udx[15] = {0}, xdu[15] = {0}, vdy[15] = {0}, ydv[15] = {0} ;

        loop(j,15) loop(m,3) {
          udx[j] += u[m]*dx[m][j];    xdu[j] += x[m]*du[m][j];    
          vdy[j] += v[m]*dy[m][j];    ydv[j] += y[m]*dv[m][j];    
        }
        
        loop(i,3) loop(j,15) {
          U[i][j] = uu*( uu*dx[i][j] - x[i]*2*udu[j] - ux*du[i][j] - u[i]*(udx[j] + xdu[j])  ) + 4*u[i]*ux*udu[j];
          V[i][j] = vv*( vv*dy[i][j] - y[i]*2*vdv[j] - vy*dv[i][j] - v[i]*(vdy[j] + ydv[j])  ) + 4*v[i]*vy*vdv[j];
          
          dPhidq[i][j] = U[i][j]/(uu*uu*uu) - V[i][j]/(vv*vv*vv);
        }
        
         double dxdu[15][15] = {0}, dydv[15][15] = {0};
        
        loop(j,15) loop(k,15)
        loop(m,3) {
          dxdu[j][k] += dx[m][j]*du[m][k];
          dydv[j][k] += dy[m][j]*dv[m][k];
        }
        
        loop(i,3) loop(k,15) loop(j,15) {
          U2[i][k][j] = uu*( dx[i][k]*2*udu[j] - 2*dx[i][j]*udu[k] - 2*x[i]*dudu[j][k]  - du[i][k]*(udx[j] + xdu[j]) 
                                         - u[i]*dxdu[k][j] - du[i][j]*udx[k] - u[i]*dxdu[j][k] - du[i][j]*xdu[k]  )
                                + (  uu*dx[i][k] - x[i]*2*udu[k] - ux*du[i][k] - u[i]*(udx[k] + xdu[k])   )*2*udu[j]
                                + 4*ux*( u[i]*dudu[j][k] + du[i][j]*udu[k] ) + 4*u[i]*udu[k]*(udx[j] + xdu[j]) ;
                        
          V2[i][k][j] = vv*( dy[i][k]*2*vdv[j] - 2*dy[i][j]*vdv[k] - 2*y[i]*dvdv[j][k]  - dv[i][k]*(vdy[j] + ydv[j]) 
                                         - v[i]*dydv[k][j] - dv[i][j]*vdy[k] - v[i]*dydv[j][k] - dv[i][j]*ydv[k]  )
                                + (  vv*dy[i][k] - y[i]*2*vdv[k] - vy*dv[i][k] - v[i]*(vdy[k] + ydv[k])   )*2*vdv[j]
                                + 4*vy*( v[i]*dvdv[j][k] + dv[i][j]*vdv[k] ) + 4*v[i]*vdv[k]*(vdy[j] + ydv[j]) ;
          
          d2Phidqq[i][j][k] = ( (uu*uu*uu)*U2[i][k][j] - U[i][k]*6*(uu*uu)*udu[j] )/ (uu*uu*uu*uu*uu*uu)
                                          - ( (vv*vv*vv)*V2[i][k][j] - V[i][k]*6*(vv*vv)*vdv[j] )/ (vv*vv*vv*vv*vv*vv) ;
        }
        
        double Gammai2[3] = {0};
      
        loop(j,15) loop(k,15) loop(m,3){
          Gammai2[m] += -d2Phidqq[m][j][k]*vsys0[j]*vsys0[k];
        }
    
        loop(m,5) {
          int loc = FN->C2_cps[c].comp[m]-1;
        
          loop(j,3) {
            gsl_matrix_set(Lambda, 3*c+3, 3*loc+j, dPhidq[0][3*m+j]);
            gsl_matrix_set(Lambda, 3*c+4, 3*loc+j, dPhidq[1][3*m+j]);
            gsl_matrix_set(Lambda, 3*c+5, 3*loc+j, dPhidq[2][3*m+j]);
          }
        }
      
        loop(j,3) gsl_matrix_set(Gamma, FN->nfdof + 3*c+3+j, 0, Gammai2[j]);
        
        if(VERBOSE1) printf("Constraint C2 applied: %d %d %d %d %d %d %d %.16f %.16f %.16f %.16f %.16f\n",FN->tstep,c,bi,ci,di,ei,fi,uu,vv, Phi[0], Phi[1], Phi[2]);
      
      }
      
      if(VERBOSE1) printf("Constraint C1 applied: %d %d %d %d %d %d %d %.16f %.16f %.16f %.16f %.16f\n",FN->tstep,c,bi,ci,di,ei,fi,uu,vv, phi[0], phi[1], phi[2]);
  }
  
  if(VERBOSE1) gsl_matrix_print(Gamma,"Gamma:");
    
  //determine and apply constraint forces
  
  gsl_matrix_set_zero(Lsys);
  int count = 0;
  
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    gsl_matrix_set(Lsys, count++, it->first-1, 1.0);
  }

  //printf("\n%d ncon = %d",count,FN->ncon);
  //gsl_matrix_print(Lambda,"Lambda");

  //gsl_matrix_print(Lambda,"Lambda");
  
  if(1) loop(i,3*ncc_active) loop(j,3*np) {
    
    double val = gsl_matrix_get(Lambda,i,j);
    gsl_matrix_set(Lsys, count+i, j, val);
    
  }
   
  gsl_vector *gsys = gsl_vector_alloc (3*np);
  gsl_vector *g2sys = gsl_vector_alloc (3*np);

  gsl_vector *Lg2sys = gsl_vector_alloc (ncon_active);

  loop(j,3*np) {
    double val = 0;
    
    loop(k,3*np) loop(l,3*np) {
      
      int j1 = j/3;
      int k1=k/3;
      int m = l/3;
      int n = l%3;
      
      if( j%3 == k%3 ) val += FN->dMdQ(j1,k1,m,n) * gsl_vector_get(vsys,k) * gsl_vector_get(vsys,l);
    }
    gsl_vector_set(gsys, j, val);
  
  }
  
  gsl_matrix_memcpy(GMtmp2, FN->GM3);
  //gsl_matrix_print(GMtmp2, "GMtmp (284) before LU:");
  
  gsl_linalg_LU_decomp (GMtmp2, p1, &s);
  gsl_linalg_LU_solve (GMtmp2, p1, gsys, g2sys);
  
  //printf("LU code (284): %d success = %d\n", res284, GSL_SUCCESS);
  //gsl_matrix_print(GMtmp, "GMtmp (283):");
  
  gsl_blas_dgemv (CblasNoTrans, 1.0, Lsys, g2sys, 0.0, Lg2sys);

  //printf ("b = \n");
  //gsl_vector_fprintf (stdout, Lg2sys, "%.24f");
  
  gsl_permutation_free(p1);
  gsl_vector_free(wsys);
  gsl_vector_free(vsys);
  gsl_matrix_free(GMtmp1);
  gsl_matrix_free(GMtmp2);
  
  //gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Lsys, FN->b2sys, 0.0, bsys);
  
  //printf("bsys parts: \n");
  loop(j,ncon_active) {
    
    double val = gsl_vector_get(Lg2sys,j) + gsl_matrix_get(Gamma ,j,0);    
    //printf("%lf %lf %lf \n",gsl_vector_get(Lg2sys,j), gsl_matrix_get(Gamma ,j,0), val);
    gsl_matrix_set(bsys,j,0, val );
  
  }
  //gsl_matrix_print(bsys,"bsys:");

  gsl_matrix_solve_LMinvLt(Lsys, FN->GM3, AMA, ncon_active, 3*np);
  
  loop(i,np) loop(j,3) gsl_matrix_set(FN->Hq, i*3+j, 0, gsl_matrix_get(FN->qH,i,j) );

  gsl_matrix *xm = gsl_matrix_alloc(3*np,1);
  gsl_matrix_solve_Ainvb(FN->GM3, xm, FN->Hq, 3*np);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Lsys, xm, 0.0, AHq);

  gsl_matrix_free(xm);

  gsl_matrix_add(bsys, AHq); //b+AHq
  
  xm = gsl_matrix_alloc(ncon_active,1);
  gsl_matrix *qc = gsl_matrix_alloc(3*np,1);
  
  gsl_matrix_solve_Ainvb2(AMA, xm, bsys, ncon_active); //uncomment
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, Lsys, xm, 0.0, qc);
  
  //gsl_matrix_print(FN->qH,"qH (before):");
  //gsl_matrix_print(qc,"qc:");
  
  loop(i,np) loop(j,3) {
    double qc0 =  gsl_matrix_get(FN->qH,i,j) - 1.0*gsl_matrix_get(qc, i*3+j,0); //gsl_matrix_get(FN->Hq,i*3+j,0);
    
    gsl_matrix_set(FN->qH, i, j, qc0);
  }
  
  //gsl_matrix_print(FN->qH,"qH (final):");
  //gsl_matrix_print(FN->pH,"pH (final):");
  
  gsl_matrix_free(xm);
  gsl_matrix_free(qc);
  
  gsl_vector_free(gsys);
  gsl_vector_free(g2sys);
  gsl_vector_free(Lg2sys);
  gsl_matrix_free(Lambda);
  gsl_matrix_free(Gamma);
  
  gsl_matrix_free(Lsys);
  gsl_matrix_free(bsys);
  
  gsl_matrix_free(AMA);
  gsl_matrix_free(AHq);
  
}

void tao_integrate(fibnetwork &fn, int nsteps, double tstep)
{
  size_t np = fn.np;
	
	fn.Psys = gsl_matrix_alloc (np, 3);
	fn.Qsys = gsl_matrix_alloc (np, 3);
	fn.Rsys = gsl_matrix_alloc (np, 3);
	fn.Fsys = gsl_matrix_alloc (np, 3);
	fn.Wsys = gsl_matrix_alloc (np, 3);
	fn.qH = gsl_matrix_alloc (np, 3);
	fn.pH = gsl_matrix_alloc (np, 3);
  
	fn.GM = gsl_matrix_alloc (np, np);
	fn.GMinv = gsl_matrix_alloc (np, np);
	fn.GMdot = gsl_matrix_alloc (np, np);
	fn.GN = gsl_matrix_alloc (np, np);
  
  fn.Lsys = gsl_matrix_alloc (fn.ncon, 3*np);
  fn.Asys = gsl_matrix_alloc (fn.ncon, 3*np);
  fn.bsys = gsl_matrix_alloc (fn.ncon, 1);
  fn.b2sys = gsl_matrix_alloc (3*np, 1);
  
  fn.GM3 = gsl_matrix_alloc (3*np, 3*np);
  fn.GM3inv = gsl_matrix_alloc (3*np, 3*np);
  
  fn.qsys = gsl_matrix_alloc (3*np, 1);
  fn.dqsys = gsl_matrix_alloc (3*np, 1);
  fn.psys = gsl_matrix_alloc (3*np, 1);
  
  fn.Hq = gsl_matrix_alloc (3*np, 1);
  fn.MA = gsl_matrix_alloc (3*np, fn.ncon);
  fn.AMA = gsl_matrix_alloc (fn.ncon, fn.ncon);
  
  fn.MAX = gsl_matrix_alloc (3*np, fn.ncon);
  fn.X = gsl_matrix_alloc (fn.ncon, fn.ncon);
  fn.AHq = gsl_matrix_alloc (fn.ncon, 1);
  
  fn.V = gsl_matrix_alloc (fn.ncon, fn.ncon);
  fn.S = gsl_vector_alloc (fn.ncon);
  
  double *X = new double[np*3];
  double *Y = new double[np*3];
  double *Q = new double[np*3];
  double *P = new double[np*3];
  
	loop(i,np) loop(j,3) {
    
    Q[i*3+j] = fn.CPx[i].comp[j];
    Y[i*3+j] = fn.CPw[i].comp[j];
    
    X[i*3+j] = fn.CPx[i].comp[j];
    P[i*3+j] = fn.CPw[i].comp[j];

  }
	fprintf(fn.flog,"\n"); 
	
  for (int s = 1; s <= nsteps; s++) {
    
    //printf("\nSTART TAO STEP:\n");
    //dH(Q,P,&fn,0);
    //fn.totE = fn.axialE + fn.bendE + fn.cohE + fn.kinE;
    //printf("\nPRE: %d %lf %lf %lf %lf %lf\n", s, fn.axialE, fn.bendE, fn.cohE, fn.kinE, fn.totE);
    
    //TAO STEP 1
    //printf("\nSTEP 1:\n");
    dH(Q,Y,&fn,0);
    loop(i,np) loop(j,3) {
      X[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      P[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }
    
    //TAO STEP 2
    //printf("\nSTEP 2:\n");
    dH(X,P,&fn,0);
    loop(i,np) loop(j,3) {
      Q[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      Y[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }

    //TAO STEP 3
    //printf("STEP 3:\n");
    loop(i,np) loop(j,3) {
      
      double a = Q[i*3+j]+X[i*3+j];
      double b = Q[i*3+j]-X[i*3+j];
      
      double c = P[i*3+j]+Y[i*3+j];
      double d = P[i*3+j]-Y[i*3+j];
      
      Q[i*3+j] = 0.5*( a + b*cos(2*tao_omega*tstep) + d*sin(2*tao_omega*tstep) );
      P[i*3+j] = 0.5*( c + d*cos(2*tao_omega*tstep) - b*sin(2*tao_omega*tstep) );
      X[i*3+j] = 0.5*( a - b*cos(2*tao_omega*tstep) - d*sin(2*tao_omega*tstep) );
      Y[i*3+j] = 0.5*( c - d*cos(2*tao_omega*tstep) + b*sin(2*tao_omega*tstep) );
    }
    
    //TAO STEP 4
    //printf("\nSTEP 4:\n");
    dH(X,P,&fn,0);
    loop(i,np) loop(j,3) {
      Q[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      Y[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }
    
    //TAO STEP 5
    //printf("\nSTEP 5:\n");

    dH(Q,Y,&fn,0);
    loop(i,np) loop(j,3) {
      X[i*3+j] += 0.5*tstep*gsl_matrix_get(fn.pH,i,j);
      P[i*3+j] -= 0.5*tstep*gsl_matrix_get(fn.qH,i,j);
    }
    
    dH(Q,P,&fn,0);
    
    fn.tstep = s;
    
    fn.printlammps_cps((char*) DUMP1,(char*) "a");
    if(s%1==0) fn.printlammps((char*) DUMP2,(char*) "a",NBEZ);
    
    fn.totE = fn.axialE + fn.bendE + fn.cohE + fn.kinE;
    fprintf(fn.flog,"Tao Integration: %6d %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n", s, fn.axialE, fn.bendE, fn.cohE, fn.kinE, fn.totE, fn.momT[0], fn.momT[1], fn.momT[2]);
  }
  
  delete[] X;
  delete[] Y;
  delete[] Q;
  delete[] P;
  
  gsl_matrix_free(fn.Psys);
  gsl_matrix_free(fn.Qsys);
  gsl_matrix_free(fn.Rsys);
  gsl_matrix_free(fn.Fsys);
  gsl_matrix_free(fn.Wsys);
  gsl_matrix_free(fn.qH);
  gsl_matrix_free(fn.pH);
  
  gsl_matrix_free(fn.GM);
  gsl_matrix_free(fn.GMinv);
  gsl_matrix_free(fn.GMdot);
  gsl_matrix_free(fn.GN);
  
  gsl_matrix_free(fn.Lsys);
  gsl_matrix_free(fn.Asys);
  gsl_matrix_free(fn.bsys);
  gsl_matrix_free(fn.b2sys);
  
  gsl_matrix_free(fn.GM3);
  gsl_matrix_free(fn.GM3inv);
  
  gsl_matrix_free(fn.qsys);
  gsl_matrix_free(fn.dqsys);
  gsl_matrix_free(fn.psys);
  
  gsl_matrix_free(fn.Hq);
  gsl_matrix_free(fn.MA);
  gsl_matrix_free(fn.AMA);
  
  gsl_matrix_free(fn.MAX);
  gsl_matrix_free(fn.X);
  gsl_matrix_free(fn.AHq);
  
  gsl_matrix_free(fn.V);
  gsl_vector_free(fn.S);
  
}
  
