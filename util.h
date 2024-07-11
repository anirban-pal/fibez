#define loop(x,n) for(int x = 0; x < n; ++x)
#define loop2(x,a,b) for(int x = a; x < b; ++x)

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

typedef std::stringstream sss;

static int singular (const gsl_matrix * LU)
{
  size_t i, n = LU->size1;

  for (i = 0; i < n; i++)
    {
      double u = gsl_matrix_get (LU, i, i);
      if (u == 0) return 1;
    }
 
  return 0;
}

struct Vector3 {
  double comp[3];
  
  Vector3(double x, double y, double z)
  {
    comp[0] = x;  comp[1] = y;  comp[2] = z;
  }
  
  Vector3()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;
  }
};

struct iVector3 {
  int comp[3];
  
  iVector3(int x, int y, int z)
  {
    comp[0] = x;  comp[1] = y;  comp[2] = z;
  }
  
  iVector3()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;
  }
};

struct Vector4 {
  double comp[4];
  
  Vector4(double t, double x, double y, double z)
  {
    comp[0] = t; comp[1] = x;  comp[2] = y;  comp[3] = z;
  }
  
  Vector4()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;  comp[3] = 0;
  }
};

struct iVector5 {
  int comp[5];
  
  iVector5(int a, int b, int c, int d, int e)
  {
    comp[0] = a;  comp[1] = b;  comp[2] = c;  comp[3] = d;  comp[4] = e;
  }
  
  iVector5()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;  comp[3] = 0;  comp[4] = 0;
  }
};

int gsl_matrix_print(const gsl_matrix *m, const char *title)
{
  int status, n = 0;
  printf("%s\n",title);
  
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      if ((status = printf("%.24f ", gsl_matrix_get(m, i, j))) < 0)
        return -1;
      n += status;
    }
  
    if ((status = printf("\n")) < 0)
      return -1;
    n += status;
  }
  
  return n;
}

double cross_mag(Vector3 u, Vector3 v)
{
  double w0 = u.comp[1]*v.comp[2] - u.comp[2]*v.comp[1];
  double w1 = u.comp[2]*v.comp[0] - u.comp[0]*v.comp[2];
  double w2 = u.comp[0]*v.comp[1] - u.comp[1]*v.comp[0];
  
  return sqrt(w0*w0 + w1*w1 + w2*w2);
}

Vector3 cross(Vector3 u, Vector3 v)
{
  double w0 = u.comp[1]*v.comp[2] - u.comp[2]*v.comp[1];
  double w1 = u.comp[2]*v.comp[0] - u.comp[0]*v.comp[2];
  double w2 = u.comp[0]*v.comp[1] - u.comp[1]*v.comp[0];
  
  Vector3 cprod(w0,w1,w2);
  
  return cprod;
}

Vector3 add(Vector3 u, Vector3 v)
{
  double w0 = u.comp[0] + v.comp[0];
  double w1 = u.comp[1] + v.comp[1];
  double w2 = u.comp[2] + v.comp[2];
  
  Vector3 sum(w0,w1,w2);
  
  return sum;
}

double dot(Vector3 u, Vector3 v)
{
  double w0 = u.comp[0] * v.comp[0];
  double w1 = u.comp[1] * v.comp[1];
  double w2 = u.comp[2] * v.comp[2];
  
  return (w0+w1+w2);
}

void gsl_matrix_solve_Ainvb(const gsl_matrix *A, gsl_matrix *xm, const gsl_matrix *bm, int m)
{
    gsl_vector *x = gsl_vector_alloc (m);
    gsl_vector *y = gsl_vector_alloc (m);
    gsl_matrix *Atmp = gsl_matrix_alloc(m,m);
    gsl_vector *b = gsl_vector_alloc (m);
  
    gsl_matrix_memcpy(Atmp, A);
    
    loop(i,m) gsl_vector_set(b, i, gsl_matrix_get(bm, i, 0) );
  
    gsl_permutation * p = gsl_permutation_alloc (m); int s;
    gsl_linalg_LU_decomp (Atmp, p, &s);
    
    //gsl_matrix_print(Atmp, "Atmp before LU:");
    gsl_linalg_LU_solve (Atmp, p, b, x);
    
    
    
    loop(i,m) gsl_matrix_set(xm, i, 0, gsl_vector_get(x, i) );
  
    
    gsl_vector_free(x);
    gsl_vector_free(y); 
    gsl_vector_free(b);
    gsl_matrix_free(Atmp);
    gsl_permutation_free(p);
}

void gsl_matrix_solve_Ainvb2(const gsl_matrix *A, gsl_matrix *xm, const gsl_matrix *bm, int m)
{
    gsl_vector *x = gsl_vector_alloc (m);
    gsl_vector *y = gsl_vector_alloc (m);
    gsl_matrix *U = gsl_matrix_alloc(m,m);
    gsl_matrix *V = gsl_matrix_alloc(m,m);
    gsl_matrix *Si = gsl_matrix_alloc(m,m);   
    
    gsl_matrix *VSi = gsl_matrix_alloc(m,m); 
    gsl_matrix *VSiU = gsl_matrix_alloc(m,m); 

    gsl_vector *S = gsl_vector_alloc (m);
    gsl_vector *work = gsl_vector_alloc (m);
  
    gsl_matrix_memcpy(U, A);
    
    gsl_vector *b = gsl_vector_alloc (m);
    loop(i,m) gsl_vector_set(b, i, gsl_matrix_get(bm, i, 0) );
  
    //compute SVD of Atmp
    
    gsl_linalg_SV_decomp(U, V, S, work);
    
    gsl_matrix_set_zero(Si);
    loop(i,m) {
      double val = gsl_vector_get(S, i);
      double vali;
      if(fabs(val) < SVD_factor*gsl_vector_get(S, 0)) { 
        vali = 0;
      } else {
        vali = 1.0/val;
      }
      
      gsl_matrix_set(Si, i, i, vali);
    }
    
    //gsl_matrix_print(Si,"Si:");
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, V, Si, 0.0, VSi);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, VSi, U, 0.0, VSiU);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, VSiU, bm, 0.0, xm);

    //gsl_linalg_SV_solve(U, V, S, b, x);
    //loop(i,m) gsl_matrix_set(xm, i, 0, gsl_vector_get(x, i) );
  
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    
    gsl_matrix_free(VSi);
    gsl_matrix_free(VSiU);
    
    gsl_vector_free(x);
    gsl_vector_free(y); 
    gsl_vector_free(b);
    //gsl_permutation_free(p);
}

void gsl_matrix_solve_LMinvLt(const gsl_matrix *L, const gsl_matrix *M, gsl_matrix *z, int m, int n)
{
    gsl_matrix *x = gsl_matrix_alloc (n,m);
    gsl_vector *xi = gsl_vector_alloc (n);
    
    gsl_vector *y = gsl_vector_alloc (n);
    gsl_matrix *Mtmp = gsl_matrix_alloc(n,n);
    gsl_matrix *Mtmp2 = gsl_matrix_alloc(n,n);

    gsl_vector *b = gsl_vector_alloc (n);
  
    gsl_matrix_memcpy(Mtmp, M);
    
    gsl_permutation * p = gsl_permutation_alloc (n); int s;
    gsl_linalg_LU_decomp (Mtmp, p, &s);
    
    
    if(1) loop(i,m) {
      
      gsl_matrix_memcpy(Mtmp2, Mtmp);
      loop(j,n) gsl_vector_set(b, j, gsl_matrix_get(L, i, j) ); //?
      
      gsl_linalg_LU_solve (Mtmp2, p, b, xi);
    
      loop(j,n) gsl_matrix_set(x, j, i, gsl_vector_get(xi, j) );
  
    }
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, L, x, 0.0, z);
    
  
    gsl_vector_free(xi);
    gsl_vector_free(y); 
    gsl_vector_free(b);
    
    gsl_matrix_free(Mtmp);
    gsl_matrix_free(Mtmp2);
    gsl_matrix_free(x);
    
    gsl_permutation_free(p);
}
