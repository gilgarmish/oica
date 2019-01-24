/* Copyright 2019, Anastasia Podosinnikova */

#include <math.h>
#include "mex.h"
#include "matrix.h"

void quadricov_routine
        (double** data, int n, int p, double** Q)
{
  int i,a,b,c,d;
  
  double** C; /* biased covariance matrix */
  C = new double*[p];
  for (a=0; a<p; a++)
    C[a] = new double[p];
  for (a=0; a<p; a++){
    for (b=0; b<p; b++)
      C[a][b] = 0;
  }
  
  double tt;
  double* temp;
  temp = new double[p*p];
  
  for (i=0; i<n; i++){
    
    for (a=0; a<p; a++){
      for (b=0; b<p; b++){
        tt = data[i][a]*data[i][b];
        C[a][b] += tt;
        temp[a+b*p] = tt;
      }
    }
    
    for (a=0; a<p*p; a++){
      for (b=a; b<p*p; b++)
        Q[a][b] += temp[a]*temp[b];
    }
    
  }
  
  for (a=0; a<p; a++){
    for (b=0; b<p; b++)
      C[a][b] = C[a][b]/n;
  }
  for (a=0; a<p*p; a++){
    for (b=a; b<p*p; b++)
      Q[a][b] = Q[a][b]/n;
  }
  
  for (a=0; a<p; a++){
    for (b=0; b<p; b++){
      for (c=0; c<p; c++){
        for (d=0; d<p; d++){
          if ( c+d*p >= a+b*p ){
            Q[a+b*p][c+d*p] = Q[a+b*p][c+d*p] - C[a][b]*C[c][d] - C[a][c]*C[b][d] - C[a][d]*C[b][c];
          }
        }
      }
    }
  }
  
}

/* Input:  data matrix of size p-by-n */
/*         p is the dimension */
/*         n is the sample size */
/*         NOTE: data matrix must be centered!!! */
/* Output: Q is the matrix of the flattened fourth-order cumulant (quadricovariance) */
/*         Q is of size p*p-by-p*p */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   /* Check the number of inputs and outputs */
  if (nrhs!=1) mexErrMsgTxt("Wrong number of inputs.");
  if (nlhs>1)  mexErrMsgTxt("Wrong number of outputs.");
  
  
  int p = mxGetM(prhs[0]); /* number of rows -- p is the dimension */
  int n = mxGetN(prhs[0]); /* number of cols -- n is the sample size */
  
  /* Get pointer to the data matrix of size p-by-n */
  double* data_pr = (double*)mxGetPr(prhs[0]);
  /* NOTE: The same matrix in C++ is of size n-by-p !!! */
  double** data = new double*[n];
  for (int i=0; i<n; i++)
    data[i] = data_pr + i*p;
  
  
  /* Initialize the output matrix Q of size p*p*p*p-by-s */
  plhs[0] = mxCreateNumericMatrix(p*p, p*p, mxDOUBLE_CLASS, mxREAL);
  double* q_pr = (double*)mxGetData(plhs[0]);
  /* NOTE: The same matrix in C++ is of size p*p-by-p*p */
  /*       However, nothing changes since this matrix is symmetric! */
  double** Q = new double*[p*p];
  for (int j=0; j<p*p; j++)
    Q[j] = q_pr + j*p*p;
  
  for (int i=0; i<p*p; i++){
    for (int j=0; j<p*p; j++)
      Q[i][j] = 0;
  }
  
  quadricov_routine(data,n,p,Q);

}


