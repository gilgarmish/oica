
/*
 * USAGE. The function is called from Matlab as follows
 *
 *      [out1 out2] = jd_in( B(:), m, n, eps, V0(:));
 *
 * Input:
 * 	 B     : is m-by-(m*n) matrix such that B = [B1 B2 ... Bn];
 * 	 B(:)  : matrix B as a vector of length m*m*n (see Matlab's docs)
 * 	 m     : size of matrices
 *   n     : number of matrices
 *	 eps   : a threshold for the stopping criterion (default 1e-8)
 *	 V0    : starting point - any orthogonal matrix (default identity)
 *	 V0(:) : matrix V0 as a vector of length m*m
 *
 * ALL inputs are compulsary!
 *
 * Output:
 *   out1 : optimal matrix V as a vector of length m*m
 *   out2 : matrix A = [A1 A2 ... An] of approximately diagonal matrices
 *             as a vector of length m*m*n
 *
 * Examples: "algorithms/joint_diagonalization.m"
 *           "example_joint_diagonalization.m"
 */

/* DESCRIPTION. Given n matrices B1, B2, ..., Bn, each of size m-by-m,
 * the function looks for such orthogonal matrix V that all matrices
 *          A1 = V'*B1*V, A2 = V'*B2*V, ..., An = V'*Bn*V 
 * are (approximately) as diagonal as possible. For that, it minimizes
 * the sum of all squared off-diagonal elements, i.e. the objective
 *      obj(V) = \sum_{r=1^n} \sum_{1<=i \ne j<=m} ( Ar[i,j] )^2.
 * The minimization procedure is based on the itertive Jacobi (=Givens)
 * rotations. At each iteration, a closed form solution is computed for 
 * the respective Jacobi rotation, making the algorithm fast.
 */

/*
 * ACKNOWLEDGMENTS. This is C++/MEX-Matlab implementation of the joint 
 * diagonalization algorithm (for real matrices) as described in the paper:
 *   J.F. Cardoso, A. Souloumiac. Jacobi angles for simultaneous 
 *   diagonalization. SIAM Jounal on Mathematical Analysis, 1996.
 * 
 * The original Matlab and R implementations can be found here:
 *   http://perso.telecom-paristech.fr/~cardoso/jointdiag.html
 */

/* Copyright 2015, Anastasia Podosinnikova */

#include <math.h>
#include "mex.h"
#include "matrix.h"


void ojd_routine(double* A, int m, int n, double eps, double* V)
{
  
  double g11, g12, g21, g22, g1i, g2i;
  double ton, toff, theta, c, s;
  double tempp, tempq;
  
  int p, q, i, t, r;
  
  bool encore = true;
  while (encore)
  {
    encore = false;
    for (p = 0; p < m - 1; p++)
    {
      for (q = p + 1; q < m; q++)
      {
        /* computation of the Jacobi (=Givens) rotations */
        g11 = 0; g22 = 0; g12 = 0; g21 = 0;
        for (i = 0; i < n; i++)
        {
          g1i = A[i*m*m+p*m+p] - A[i*m*m+q*m+q];
          g2i = A[i*m*m+q*m+p] + A[i*m*m+p*m+q];
          g11 += g1i*g1i; g21 += g2i*g1i;
          g12 += g1i*g2i; g22 += g2i*g2i;
        }
        ton = g11 - g22; toff = g12 + g21;
        theta = 0.5 * atan2(toff, ton + sqrt(ton*ton+toff*toff));
        c = cos(theta); s = sin(theta);
        if (encore || (fabs(s)>eps)) {encore=true;}
   
        /* update of A and V */
        if (fabs(s) > eps)
        {
          for (t = 0; t < m; t++)
          {
            for (i = 0; i < n; i++)
            {
              tempp = A[i*m*m+p*m+t];
              tempq = A[i*m*m+q*m+t];
              A[i*m*m+p*m+t] = c*tempp+s*tempq;
              A[i*m*m+q*m+t] = c*tempq-s*tempp;
            }
          }
          
          for (i = 0; i < n; i++)
          {
            for (r = 0; r < m; r++)
            {
              tempp=A[i*m*m+r*m+p];
              tempq=A[i*m*m+r*m+q];
              A[i*m*m+r*m+p]=c*tempp+s*tempq;
              A[i*m*m+r*m+q]=c*tempq-s*tempp;
            }
          }
          
          for (t = 0; t < m; t++)
          {
            tempp = V[p*m+t]; 
            tempq = V[q*m+t];
            V[p*m+t] = c*tempp + s*tempq;
            V[q*m+t] = c*tempq - s*tempp;
          }
          
        } /* end of update */
        
      }
    }
  }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
  /* checking if the number of inputs/outputs is correct */
  if (nrhs!=7) {mexErrMsgTxt("Wrong number of inputs.");}
  if (nlhs>2)  {mexErrMsgTxt("Wrong number of outputs.");}
  
  /* matrix B as a m*m*n-vector, i.e. B(:) */
  double* B  = mxGetPr(prhs[0]);
  /* m is the size of each matrix (different from the paper's notation) */
  int m = (int)mxGetScalar(prhs[1]);
  /* n is the number of matrices (different from the paper's notation) */
  int n = (int)mxGetScalar(prhs[2]);
  /* checking if the size of B is correct */
  if ( ((int)mxGetM(prhs[0])!=m*m*n) || ((int)mxGetN(prhs[0])!=1) )
  {
    mexErrMsgTxt("Wrong size of input B. Did you use B(:)?");
  }
  /* the value of the treshold used for the stopping criterion */
  double eps = (double)mxGetScalar(prhs[3]);
  /* matrix V0 as a m*mvector, i.e. V0(:) */
  double* V0 = mxGetPr(prhs[4]);
  /* checking if the size of V0 is correct */
  if ( ((int)mxGetM(prhs[4])!=m*m) || ((int)mxGetN(prhs[4])!=1) )
  {
    mexErrMsgTxt("Wrong size of input V0. Did you use V0(:)?");
  }
  int kmax = (int)mxGetScalar(prhs[5]);
  int isdebug = (int)mxGetScalar(prhs[6]);
  
  /* initializing output V with V0 */
  plhs[0] = mxCreateDoubleMatrix(m*m, 1, mxREAL);
  double* V = mxGetPr(plhs[0]);
  for (int i = 0; i < m*m; i++) {V[i] = V0[i];}
  /* initializing output A with B (different from the paper's notation) */
  plhs[1] = mxCreateDoubleMatrix(m*m*n, 1, mxREAL);
  double* A = mxGetPr(plhs[1]);
  for (int i = 0; i < m*m*n; i++) {A[i] = B[i];}
  /* WARNING: congruence transform of A with V is mandatory in the caller
   *          function, unless V is the identity */
  
  ojd_routine(A, m, n, eps, V);
  
}

