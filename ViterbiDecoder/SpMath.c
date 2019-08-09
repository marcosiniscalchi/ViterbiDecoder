/*-------------------------------------------------------*/
/*                       SpMath.c                        */
/*             Module: Related Math Operations           */
/*                                                       */
/* (1) Naive Matrix Inversion                            */
/* (2) Naive Gauss Elimination                           */
/* (3) Naive Solving Systems of Linear Equations         */
/* (4) Calculate likelihood/log-likelihood               */
/* (5) Calculate divergence between 2 Gaussians          */
/* (6) Calculate Bhattacharyya distance between          */
/*     2 Gaussians                                       */
/*                                                       */
/*                             Version: 0.10             */
/*                             Kenneth, Kuan-Ting Chen   */
/* Bugs Report: please contact kenneth@iis.sinica.edu.tw */
/*-------------------------------------------------------*/
#include "BasicDef_MAP.h"



/*--------------------------- Take Log-value -----------------------------*/
float TakeLog(float value)
{
	float log_value;
	
	if (value<=0.0)
		log_value = MinProb;
	else
		log_value = log(value);

	return(log_value);
}

/*--------------------------- Gauss Elimination --------------------------*/
void GaussEliminate(Matrix A, int dim, int *perm, int *sign) {
   Vector s;    /* scale vector */
   int i, j, k, imax, li, lk;
   float r, rmax, smax, xmult, temp;

   *sign = 1;
   s = CreateVector(dim);
   /* initialize perm[], s[] */
   for (i=0; i<dim; i++) {
      perm[i] = i;
      smax = 0.0;
      for (j=0; j<dim; j++) {
	 if ((temp = fabs(A[i][j])) > smax) {
	    smax = temp;
	 }
      }
      if (smax == 0.0) {
	 printf("ERROR! Matrix is Singular!\n");
	 exit (-1);
      }
      s[i] = smax;
      /*printf("%d %f\n", perm[i], s[i]);*/
   } /* end initialization */

   for (k=0; k<dim-1; k++) {  /* index of Gauss Elimination procedure */
      rmax = 0.0;
      for (i=k; i<dim; i++) { /* choose pivot */
	 /*printf("start choosing pivot\n");*/
	 li = perm[i]; /*printf("li=%d", li);*/
	 r = fabs(A[li][k]/s[li]);
	 if (r > rmax) {
	    rmax = r;
	    j = i;    /* j indicates the pivot row */
	 }
	 if (rmax == 0.0) {
	    printf("\nERROR! rmax=0.0,Matrix is Singular!\n");
	    exit (-1);
	 }
      } /* end choosing pivot */
      imax = perm[j];
      perm[j] = perm[k];    /* interchange perm[j] & perm[k] */
      perm[k] = imax;
      if (((j-k)%2) == 1)
	 *sign = -(*sign);

      for (i=k+1; i<dim; i++) { /* for row i below row k */
	 lk = perm[k];
	 li = perm[i];
	 xmult = A[li][k]/A[lk][k];
	 A[li][k] = xmult;   /* store xmult in A[li][k] */
	 for (j=k+1; j<dim; j++) {  /* for columns in the right-side of k */
	    A[li][j] = A[li][j] - xmult*A[lk][j];
	 }
      } /* end elimination for row i */
   } /* end Gauss elimination */
} /* end of Gauss_eliminate */

/*------------------------- Solve Linear Equation ------------------------*/
void LinearSolve(Matrix A, int dim, int *perm, float *b, float *x) {
   int i, j, k, li, lk, ln;
   float sum;

   /* forward elimination of b[] */
   for (k=0; k<dim-1; k++) {
      lk = perm[k];
      for (i=k+1; i<dim; i++) { /* for the rows below k */
	 li = perm[i];
	 b[li] = b[li] - A[li][k]*b[lk];
      }
   }

   /* backward substraction */
   ln = perm[dim-1];
   x[dim-1] = b[ln]/A[ln][dim-1];

   for (i=dim-2; i>=0; i--) {
      li = perm[i];
      sum = b[li];
      for (j=i+1; j<dim; j++) {
	 sum = sum - A[li][j]*x[j];
      }
      x[i] = sum/A[li][i];
   }
}

/*---------------------------- Invert Matrix -----------------------------*/
int Invert(Matrix C, Matrix C_inv, int dim) {
   Matrix A;
   float col[100], x[100], xx;
   float log_det;
   int sign;
   int i, j, perm[100];
   int li;

   A = CreateMatrix(dim, dim);
   for (i=0; i<dim; i++)    /* make a copy of C */
      for (j=0; j<dim; j++)
	 A[i][j] = C[i][j];

   GaussEliminate(A, dim, perm, &sign);
   for (j=0; j<dim; j++) {
      for (i=0; i<dim; i++) { /* reset col[] */
	 col[i] = 0.0;
      }
      col[j] = 1.0;
      LinearSolve(A, dim, perm, col, x);   /* A*col'=col, solve column by column*/
      for (i=0; i<dim; i++) {
	 C_inv[i][j] = x[i];
      }
   }

   log_det = 0.0;  /* calculate log(det(C)) */
   
   for (i=0; i<dim; i++) {
      li = perm[i];
      xx = A[li][i];
      if (fabs(xx) == 0.0) {
         return 0;
      }
      if (xx < 0.0)
	 sign = -sign;
      log_det += log(fabs(xx));
   }
   
   FreeMatrix(A, dim);
   return 1;
}

/*-------------------------- Calculate bj(ot) ----------------------------*/
float CalStateLikeli(int LogKey, int M, Vector data, Matrix Mean, Matrix Var, Vector C, Vector GConst)
{		            /* mixture no */        /* mean[m][k] */            /* C[m], GConst[m] */	
	int m, i;
   	float tmp1, tmp2, mix_likeli;
   	float sum = 0.0;
   	float state_likeli;

   	for (m=0; m<M; m++) 
	{
      		tmp1 = (-1.0)*GConst[m];
      		for (i=0; i<FeaDim; i++) 
		{
	 		tmp2 = data[i] - Mean[m][i];   /* o-mu */
	 		tmp1 = tmp1 - (tmp2*tmp2)*Var[m][i];
      		}
      		mix_likeli = exp((0.5)*tmp1);
      		if (mix_likeli <= 0.0)
	 		mix_likeli = 1e-30;
      		sum += C[m]*mix_likeli;
   	} /* end loop for mixture component */
   	state_likeli = log(sum);
   	if (LogKey == 0)
      		return (sum);
   	else
      		return (state_likeli);
}

/*-------------------------- Calculate bjk(ot) ---------------------------*/
float CalMixLikeli(int LogKey, Vector data, Vector Mean, Vector Var, float GConst)
{
   	int i;
   	float tmp1, tmp2;
   	float sum = 0.0;
   	float mix_likeli;

   	tmp1 = (-1.0) * GConst;
   	for (i=0; i<FeaDim; i++) 
	{
      		tmp2 = data[i] - Mean[i];
      		tmp1 = tmp1 - (tmp2*tmp2)*Var[i];
   	}
   	sum = exp((0.5)*tmp1);
   	if (sum <= 0.0)
      		sum = 1e-30;
   	mix_likeli = log(sum);
   	if (LogKey == 0)
      		return (sum);
   	else
      		return (mix_likeli);
}

float CalDivergence(Vector mean1, Vector sigma1, Vector mean2, Vector sigma2, int dim)
{
   int i, j, k;
   float divergence = 0.0;
   float diff = 0.0;
   float delta_kmi2 = 0.0;
   float term = 0.0;

   for (i=0; i<dim; i++) {
      diff = mean1[i] - mean2[i];
      delta_kmi2 = diff * diff;
      term = ((sigma1[i]+delta_kmi2)/sigma2[i])+((sigma2[i]+delta_kmi2)/sigma1[i]);
      divergence += term;
   }

   return (divergence);
}

/* kenneth_990524: modify to handle stream-based case */
float CalBhattaDistance(Vector mean1, Vector sigma1, Vector mean2, Vector sigma2, int begin, int end)
{
   int i, j, k;
   float distance = 0.0;
   float diff = 0.0;
   float log_det_sigma1 = 0.0;
   float log_det_sigma2 = 0.0;
   float log_det_sigma1plus2 = 0.0;
   float term1 = 0.0;
   float term2 = 0.0;

   for (i=begin; i<=end; i++) {
      diff = mean1[i] - mean2[i];
      term1 += diff * diff * (1.0)/(0.5*(sigma1[i]+sigma2[i]));
   }
   term1 *= (1.0)/(8.0);

   for (i=begin; i<=end; i++) {
      log_det_sigma1 += log(sigma1[i]);
      log_det_sigma2 += log(sigma2[i]);
      log_det_sigma1plus2 += log(((sigma1[i]+sigma2[i])/2.0));
   }
   term2 = (0.5) * (log_det_sigma1plus2 - (0.5)*log_det_sigma1 - (0.5)*log_det_sigma2);

   distance = term1 + term2;
   return (distance);
}