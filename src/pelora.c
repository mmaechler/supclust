
    /*******************************************************************/
    /*                                                                 */
    /*  Supervised Variable Clustering with Logistic Ridge Regression  */
    /*  -------------------------------------------------------------  */
    /*                                                                 */
    /*                   Written by Marcel Dettling                    */
    /*                                                                 */
    /*******************************************************************/

/* --- M.Maechler: --
 *
 o	Using LAPACK
 o	exp() and product instead of R_pow()
 o	much clean up
 */

/*  Weitere Hilfsfunktionen aus R */
#include <Rmath.h>

/* BLAS (incl. BLAS 3) -- */
#include <R_ext/BLAS.h>
/* Lapack : */
#include <R_ext/Lapack.h>

#include "supclust.h"


/* NOTA BENE: do  NOT use global variables !
 * --------- */

/* Sum of a vector of length n */
static double sum(double *x, int n)
{
  int i;
  double total=0;
  for(i=0; i < n; i++)
      total += x[i];
  return(total);
}

/*  Diese Funktion berechnet das arith. Mittel eines Vektors der Länge n */
static double mean(double *x, int n)
{
  return(sum(x,n)/n);
}

/*  Auxiliary function: computes the variance of a vector */
static double var(double *x, int n)
{
  /*  Initialization of local variables */
  int    i;
  double mu       = mean(x, n);
  double variance = 0.0;

  /*  Computation of the variance */
  for (i=0; i < n; i++)
    variance += (x[i]-mu)*(x[i]-mu);

  /*  Returning the variance */
  return variance/(n-1);
}


/*  Diese Funktion berechnet die Matrixmultiplikation */
static
void txwx (double *W, double *X, double *tXWX, double *WX, int m, int n)
{
  /*  Vektorendefinitionen */
  static double one  = 1.;
  static double zero = 0.;
  char  *trW = "N", *trX  = "N", *trtX = "T";

  /*  Matrixmultiplikation, Ergebnis wird in C gespeichert */
  F77_CALL(dgemm)(trW, trX, &m, &n, &m, &one, W, &m, X, &m, &zero, WX, &m);
  F77_CALL(dgemm)(trtX,trW, &n, &n, &m, &one, X, &m, WX,&m, &zero, tXWX,&n);
}

/*  Lösen eines Gleichungssystems */
static
void loes_gls (double *A, double *b, int n, int *pivot, int info)
{
  static int one = 1;
  /* Lapack : */
  F77_CALL(dgesv)(&n, &one, A, &n, pivot, b, &n, &info);
}


/*  This function computes several goodness of fit measures */
static
double loss(double *prob, double *y, double *yminp, int n, int crit)
{

  /*  Initialisation of local variables */
  static int    nrv  =  1; /* for daxpy(): Not in ReVerse order */
  static double one  =  1.;
  static double min1 = -1.;
  double gof;
  int     i;

  /*  Computing the goodness of fit measure */
  if (crit==1) {
      F77_CALL(dcopy)(&n, prob, &nrv, yminp, &nrv);
      F77_CALL(dscal)(&n, &min1, yminp, &nrv);
      F77_CALL(daxpy)(&n, &one, y, &nrv, yminp, &nrv);
      gof = F77_CALL(dnrm2)(&n, yminp, &nrv);
      gof *= gof;/* i.e.  = dnrm2() ^ 2 */
  }
  else {
      gof = 0.;
      for (i=0; i < n; i++)
	  gof += -y[i]*log(prob[i]) - (1-y[i])*log1p(-prob[i]);
  }
  /*  Returning the goodness of fit measure */
  return(gof);
}


/*  Compute the parameter vector of logit ridge via Newton-Raphson iterations */
void ridgecoef(double *X, double *W, double *P, double *WX, double *XWX,
	       double *A, double *y, double *prob, double *theta, double *aux,
	       double *yminp, int *pivot,
	       int n, int p, int ap, int info)
{
  /*  Meaning of the input variables */

  /*  X       (n x p) Design matrix:  intercept & cluster representatives */
  /*  W       (n x n) Weight matrix:  contains initial weights */
  /*  P       (p x p) Penalty matrix: (non)-std ridge penalty, var dependent */
  /*  WX      (n x p) Auxiliary matrix: empty */
  /*  XWX     (p x p) Auxiliary matrix: empty */
  /*  A       (p x p) Auxiliary matrix: has the same content as P */
  /*  y       (n x 1) Binary response vector */
  /*  prob    (n x 1) Conditional probability vector: contains initial prob's */
  /*  theta   (p x 1) Parameter vector: contains initial values */
  /*  aux     (p x 1) Auxiliary vector: empty */
  /*  yminp   (n x 1) Auxiliary vector: contains -prob */
  /*  pivot   (p x 1) Pivot vector: empty */
  /*  n       (1)     The number of observations */
  /*  p       (1)     The final number of parameters in the model */
  /*  ap      (1)     The actual number of parameters in the model */
  /*  info    (1)     Auxiliary variable for solving the linear equations */

  int i,j;
  /*  Initialisation of variables */
  int     k       =  0;
  int     psq     =  p*p;
  int     nrv  	  =  1; /* for daxpy(): Not in ReVerse order */
  double  one   =  1.0;
  double  zero  =  0.0;
  double  min1  = -1.0;
  char   *notrans = "N";       /* for dgemv(): Matrix is not transposed */
  char   *trans   = "T";       /* for dgemv(): Matrix is transposed */

  /*  Computation of trans(X)WX */
  txwx(W, X, XWX, WX, n, p);

  /*  Addition of matrices: A <- XWX + P */
  F77_CALL(daxpy)(&psq, &one, XWX, &nrv, A, &nrv);

  /*  Subtraction of vectors: yminp <- y - prob */
  F77_CALL(daxpy)(&n, &one, y, &nrv, yminp, &nrv);

  /*  Computation of: aux <- XWX*theta */
  F77_CALL(dgemv)(notrans,&p,&p,&one, XWX, &p, theta, &nrv, &zero, aux, &nrv);

  /*  Computation of: aux <- trans(X)*yminp + aux */
  F77_CALL(dgemv)(trans, &n, &p, &one, X, &n, yminp, &nrv, &one, aux, &nrv);

  /*  Elimination of zeroes in matrix A */
  if (ap < p) {
      k = 0;
      for (i=0; i < ap; i++) {
	  for (j=0; j < ap; j++) {
	      A[k] = A[(i*p)+j];
	      k++;
	  }
      }
  }

  /*  Solving the linear equation system: aux <- inv(A)*theta */
  loes_gls(A, aux, ap, pivot, info);

  /*  The solution is theta, assign the solution to it: theta <- aux */
  F77_CALL(dcopy)(&p, aux, &nrv, theta, &nrv);

  /*  Compute the linear predictor and store it in yminp: yminp <- X*theta */
  F77_CALL(dgemv)(notrans,&n,&p,&one, X, &n, theta, &nrv, &zero, yminp, &nrv);

  /*  Update prob & W: prob <- 1/(1+exp(-X*theta)), W <- diag(prob*(1-prob)) */
  for (i=0; i < n; i++) {
      prob[i]    = 1/(1+ exp(-yminp[i]));
      W[(i*n)+i] = prob[i] * (1-prob[i]);
  }

  /*  Computation of trans(X)WX */
  txwx(W, X, XWX, WX, n, p);

  /*  Copy the contents of P to A: A <- P */
  F77_CALL(dcopy)(&psq, P, &nrv, A, &nrv);

  /*  Addition of matrices: A <- XWX + P */
  F77_CALL(daxpy)(&psq, &one, XWX, &nrv, A, &nrv);

  /*  Copy the contents of -prob to yminp */
  F77_CALL(dcopy)(&n, prob, &nrv, yminp, &nrv);
  F77_CALL(dscal)(&n, &min1, yminp, &nrv);

  /*  Subtraction of vectors: yminp <- y - prob */
  F77_CALL(daxpy)(&n, &one, y, &nrv, yminp, &nrv);

  /*  Computation of: aux <- XWX*theta */
  F77_CALL(dgemv)(notrans,&p,&p,&one, XWX, &p, theta, &nrv, &zero, aux, &nrv);

  /*  Computation of: aux <- trans(X)*yminp + aux */
  F77_CALL(dgemv)(trans, &n, &p, &one, X, &n, yminp, &nrv, &one, aux, &nrv);

  /*  Elimination of zeroes in matrix A */
  if (ap < p) {
      k = 0;
      for (i=0; i < ap; i++) {
	  for (j=0; j < ap; j++) {
	      A[k] = A[(i*p)+j];
	      k++;
	  }
      }
  }

  /*  Solving the linear equation system: aux <- inv(A)*theta */
  loes_gls(A, aux, ap, pivot, info);

  /*  The solution is theta, assign the solution to it: theta <- aux */
  F77_CALL(dcopy)(&p, aux, &nrv, theta, &nrv);
}

/*  Compute the parameter vector of logit ridge via Newton-Raphson iterations */
double ridgecrit(double *X, double *W, double *P, double *WX, double *XWX,
		 double *A, double *y, double *prob, double *theta,double *aux,
		 double *yminp, double *penalty,
		 int *pivot, int n, int p, int ap, int info, int crit)
{
  /*  Meaning of the input variables */

  /*  X       (n x p) Design matrix:  intercept & cluster representatives */
  /*  W       (n x n) Weight matrix:  contains initial weights */
  /*  P       (p x p) Penalty matrix: (non)-std ridge penalty, var dependent */
  /*  WX      (n x p) Auxiliary matrix: empty */
  /*  XWX     (p x p) Auxiliary matrix: empty */
  /*  A       (p x p) Auxiliary matrix: has the same content as P */
  /*  y       (n x 1) Binary response vector */
  /*  prob    (n x 1) Conditional probability vector: contains initial prob's */
  /*  theta   (p x 1) Parameter vector: contains initial values */
  /*  aux     (p x 1) Auxiliary vector: empty */
  /*  yminp   (n x 1) Auxiliary vector: contains -prob */
  /*  pivot   (p x 1) Pivot vector: empty */
  /*  n       (1)     The number of observations */
  /*  p       (1)     The final number of parameters in the model */
  /*  ap      (1)     The actual number of parameters in the model */
  /*  info    (1)     Auxiliary variable for solving the linear equations */
  /*  crit    (1)     Stands for the loss function to choose

  /*  Initialisation of variables */
  static int     nrv  = 1; /* for daxpy(): Not in ReVerse order */
  static double  one  = 1.;
  static double  zero = 0.;
  static char *notrans = "N";/* for dgemv(): Matrix is not transposed */
  int i;

  /*  Computing the logit ridge coefficients */
  ridgecoef(X,W,P,WX,XWX,A,y,prob,theta,aux,yminp,pivot,n,p,ap,info);

  /*  Compute the linear predictor and store it in yminp: yminp <- X*theta */
  F77_CALL(dgemv)(notrans,&n,&ap,&one, X, &n, theta, &nrv, &zero, yminp,&nrv);

  /*  Update prob: prob <- 1/(1+exp(-X*theta)) */
  for (i=0; i < n; i++)
      prob[i] = 1/(1 + exp(-yminp[i]));

  /*  Compute the penalty: penalty <- trans(theta)*P*theta */
  txwx(P, theta, penalty, aux, p, nrv);

  /*  Returning the penalized goodness of fit measure */
  return(loss(prob, y, yminp, n, crit) + *penalty);
}



/*  This is the supervised clustering algorithm */
void clusterer(double *E, double *X, double *Wini, double *P,
	       double *y, double *probini, double *thetaini,
	       double lambda,
	       int n, int g, int m, int p,
	       int critflag, int penflag, int blockflag,
	       int valiflag, int traceflag, int *genliste, double *kriterium)
{
/*  Explaining the arguments:
 *
 *  Name      Dim.	Content
 *
 *  E         (n x g) Gene expression matrix plus clinical covariates
 *  X         (n x p) Design matrix:  intercept & cluster representatives
 *  Wini      (n x n) Weight matrix:  contains initial weights
 *  P         (p x p) Penalty matrix: (non)-std ridge penalty, variance dep
 *  y         (n x 1) Binary response vector
 *  probini   (n x 1) Probability vector: contains initial probabilities
 *  thetaini  (p x 1) Parameter vector: contains initial values
 *  lambda       -    The penalty parameter lambda
 *  n            -    The number of observations
 *  g            -    The number of genes AND clinical covariates
 *  m            -    The number of clinical variables
 *  p            -    The final number of parameters in the model
 *  critflag     -    Loss function to choose: 0 is log-likeli, 1 is l2-loss
 *  penflag      -    Penalty to choose: 0 is standard, 1 is non-standard,
 *  blockflag    -    Blocking genes for the same cluster: 0 no, 1 yes
 *  valiflag     -    Should backward deletion be done: 0 no, 1 yes
 *  traceflag    -    Controls printed output: 0 none, 1 little, 2 detailed
 *  genliste  (q x 1) Long empty vector for recording the clustered genes
 *  kriterium (q x 1) Long empty vector for recording the criterion
 */

/* Matrix access macros */

/* Column _i_ of matrix X with  n  rows: */
#define Xc(_i_)  X + (_i_)*n
/* Diagonal element _i_ of matrix P[,] with  p  rows: */
#define Diag(P, _i_) P[(_i_)*(p+1)]
/*		      == _i_*p + _i_ */

  /*  Initialization of local integer variables */
  static int info = 0;
  static int nrv  = 1; /* for daxpy(): Not in ReVerse order */
  static double min1 = -1.0;

  int psq = p*p;
  int nsq = n*n;
  int noc = p-1;

  int counter = 0;
  int exclude = -1;
  int i, j, k, a, b, j_opt = -1;/*Wall*/

  double bval, bvalold;
  double rk, bk, crit;
  double penalty;

  /*  Memory allocation for local integer arrays */
  int *pivot = (int *) R_alloc(p, sizeof(int));
  int *allow = (int *) R_alloc(g, sizeof(int));
  int *clalw = (int *) R_alloc(g, sizeof(int));
  int *gls   = (int *) R_alloc(g, sizeof(int));

  /*  Memory allocation for local double arrays */
  double *W       = (double *) R_alloc(nsq,   sizeof(double));
  double *WX      = (double *) R_alloc(n*p,   sizeof(double));
  double *XWX     = (double *) R_alloc(psq,   sizeof(double));
  double *A       = (double *) R_alloc(psq,   sizeof(double));
  double *clrep   = (double *) R_alloc(n,     sizeof(double));
  double *prob    = (double *) R_alloc(n,     sizeof(double));
  double *theta   = (double *) R_alloc(p,     sizeof(double));
  double *aux     = (double *) R_alloc(p,     sizeof(double));
  double *yminp   = (double *) R_alloc(n,     sizeof(double));


  /*  Initialization of the allowance control vector for clinical covariates */
  for (j=0; j < g; j++)
      clalw[j] = 1;

  /*  Starting the loop for getting the clusters */
  for (i=0; i < noc; i++)
    {
      k = 1;
      rk   = 1.0/k;
      bk   = k/(k+1.0);
      bval = HUGE_VAL;

      /*  Getting the first gene or clinical covariate */
      for (j=0; j < g; j++)
	  /*  Only needs to be done if clinical covariate wasn't used before */
	  if (clalw[j]) {
	      /*  Initialization of the vector that records clustered genes */
	      allow[j] = 1;

	      /*  Updating of X */
	      F77_CALL(dcopy)(&n, E+(j*n), &nrv, Xc(i+1), &nrv);

	      /*  Updating of P */
	      Diag(P, i+1) = lambda*var(Xc(i+1), n);

	      /*  Updating of A, i.e. copy P to A */
	      F77_CALL(dcopy)(&psq, P, &nrv, A, &nrv);

#define UPDATE_and_RIDGE_CRIT						\
	      /*  Updating of A, i.e. copy P to A */			\
	      F77_CALL(dcopy)(&psq, P, &nrv, A, &nrv);			\
									\
	      /*  Copy the initial weights to W */			\
	      F77_CALL(dcopy)(&nsq, Wini, &nrv, W, &nrv);		\
									\
	      /*  Copy the initial probabilities to prob */		\
	      F77_CALL(dcopy)(&n, probini, &nrv, prob, &nrv);		\
									\
	      /*  Copy the initial parameter vector to theta */		\
	      F77_CALL(dcopy)(&p, thetaini, &nrv, theta, &nrv);		\
									\
	      /*  Initialize yminp */					\
	      F77_CALL(dcopy)(&n, probini, &nrv, yminp, &nrv);		\
	      F77_CALL(dscal)(&n, &min1, yminp, &nrv);			\
									\
	      /*  Computation of the loss criterion */			\
	      crit = ridgecrit(X, W, P, WX, XWX, A, y, 			\
			       prob, theta, aux, yminp, &penalty, 	\
			       pivot, n, p,(i+2), info, critflag)


	      UPDATE_and_RIDGE_CRIT;

	      if (crit < bval) { /* the current gene is the best */
		  j_opt = j;
		  bval = crit;
	      }

	      /*  Update the cluster representative clrep */
	      /*  HIER EVTL UNTERSCHEIDEN ZWISCHEN GENEN UND KLINVAR */
	      F77_CALL(dcopy)(&n, E+(j_opt*n), &nrv, clrep, &nrv);
	  }

      /*  Printed output for run-time control  */
      if (traceflag == 1) Rprintf(".");

      /*  Update gene and criterion lists : */
#define UPDATE_GLIST				\
      bvalold = bval;				\
      gls[k-1] = j_opt;				\
      genliste[counter]  = -1;			\
      kriterium[counter] = -1;			\
      counter++;				\
      genliste[counter]  = j_opt+1;		\
      kriterium[counter] = bval;		\
      counter++


      UPDATE_GLIST;

      if (j_opt >= (g-m)) { /* Updating if a clinical covariate was picked */

	  /*  Print the best covariate */
	  if (traceflag >= 2)
	      Rprintf("Added covariate %i, crit %f\n", j_opt+1, bval);

	  /*  Block the chosen covariate  */
	  clalw[j_opt] = 0;
      }
      else { /* Run clustering if a gene (not clinical covariate) was picked */

	  /*  Print the best gene */
	  if (traceflag >= 2)
	      Rprintf("Added gene %i, criterion value %f\n", j_opt+1, bval);

	  /*  Block the chosen gene (if blocking is activated) */
	  if(blockflag) allow[j_opt] = 0;

	  /*  Getting further genes */
	  while (1)
	  {
	      /*  Reinitialize bval for identification of the best gene */
	      bval = HUGE_VAL;

	      /*  The loop over all (non-blocked) genes */
	      for (j=0; j < (g-m); j++)
		  /*  Computations are only necessary if gene is non-blocked */
		  if (allow[j])
		  {
		      /*  Copy clrep to the design matrix X */
		      F77_CALL(dcopy)(&n, clrep, &nrv, Xc(i+1), &nrv);

		      /*  Updating of the design matrix X */
		      F77_CALL(daxpy)(&n, &rk, E+(j*n), &nrv, Xc(i+1), &nrv);
		      F77_CALL(dscal)(&n, &bk, Xc(i+1), &nrv);

		      /*  Updating of P */
		      Diag(P, i+1) = (penflag==0) ?
			  (lambda*var(Xc(i+1), n)) :
			  (lambda*var(Xc(i+1), n))/(k+1);

		      UPDATE_and_RIDGE_CRIT;

		      if (crit < bval) { /* the current gene is the best */
			  j_opt = j;
			  bval = crit;
		      }
		  } /* for(j .. ) if( allow[j] ) */

	      /*  Printed output for run-time control  */
	      if (traceflag==1) Rprintf(".");

	      /* Leave loop if the loss criterion has not improved : */
	      if (bval >= bvalold)
		  break;

	      /*-- else: Updating -- */

	      /*  Print the best gene and the criterion */
	      if (traceflag >= 2)
		  Rprintf("Added gene %i, criterion value %f\n", j_opt+1, bval);

	      /*  Updating of the cluster representative clrep */
	      F77_CALL(daxpy)(&n, &rk, E+(j_opt*n), &nrv, clrep, &nrv);
	      F77_CALL(dscal)(&n, &bk, clrep, &nrv);

	      /*  Updating of clustersize, genelist and criterion */
	      k++;
	      UPDATE_GLIST;

	      /*  Block the chosen gene (if activated) */
	      if(blockflag) allow[j_opt] = 0;

	      /*  VALIDATION STEP */
	      if (valiflag)
	      {
		  for (a=0; a < (k-2); a++)
		  {
		      /*  Initialization, copy the cluster mean to X */
		      F77_CALL(dcopy)(&n, E+(gls[(a == 0) ? 1 : 0]*n),
				      &nrv, Xc(i+1), &nrv);

		      /*  Still building the cluster mean in X */
		      for (b=1; b < k-1; b++)
		      {
			  rk = 1.0/b;
			  bk = b/(b+1.0);
			  F77_CALL(daxpy)(&n, &rk, E+(gls[(b < a)? b : b+1]*n),
					  &nrv, Xc(i+1), &nrv);
			  F77_CALL(dscal)(&n, &bk, Xc(i+1), &nrv);
		      }

		      /*  Updating of P */
		      Diag(P, i+1)=(penflag==0) ?
			  (lambda*var(Xc(i+1),n)) :
			  (lambda*var(Xc(i+1),n))/(k-1);

		      UPDATE_and_RIDGE_CRIT;

		      if (crit < bval) { /* exclusion of the gene improves */
			  exclude = a;
			  bval    = crit;
			  j_opt   = gls[exclude];
		      }
		  }

		  /*  Exclude a gene (the best found above) from the cluster */
		  if (exclude > (-1))
		  {
		      /*  Print which gene it was */
		      if (traceflag >= 2)
			  Rprintf("Delete gene %i, criterion value %f\n",
				  j_opt+1, bval);

		      /*  Storing the information & updating cluster size */
		      k--;
		      genliste[counter]  = -2;
		      kriterium[counter] = -2;
		      counter++;
		      genliste[counter]  = j_opt+1;
		      kriterium[counter] = bval;
		      counter++;

		      /*  Updating of the genelist */
		      for (b=0; b < k; b++)
			  gls[b] = gls[(exclude > b) ? b : (b+1)];

		      /*  Recompute the cluster representative */
		      F77_CALL(dcopy)(&n, E+(gls[0]*n), &nrv, clrep,&nrv);
		      for (b=1; b < k; b++)
		      {
			  rk = 1.0/b;
			  bk = b/(b+1.0);
			  F77_CALL(daxpy)(&n, &rk, E+(gls[b]*n),
					  &nrv, clrep, &nrv);
			  F77_CALL(dscal)(&n, &bk, clrep, &nrv);
		      }

		      /*  Reinitialization  */
		      rk      = 1.0/k;
		      bk      = k/(k+1.0);
		      bvalold = bval;
		      exclude = -1;

		      /*  Unblock the gene (if activated) */
		      if (blockflag) allow[j_opt] = 1;
		  }
	      }
	      rk = 1.0/k;
	      bk = k/(k+1.0);

	  } /* while(1) - getting more genes */

      } /* if / else "gene" (not clinical) */

      /* Updating for clinical && gene variables : */
      F77_CALL(dcopy)(&n, clrep, &nrv, Xc(i+1), &nrv);
      Diag(P, i+1) = (penflag == 0) ?
	  (lambda*var(Xc(i+1), n)) :
	  (lambda*var(Xc(i+1), n))/k;

      /*  The cluster is terminated */
      if (traceflag >= 1) Rprintf("\nCluster %i terminated \n", i+1);
      genliste[counter]  = -4;
      kriterium[counter] = -4;
      counter++;

    } /* for( cluster i ) */

}/* clusterer() */

void R_clusterer(double *E, double *X, double *Wini, double *P, double *y,
		 double *probini, double *thetaini, double *lambda, int *n,
		 int *g, int *m, int *p, int *critflag, int *penflag,
		 int *blockflag, int *valiflag, int *traceflag, int *genliste,
		 double *kriterium)
{
  clusterer(E, X, Wini, P, y, probini, thetaini, *lambda, *n, *g, *m, *p,
	    *critflag, *penflag, *blockflag, *valiflag, *traceflag,
	    genliste, kriterium);
}
