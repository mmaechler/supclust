/*-*- mode: c; kept-old-versions: 12; kept-new-versions: 20 -*-*/


	   /************************************************/
	   /*						   */
	   /*  Supervised Clustering mit Score und Margin  */
	   /*  ------------------------------------------  */
	   /*						   */
	   /*	      Written by Marcel Dettling	   */
	   /*						   */
	   /************************************************/

/* M.Mächler (2003-06):
 *	o  Improve margin() ; much more clean up,
 *	o  final for(), and error checking, not while()
 *	o  "verbose" printing
 *	o  using "small epsilon" (relErr_margin) for margin comparison
 *	o  fixed get_new_gene() bug which did not always consider used[]
 *	o  return more to R and get rid of "-1" for end coding
 */

#include <R.h>
/* for R_alloc, sorting, ... */


/* ------ EXPORTS ---------*/
double margin(double x[], int n1, int n2);
void R_margin(double *x, int *n1, int *n2, double *result);

int score(double x[], double x_srt[], int indres[], int ind_srt[], int n);
void R_score(double *x, int *indres, int *n, double *result);

double get_new_gene(double x[], double y[], int indres[], int cluster_size,
		    int n, int n1, int n2, int p, int used[],
		    double x_srt[], int ind_srt[],
		    int gen_ids[], double cluster_mittel[],
		    double score_v[], double margin_v[],
		    int verbose);

void R_multicluster(double *y, int *resp,
		    int *n, int *n1, int *n2,
		    int *p, int *used,
		    double *cluster_mean, int *startsize,
		    int *genes_in_cluster,
		    int *scores, double *margins,
		    int *once_per_clust, int *c_verbose);
/* ------ ------- ---------*/


static double maximum(double *x, int n)
{
/* return( max(x[0..(n-1)] ) */
  int i;
  double maxi = x[0];
  for(i = 1; i < n; i++)
      if(x[i] > maxi) maxi = x[i];
  return(maxi);
}

static double minimum(double *x, int n)
{
/* return( min(x[0..(n-1)] ) */
  int i;
  double mini = x[0];
  for(i = 1; i < n; i++)
      if(x[i] < mini) mini = x[i];
  return(mini);
}

/*  The margin() of x[] where x[1:n1] ~ K=0, x[(n1+1):(n1+n2)] ~ K=1 */
double margin(double *x, int n1, int n2)
{
    /* min(x[n1 .. (n1+n2 -1)]) - max(x[0 .. (n1 - 1)]) */
    return(minimum(x+n1, n2) - maximum(x, n1));
}

/*  Margin() to be callable from R : */
void R_margin(double *x, int *n1, int *n2, double *result)
{
    *result = margin(x, *n1, *n2);
}


/*  Compute "Wilcoxon" Score of a vector of length n */
int score(double x[], double  x_srt[],
	  int indres[], /* the 0 / 1 `INDicators' of the `RESponse' y */
	  int ind_srt[], int n)
{
  /*  Variablendefinition */
  int i, j, nlarger = 0;

  /*  Berechnung des Score */
  for (i=0; i < n; i++)	 {
    ind_srt[i] = indres[i];
      x_srt[i] = x[i];
  }
  rsort_with_index(x_srt, ind_srt, n);

  for (i=0; i < n; i++)
    for (j=i; j < n; j++)
	nlarger += (ind_srt[i] > ind_srt[j]);

  return(nlarger);
}

/*  Score() to be callable from R : */
void R_score(double *x, int *indres, int *n, double *result)
{
  /*  Auxiliary needed (for sorting) */
  int	 *ind_srt =    (int *) R_alloc(*n, sizeof(int));
  double *x_srt	  = (double *) R_alloc(*n, sizeof(double));

  /*  Funktionsaufruf */
  *result = score(x, x_srt, indres, ind_srt, *n);
}


/*  Auxiliary function for  get_new_gene(): */
static
int gene_finder(int gen_ids[], double margin_v[], int nr_j)
{
  int j, j_max;
  double marg, max_marg;
  /* Search for the margin-best variable ("gene")
   * only among those with (optimal score, i.e.,) gen_ids[0..(nr_j-1)] */
  max_marg = margin_v[gen_ids[0]];
  j_max = 0;
  for (j = 1; j < nr_j; j++) {
      marg = margin_v[gen_ids[j]];
      if(max_marg < marg) {
	  max_marg = marg; j_max = j;
      }
  }
  return gen_ids[j_max];
}

/* Diese Funktion identifiziert jenes Gen, das Score und Margin optimiert */
double get_new_gene(double x[], double y[], int indres[], int cluster_size,
		    int n, int n1, int n2, int p, int used[],
		    /* auxiliary: */
		    double x_srt[], int ind_srt[],
		    /* output:*/
		    int gen_ids[], double cluster_mittel[],
		    double score_v[], double margin_v[],
		    /* diagnostic printing: no new line in this function! */
		    int verbose)
{
    /*	Variablendefinition */
    int i, j, k, j_max = -1/* -Wall*/;
    double Isiz1 = 1/(cluster_size+1.), marg, max_marg = R_NegInf;

    /* New averaged ("expression") matrix, margin and max.margin
     *  only for un-used[] variables : */
    for (j=0, k = 0; j < p; j++) {
	if (!used[j]) {
	    for (i=0; i < n; i++)
		cluster_mittel[i] = Isiz1*(y[(j*n)+i] + cluster_size* x[i]);

	    marg = margin(cluster_mittel, n1, n2);
	    margin_v[k] = marg;
	    gen_ids [k] = j;/* Index for genes in 0:(p-1) */
	    if(max_marg < marg) {
		max_marg = marg; j_max = j;
	    }
	    k++;
	}
    }
    /* now,  k == number_of{ unused } */

    /* Best margin of un-used ones: if > 0, we've found the best variable/gene*/
    if (max_marg > 0)
    {
	j_max++;/* 0:(p-1) |--> 1:p */
	if(verbose) Rprintf("g_new_g(): best margin > 0 at %d", j_max);
	return j_max;
    }
    else { /* der Score für die ganze Matrix berechnet	*/
	for (j=0, k = 0; j < p; j++) {
	    if (!used[j]) {
		for (i=0; i < n; i++)
		    cluster_mittel[i] = Isiz1*(y[(j*n)+i] + cluster_size* x[i]);

		score_v[k] = score(cluster_mittel, x_srt, indres, ind_srt,  n);
		/* have gen_ids [k] == j */
		k++;
	    }
	}
	/* again, k == number_of{ unused } -- same as above */

	rsort_with_index(score_v, gen_ids, k);

	/* Falls bei einem Gen der tiefste Score-Wert angenommen wird: fertig */
	if (score_v[0] < score_v[1]) {
	    gen_ids[0]++;
	    if(verbose)
		Rprintf("g_new_g(): unique lowest score at %d", gen_ids[0]);
	    return(gen_ids[0]);
	}
	else {
	    /*	Sonst betrachte die Margins aller Gene mit minimalem Score  */
	    j=1;
	    while (score_v[0] == score_v[j])
		j++;

	    if(verbose)
		Rprintf("g_new_g(): j=%d > 1 minimal scores -> g_finder", j);
	    return(gene_finder(gen_ids, margin_v, j) + 1);
	}
    }
} /* get_new_gene() */


/*  Called from R : Does "forward search" (for 1 cluster only) */
void R_multicluster(double *y, int *resp,	/* [1:2] */
		    int *n, int *n1, int *n2,	/* [3:5] */
		    int *p,			/*   [6] */
		    int *used,			/*   [7] input AND output !*/
		    double *cluster_mean,	/*   [8] */
		    int *gl_size,		/*   [9] */
		    int *genes_in_cluster,	/*  [10] : `gic' with result!
						 * length = gl_size + p */
		    int *scores,		/*  [11] */
		    double *margins,		/*  [12] */
		    int *once_per_clust,	/*  [13] : */
		    /* `logical: if `true', each gene appears only once in a
		     *  cluster, i.e., only weights +/- 1 allowed */
		    int *c_verbose)		/*  [14] */
{
#define relErr_margin 1e-14

    /*	Variablendefinition  */
    int	 i, j, c, size;
    int	 t_score, old_score;
    Rboolean foundBest;

    int max_size = *gl_size + *p;
    int verbose	 = *c_verbose;

    double Isiz1, t_margin, old_margin;

    /* Allocate auxiliary vectors */
    int	 *ind_srt		= (int *) R_alloc(*n, sizeof(int));
    int	 *gen_ids		= (int *) R_alloc(*p, sizeof(int));

    double *x_srt		= (double *) R_alloc(*n, sizeof(double));
    double *cluster_mittel	= (double *) R_alloc(*n, sizeof(double));
    double *score_v		= (double *) R_alloc(*p, sizeof(double));
    double *margin_v		= (double *) R_alloc(*p, sizeof(double));

    /*	Vergrösserung des Clusters  */
    size = *gl_size;

    /*	Identifizieren des ersten Gens	*/
    if (size==0)
    {
	if(verbose)
	    Rprintf("R_multicluster(*, gl_size = 0 [_zero_])\n");

#define G_NEW(SIZE)							\
	j = get_new_gene(cluster_mean, y, resp, SIZE, *n, *n1, *n2, *p,	\
			 used, x_srt, ind_srt, gen_ids,			\
			 cluster_mittel, score_v, margin_v, verbose);	\
	genes_in_cluster[SIZE] = j;					\
									\
        j--;/* {1:p} |--> {0:(p-1)} */					\
	if(*once_per_clust) {						\
	    if(verbose) Rprintf(" used: %d", j+1);			\
	    used[j] = 1;						\
	}

	G_NEW(0);

	/* Cluster "mean" dieses Gens: */
	j *= (*n);
	for (i=0; i < *n; i++)
	    cluster_mean[i] = y[i+j];
    }
    else {
	if(verbose)
	    Rprintf("R_multicluster(*, gl_size = %d > 0)", size);
	/* leaves scores[0..(size-2)] and margins[..] undefined */
	size--;
    }

    t_score  = score (cluster_mean, x_srt, resp, ind_srt, *n);
    t_margin = margin(cluster_mean, *n1, *n2);
    scores [size] = t_score;
    margins[size] = t_margin;
    size++;
    if(verbose/* >= 2 */)
	Rprintf(";  sc()= %d, marg()= %17.15g\n", t_score, t_margin);

    /*	Cluster-Sequenz: Identifizieren weiterer Gene  */

    /* was  while(1)  { ... } */
    foundBest = FALSE;
    for(c = size; c < max_size; c++) {/* usually break LONG before c = max! */

        G_NEW(c);

	Isiz1 = 1/(c + 1.);
	j *= (*n);
	for (i=0; i < *n; i++)
	    cluster_mean[i] = Isiz1*(y[i + j] + c * cluster_mean[i]);

	old_score = t_score;
	old_margin= t_margin;
	t_score	 = score (cluster_mean, x_srt, resp, ind_srt, *n);
	t_margin = margin(cluster_mean, *n1, *n2);
	scores [c] = t_score;
	margins[c] = t_margin;
	if(verbose/* >= 2 */)
	    Rprintf(";  sc()= %d, marg()= %17.15g\n", t_score, t_margin);

	if (t_score  >	old_score ||
	    (t_score == old_score &&
	     old_margin - t_margin > - relErr_margin * t_margin))
	    /* t_margin < old_margin  is  not sufficient! */
	{
	    /* adding this one would not improve : */
	    foundBest = TRUE;
	    if(*once_per_clust) /* reset to "un-used" */
		used[genes_in_cluster[c] - 1] = 0;
	    break;
	}

    } /* for */

    if(!foundBest)
	REprintf("R_multicluster() __BUG__ : not foundBest!!");
    if(! *once_per_clust)
	/* Sperren schon benutzter Gene für den nächsten Cluster */
	for (i=0; i < c; i++)
	    used[genes_in_cluster[i] - 1] = 1;
    *gl_size = c;/* return the size */

    return;
} /* R_multicluster() */
