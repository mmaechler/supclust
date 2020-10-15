
#include <R.h>
/* for R_alloc, sorting, ... */

// In  ./pelora.c : ---------------------------------------------------------------------

void ridgecoef(double *X, double *W, double *P, double *WX, double *XWX,
	       double *A, double *y, double *prob, double *theta, double *aux,
	       double *yminp, int *pivot,
	       int n, int p, int ap, int info);

double ridgecrit(double *X, double *W, double *P, double *WX, double *XWX,
		 double *A, double *y, double *prob, double *theta,double *aux,
		 double *yminp, double *penalty,
		 int *pivot, int n, int p, int ap, int info, int crit);

void clusterer(double *E, double *X, double *Wini, double *P,
	       double *y, double *probini, double *thetaini,
	       double lambda,
	       int n, int g, int m, int p,
	       int critflag, int penflag, int blockflag,
	       int valiflag, int traceflag, int *genliste, double *kriterium);

// Called from R :
void R_clusterer(double *E, double *X, double *Wini, double *P, double *y,
		 double *probini, double *thetaini, double *lambda, int *n,
		 int *g, int *m, int *p, int *critflag, int *penflag,
		 int *blockflag, int *valiflag, int *traceflag, int *genliste,
		 double *kriterium);


// In  ./wilma.c : ----------------------------------------------------------------------

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
