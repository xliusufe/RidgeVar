#ifndef RBS_H_INCLUDED
#define RBS_H_INCLUDED

#define FMAX(a,b) ((a) > (b) ? (a) : (b))

int LowTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

void freematrix(double **x, int n);

double** requirematrix(int n, int p, int zero);

double pythag(double a, double b);

double SIGN(const double a, const double b);

void tred2(double **z, double *d, double *e, int n, int yesvecs);

void tqli(double **z, double *d, double *e, int n, int yesvecs);

void eigencmp(double **z, double *d, int n, int yesvecs);

void matrixprod(double **x, double **y, double **xy, int n, int p, int q, int transpose);

void solve_(double **Sn, double eta_k, int p, double **Sn_inv);

double setup_eta(double *y, double *x, double *alpha, int *para);

double setup_eta_cov(double *y, double *x, double *alpha, int *para);

double var(double *y, int n, int flag);

#endif // RBS_H_INCLUDED
