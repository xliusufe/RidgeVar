#include <stdio.h>  	
#include <stdlib.h> 	
#include <string.h> 	
#include <math.h>
#include <float.h>
#include "var_hdh.h"

int LowTriangularInv(double *B, int n, double *A){ 
	int i,j,k;
	const double EPS=DBL_EPSILON;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<EPS)	return(0);
	for(i=0;i<n;i++)	B[i*n+i] = 1.0;
	for(j=1;j<n;j++)
		for(i=0;i<j;i++)	B[j*n+i] = 0.0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(fabs(A[i*n+i]-1)>EPS)
			for(j=i;j<n;j++)
				B[j*n+i] /= A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] -= A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){	
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;
	
	Z = (double*)malloc(sizeof(double)*n*p);
	znorm = (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];		
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p	
	for(j=1;j<p;j++){		
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) E[j*n+i] = Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) tmp1 += E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

void freematrix(double **x, int n)
{
	for (int i = 0; i < n; ++i)  free(x[i]);
	free(x);
}

double** requirematrix(int n, int p, int zero)
{
	double **x;
	x = (double**)malloc(sizeof(double*)*n);
	if (zero == 1) {
		for (int i = 0; i < n; ++i) {
			x[i] = (double*)calloc(p, sizeof(double));
			if (x[i] == NULL) fprintf(stderr, "Unable to allocate enough memory for array!\n");
		}
	} else {
		for (int i = 0; i < n; ++i) {
			x[i] = (double*)malloc(p*sizeof(double));
			if (x[i] == NULL) fprintf(stderr, "Unable to allocate enough memory for array!\n");
		}
	}
	return x;
}

double pythag(double a, double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

double SIGN(const double a, const double b)
{return (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

void tred2(double **z, double *d, double *e, int n, int yesvecs)
{
    int l,k,j,i;
    double scale,hh,h,g,f;
    for (i=n-1;i>0;i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0;k<i;k++)
                scale += fabs(z[i][k]);
            if (scale == 0.0)
                e[i]=z[i][l];
            else {
                for (k=0;k<i;k++) {
                    z[i][k] /= scale;
                    h += z[i][k]*z[i][k];
                }
                f=z[i][l];
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                z[i][l]=f-g;
                f=0.0;
                for (j=0;j<i;j++) {
                    if (yesvecs)
                        z[j][i]=z[i][j]/h;
                    g=0.0;
                    for (k=0;k<j+1;k++)
                        g += z[j][k]*z[i][k];
                    for (k=j+1;k<i;k++)
                        g += z[k][j]*z[i][k];
                    e[j]=g/h;
                    f += e[j]*z[i][j];
                }
                hh=f/(h+h);
                for (j=0;j<i;j++) {
                    f=z[i][j];
                    e[j]=g=e[j]-hh*f;
                    for (k=0;k<j+1;k++)
                        z[j][k] -= (f*e[k]+g*z[i][k]);
                }
            }
        } else
            e[i]=z[i][l];
        d[i]=h;
    }
    if (yesvecs) d[0]=0.0;
    e[0]=0.0;
    for (i=0;i<n;i++) {
        if (yesvecs) {
            if (d[i] != 0.0) {
                for (j=0;j<i;j++) {
                    g=0.0;
                    for (k=0;k<i;k++)
                        g += z[i][k]*z[k][j];
                    for (k=0;k<i;k++)
                        z[k][j] -= g*z[k][i];
                }
            }
            d[i]=z[i][i];
            z[i][i]=1.0;
            for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
        } else {
            d[i]=z[i][i];
        }
    }
}

void tqli(double **z, double *d, double *e, int n, int yesvecs)
{
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double EPS=DBL_EPSILON;
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                if (fabs(e[m]) <= EPS*dd) break;
            }
            if (m != l) {
                if (iter++ == 50) {
                    printf("Too many iterations in eigencmp.");
                    //exit(1);
                }
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    if (yesvecs) {
                        for (k=0;k<n;k++) {
                            f=z[k][i+1];
                            z[k][i+1]=s*z[k][i]+c*f;
                            z[k][i]=c*z[k][i]-s*f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}

void eigencmp(double **z, double *d, int n, int yesvecs)
{
    double *e = (double*)calloc(n, sizeof(double));
    tred2(z, d, e, n, yesvecs);
    tqli(z, d, e, n, yesvecs);
    free(e);
}

void matrixprod(double **x, double **y, double **xy, int n, int p, int q, int transpose) 
{
  int i,j,k;
  if (transpose == 0) {
  	for (i = 0; i < n; ++i) {
  		for (j = 0; j < q; ++j) {
  			xy[i][j] = 0;
  			for (k = 0; k < p; ++k) xy[i][j] += x[i][k]*y[k][j];
  		}
  	}
  } else if (transpose == 1) {
  	double **tx = requirematrix(n, p, 0);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < p; ++j) tx[i][j] = x[j][i];
	}
	matrixprod(tx, y, xy, n, p, q, 0);
	freematrix(tx, n);
  } else if (transpose == 2) {
  	double **ty = requirematrix(p, q, 0);
	for (i = 0; i < p; ++i) {
		for (j = 0; j < q; ++j) ty[i][j] = y[j][i];
	}
	matrixprod(x, ty, xy, n, p, q, 0);
	freematrix(ty, p);
  }
}

void solve_(double **Sn, double eta_k, int p, double **Sn_inv)
{
	int i,j,k;
	double **U = requirematrix(p, p, 0);
	double *singular = (double*)malloc(sizeof(double)*p);

	for (j = 0; j < p; ++j) {
		U[j][j] = Sn[j][j] + eta_k;
		for (i = j+1; i < p; ++i) 
			U[j][i] = U[i][j] = Sn[i][j];
	}
	eigencmp(U, singular, p, 1);
	for (j = 0; j < p; ++j) {
		for (i = j; i < p; ++i) {
			Sn_inv[i][j] = 0;
			for (k = 0; k < p; ++k)
				Sn_inv[i][j] += U[i][k]*U[j][k]/singular[k];
			if (i != j) Sn_inv[j][i] = Sn_inv[i][j];
		}
	}

	freematrix(U, p);
	free(singular);
}

double setup_eta(double *y, double *x, double *alpha, int *para)
{
	int n=para[0], p=para[1];
	double txy, maxtxy;
	maxtxy = 0.0;
	for (int i = 0; i < p; ++i) {
		txy = 0.0;
		for (int j = 0; j < n; ++j) txy += y[j]*x[j*p+i];
		maxtxy = FMAX(fabs(txy), maxtxy);
	}
	maxtxy = maxtxy/n/p*alpha[0];

	return maxtxy;
}

double setup_eta_cov(double *y, double *x, double *alpha, int *para)
{
	int i,j,k,n=para[0], p=para[1], q=para[2];
	double txy, maxtxy=0.0;
	for(k=0; k<q; k++){		
		for (i = 0; i < p; i++) {
			txy = 0.0;
			for (j = 0; j < n; j++) txy += y[k*n+j]*x[j*p+i];
			maxtxy = FMAX(fabs(txy), maxtxy);
		}		
	}
	maxtxy = maxtxy/n/p*alpha[0];
	return maxtxy;
}

double var(double *y, int n, int flag){
	int i;
	double meany = 0.0, variance = 0.0;
	for (i = 0; i < n; i++){ 
		meany += y[i];
		variance += y[i]*y[i];
	}
	meany /= n;
	variance /= n;
	variance -= meany*meany;
	if(flag) variance *= n/(n-1);
	return variance;
}


