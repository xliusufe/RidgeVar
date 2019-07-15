#include <math.h> // required for sqrt(), fabs();
#include <stdio.h>  // required for exit
#include <stdlib.h> // required for malloc(), free();
#include <Rinternals.h>  // required for SEXP et.al.;
#include <float.h>  // required for DBL_EPSILON
#include <string.h> // required for memcpy()
#include <time.h>

#define FMAX(a,b) ((a) > (b) ? (a) : (b))

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

void matrixprod(double **x, double **y, double **xy, 
	int n, int p, int q, int transpose) 
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

double setup_eta(double *y, double **x, double *alpha, int *para)
{
	int n=para[0], p=para[1];
	double txy, maxtxy;
	maxtxy = 0.0;
	for (int i = 0; i < p; ++i) {
		txy = 0.0;
		for (int j = 0; j < n; ++j) txy += y[j]*x[j][i];
		maxtxy = FMAX(fabs(txy), maxtxy);
	}
	maxtxy = maxtxy/n/p*alpha[0];

	return maxtxy;
}

void var_rr(double *y, double **x, double eta_k, int *para, double *sigma2, int *debias)
{
	int n=para[0], p=para[1];
	int i, j;

	if (p<=n) {
		double **Sn = requirematrix(p, p, 0);
		double **Sn_inv = requirematrix(p, p, 0);
		double **temp = requirematrix(p, n, 0);
		double **An1 = requirematrix(n, n, 0);
		double An1_tr = 0.0, ynorm = 0.0, temp2 = 0.0;
		double *temp1 = (double*)calloc(n, sizeof(double));


		matrixprod(x, x, Sn, p, n, p, 1);
		for (i = 0; i < p; ++i) {
			for (j = 0; j < p; ++j) Sn[i][j] /= n;
		}

		solve_(Sn, eta_k, p, Sn_inv);
		matrixprod(Sn_inv, x, temp, p, p, n, 2);
		matrixprod(x, temp, An1, n, p, n, 0);
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				An1[i][j] /= n;
				temp1[i] += An1[i][j]*y[j];
			}
			An1_tr += An1[i][i];
			ynorm += pow(y[i], 2);
			temp2 += temp1[i]*y[i];
		}
		ynorm /= n;
		sigma2[0] = (ynorm - temp2/n) / (1-An1_tr/n);

		free(temp1);
		freematrix(Sn, p);
		freematrix(Sn_inv, p);
		freematrix(temp, p);
		freematrix(An1, n);
	} else {
		double **xtx = requirematrix(n, n, 0);
		double *lam = (double*)malloc(sizeof(double)*n);
		double *vy = (double*)calloc(n, sizeof(double));
		double An1_tr = 0.0, ynorm = 0.0, YAn1Y = 0.0, bias = 0.0;

		matrixprod(x, x, xtx, n, p, n, 2);
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) xtx[i][j] /= n;
		}

		eigencmp(xtx, lam, n, 1);

		for (j = 0; j < n; ++j) {
			for (i = 0; i < n; ++i) vy[j] += xtx[i][j]*y[i];
			An1_tr += 1/(eta_k+lam[j]);
			YAn1Y += pow(vy[j], 2)/(eta_k+lam[j]);
			ynorm += pow(y[j], 2);
		}
		if (debias[0] == 1) bias = eta_k*YAn1Y/An1_tr;
		sigma2[0] = YAn1Y/An1_tr - bias;

		freematrix(xtx, n);
		free(lam);
		free(vy);
	}
}

SEXP VAR_RR(SEXP Y_, SEXP X_, SEXP ETAK_, SEXP PARA_, SEXP DEBIAS_, SEXP ALPHA_, SEXP ISNULLETA_)
{
	// dimensions
	int *para = INTEGER(PARA_);
	int n     = para[0];
	int p     = para[1];

	// Pointers
	double **x = (double**)malloc(sizeof(double*)*n);
  	for (int i = 0; i < n; ++i) x[i] = &REAL(X_)[i*p];
	double *y  = REAL(Y_);
	double eta;

	// Outcome
	SEXP _output, _sigma2, _r_names;
  	PROTECT(_output = allocVector(VECSXP, 1));
  	PROTECT(_sigma2 = allocVector(REALSXP, 1));
  	PROTECT(_r_names   = allocVector(STRSXP, 1));

  	if (INTEGER(ISNULLETA_)[0] == 1) {
  		eta = setup_eta(y, x, REAL(ALPHA_), para);
  	} else {
  		eta = REAL(ETAK_)[0];
  	}
  	var_rr(y, x, eta, para, REAL(_sigma2), INTEGER(DEBIAS_));
  	free(x);

  	SET_STRING_ELT(_r_names, 0,  mkChar("sigma2"));
	SET_VECTOR_ELT(_output, 0, _sigma2);
	setAttrib(_output, R_NamesSymbol, _r_names); 

	UNPROTECT(3);
	return _output;
}



SEXP VAR_MM(SEXP Y_, SEXP X_, SEXP SIGMA_, SEXP PARA_, SEXP IDENTITY_, SEXP ISNULL_)
{
	// dimensions
	int *para = INTEGER(PARA_);
	int n     = para[0];
	int p     = para[1];

	// Pointers
	double **x = (double**)malloc(sizeof(double*)*n);
	double *y  = REAL(Y_);

	// Intermediate quantities
	int i, j;
	double ynorm = 0.0, xynorm = 0.0;
	double *temp = (double*)calloc(p, sizeof(double));

	// Outcome
	SEXP _output, _sigma2, _r_names;
  	PROTECT(_output = allocVector(VECSXP, 1));
  	PROTECT(_sigma2 = allocVector(REALSXP, 1));
  	PROTECT(_r_names   = allocVector(STRSXP, 1));

  	for (i = 0; i < n; ++i) {
  		x[i] = &REAL(X_)[i*p];
  		ynorm += pow(y[i], 2);
  	}
  	for (j = 0; j < p; ++j) {
  		for (i = 0; i < n; ++i) temp[j] += x[i][j]*y[i];
  	}
  	ynorm /= n;

    if (INTEGER(ISNULL_)[0] == 0) {
        double temp1;
        for (i = 0; i < p; ++i) {
            for (temp1 = j = 0; j < p; ++j) temp1 += REAL(SIGMA_)[j*p+i]*temp[j];
            xynorm += temp1*temp[i];
        }
        REAL(_sigma2)[0] = ((p+n+1)*ynorm - xynorm/n)/(n+1);
  	} else {
        if (INTEGER(IDENTITY_)[0] != 0) {
            for (j = 0; j < p; ++j) xynorm += temp[j]*temp[j];
            REAL(_sigma2)[0] = ((p+n+1)*ynorm - xynorm/n)/(n+1);
  		} else {
  			double m1 = 0.0, m2 = 0.0;
  			for (j = 0; j < p; ++j) xynorm += pow(temp[j], 2);
  			double **xtx = requirematrix(n, n, 0);
  			double *lam = (double*)malloc(sizeof(double)*n);

  			matrixprod(x, x, xtx, n, p, n, 2);
  			eigencmp(xtx, lam, n, 0);
  			for (i = 0; i < n; ++i) {
  				lam[i] /= n;
  				m1 += lam[i];
  				m2 += pow(lam[i], 2);
  			}
  			m2 = m2/p - pow(m1, 2)/p/n;
  			m1 /= p;
  			REAL(_sigma2)[0] = (1+p*pow(m1, 2)/m2/(n+1))*ynorm - m1*xynorm/n/m2/(n+1);
  			freematrix(xtx, n);
  			free(lam);
  		}
  	}

  	SET_STRING_ELT(_r_names, 0,  mkChar("sigma2"));
	SET_VECTOR_ELT(_output, 0, _sigma2); 
	setAttrib(_output, R_NamesSymbol, _r_names); 

  	free(x);
  	free(temp);
	UNPROTECT(3);
	return _output;
}

double var(double *y, int n)
{
	int i;
	double meany = 0.0, variance = 0.0;
	for (i = 0; i < n; ++i) meany += y[i];
	meany /= n;
	for (i = 0; i < n; ++i) variance += pow(y[i] - meany, 2);
	return variance/(n-1);
}

SEXP VAR_MLE(SEXP Y_, SEXP X_, SEXP PARA_, SEXP MAXITER_, SEXP TOL_)
{
	// dimensions
	int *para = INTEGER(PARA_);
	int n     = para[0];
	int p     = para[1];

	// Pointers
	double **x = (double**)malloc(sizeof(double*)*n);
	double *y  = REAL(Y_);

	// Intermediate quantities
	int i, j, step = 1;
	double theta0[2], theta[2], grad[2], sigma2, eta2, likelih0, likelih;
	double yUy, trUxtx, yUxtxUy, trUxtxUxtx, yUxtxUxtxUy, theta1;
	double **H = requirematrix(2, 2, 0);
	double **xtx = requirematrix(n, n, 0);
  	double *lam = (double*)malloc(sizeof(double)*n);
  	double *vy = (double*)calloc(n, sizeof(double));
  	double *vy2 = (double*)malloc(sizeof(double)*n);
  	double *temp = (double*)malloc(sizeof(double)*n);

	for (i = 0; i < n; ++i) x[i] = &REAL(X_)[i*p];
	theta0[0] = log(var(y, n));
	theta0[1] = 0.01;
	sigma2 = theta0[0];
	eta2 = theta0[1];

	matrixprod(x, x, xtx, n, p, n, 2);
	eigencmp(xtx, lam, n, 1);
	for (i = 0; i < n; ++i) {
		lam[i] /= p;
		for (j = 0; j < n; ++j) vy[i] += xtx[j][i]*y[j];
		vy2[i] = pow(vy[i], 2);
		temp[i] = eta2*lam[i]+1;
	}
	for (likelih0 = i = 0; i < n; ++i) likelih0 += vy2[i]/(temp[i]);
	likelih0 *= exp(-1.0*sigma2);
	likelih0 += n*sigma2;
	for (i = 0; i < n; ++i) likelih0 += log(temp[i]);

	while(step < INTEGER(MAXITER_)[0]) {
		step += 1;
		yUy = trUxtx = yUxtxUy = trUxtxUxtx = yUxtxUxtxUy = 0.0;
		for (i = 0; i < n; ++i) {
			yUy += vy2[i]/temp[i];
			trUxtx += lam[i]/temp[i];
			yUxtxUy += vy2[i]*lam[i]/pow(temp[i], 2);
			trUxtxUxtx += pow(lam[i]/temp[i], 2);
			yUxtxUxtxUy += pow(vy[i]*lam[i], 2)/pow(temp[i], 3);
		}
		theta1 = exp(-1.0*sigma2);
		grad[0] = n-theta1*yUy;
		grad[1] = trUxtx-theta1*yUxtxUy;
		H[0][0] = theta1*yUy;
		H[0][1] = theta1*yUxtxUy;
		H[1][1] = -1.0*trUxtxUxtx + 2*theta1*yUxtxUxtxUy;
		grad[1] -= grad[0]*yUxtxUy/yUy;
		H[1][1] -= H[0][1]*yUxtxUy/yUy;
		theta[1] = theta0[1] + grad[1]/H[1][1];
		theta[0] = theta0[0] + (grad[0]-H[0][1]*grad[1]/H[1][1])/H[0][0];
		sigma2 = theta[0];
		eta2 = theta[1];
		if (eta2 <= 0) {
			sigma2 = theta0[0]; 
			eta2 = theta0[1];
			break;
		}
		for (i = 0; i < n; ++i) temp[i] = eta2*lam[i]+1;
		for (likelih = i = 0; i < n; ++i) likelih += vy2[i]/(temp[i]);
		likelih *= exp(-1.0*sigma2);
		likelih += n*sigma2;
		for (i = 0; i < n; ++i) likelih += log(temp[i]);
		if (likelih > likelih0) {
			sigma2 = theta0[0]; 
			eta2 = theta0[1];
			break;
		} else {
			if (FMAX(fabs(theta[0]-theta0[0]), fabs(theta[1]-theta0[1])) < REAL(TOL_)[0]) {
				break;
			} else {
				theta0[0] = theta[0];
				theta0[1] = theta[1];
				likelih0 = likelih;
			}
		}
	}

	freematrix(xtx, n);
	freematrix(H, 2);
	free(lam);
	free(vy);
	free(vy2);
	free(temp);
	free(x);

	// Outcome
	SEXP _output, _sigma2, _r_names;
  	PROTECT(_output = allocVector(VECSXP, 1));
  	PROTECT(_sigma2 = allocVector(REALSXP, 1));
  	PROTECT(_r_names   = allocVector(STRSXP, 1));
  	REAL(_sigma2)[0] = exp(sigma2);

  	SET_STRING_ELT(_r_names, 0,  mkChar("sigma2"));
	SET_VECTOR_ELT(_output, 0, _sigma2);
	setAttrib(_output, R_NamesSymbol, _r_names); 

	UNPROTECT(3);
	return _output;
}

