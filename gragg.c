/*
 * gragg.c
 * openmp algorithm for Gragg method
 * Authors: Ibrahim BAMBA
 * 	    Paul NITCHEU
 */

#include "gragg.h"

#define APPROX pow(10, -6)

/**
 * For printing execution time
 */
double my_gettimeofday(){
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

/* Functions for tests */

/**
 * Compute the min of two values
 */
double min(double a, double b){
	return a<b ? a : b;
}

/**
 * Compute the norm of a vector
 */
double norm_2(double *v, int n) {
	double sum=0;
	for(int i=0; i<n; i++)
		sum += v[i]*v[i];
	return sqrt(sum);
}

/**
 * Normalize the vector v
 */
void normalize(int n, double *v){
	double norm=norm_2(v, n);
	for(int i=0; i<n; i++)
		v[i] = v[i]/norm;
}

/**
 * Generate values of u
 */
void generate_U(int n, double* u) {
	srand(time(NULL));
	for(int i=0; i<n ; i++) {
		u[i]=rand();
	}
	normalize(n, u);
}

/**
 * Compute delta n+1 when rho > 0 or delta0 when rho < 0
 */
double compute_complete_delta(double delta, double *u, double rho, int N) {
	double res=0;
	for(int i=0; i<N; i++) {
		res+=u[i]*u[i];
	}
	return delta + res*rho;
}

/**
 * Generate deltas such as all delta_i >= 0 and delta_i < delta_i+1 for all i
 */
void generate_deltas(int n, double *u, double *d, double beta, double rho) {
	double x;
	srand(time(NULL));
	int beg, end;
	if(rho > 0) {
		beg=0;
		end=n-2;
	} else {
		beg=1;
		end=n-1;
	}
	for(int i=beg; i<=end; i++) {
		x=(rand()/(double)RAND_MAX) + i;
		d[i]=x+beta;
	}
	if(rho > 0) d[n-1]=compute_complete_delta(d[n-2], u, rho, n-1);
	else if(rho < 0) d[0]=compute_complete_delta(d[1], u, rho, n-1);
	else {
		perror("generation des deltas : rho null\n");
		exit(1);
	}
}

void print_test_result(results_t results, char* method, int numtest, double *u, double *d, double rho, int N, bool print_roots) {
	double y;
	double error=0;

	printf("/*************************Tests for %s***************************/\n", method);
	printf("\n--------------------\n"
			"TEST NUMBER %d :\n"
			"rho = %0.6g ", numtest, rho);
	rho<0 ? printf("< 0\n") : printf("> 0\n");
	printf("N = %d\n"
			"--------------------\n", N);

	for(int i=0; i<N; i++) {
		if(print_roots) printf("> root on (%.6g, %.6g) find by %s : %.6g\n", d[i], d[i+1], method, results.lambdas[i]);
		y=compute_f(u, d, fabs(rho), results.lambdas[i], N);
		if(print_roots) printf("  Compute f on %.6g : f(%.6g) = %.6g\n", results.lambdas[i], results.lambdas[i], y);
		if(!isnan(y)) error+=y;
	}
	printf("\nErreur moyen sur f : %.6g\n", error/N);
	printf("Nombre moyen d'itération : %d\n", results.iter);
}

/* COMPUTE F AND HIS DERIVATIVES */

double compute_f(double *u, double *d, double rho, double lambda_j, int N) {
	double res1=0, res2=0;
	int j=0;
	/* somme des termes négatifs */
	while(d[j] < lambda_j && j<N) {
		res1+=u[j]*u[j] / (d[j] - lambda_j);
		j++;
	}
	/* somme des termes positifs */
	for(; j<N; j++) {
		res2+=u[j]*u[j] / (d[j] - lambda_j);
	}
	return 1 + rho * (res1+res2);
}

double compute_fprim(double *u, double *d, double rho, double lambda_j, int N) {
	double res1=0, res2=0;
	int j=0;
	/* somme des termes négatifs */
	while(d[j] < lambda_j && j<N) {
		res1+=u[j]*u[j] / ((d[j] - lambda_j)*(d[j] - lambda_j));
		j++;
	}
	/* somme des termes positifs */
	for(; j<N; j++) {
		res2+=u[j]*u[j] / ((d[j] - lambda_j)*(d[j] - lambda_j));
	}
	return rho * (res1+res2);
}

double compute_fsnd(double *u, double *d, double rho, double lambda_j, int N) {
	double res1=0, res2=0;
	int j=0;
	/* somme des termes négatifs */
	while(d[j] < lambda_j && j<N) {
		res1+=2*u[j]*u[j] / ((d[j] - lambda_j)*(d[j] - lambda_j)*(d[j] - lambda_j));
		j++;
	}
	/* somme des termes positifs */
	for(; j<N; j++) {
		res2+=2*u[j]*u[j] / ((d[j] - lambda_j)*(d[j] - lambda_j)*(d[j] - lambda_j));
	}
	return rho * (res1+res2);
}

double g(double *d, double *u, double rho, double x, int N, int k){
	double res1=0, res2=0;
	int j=0;
	/* sum of positives terms */
	while(d[j]<x && j<N) {
		if(j!=k && j!=k+1)
			res1+=u[j]*u[j]/(d[j]-x);
		j++;
	}
	/* sum of negatives terms */
	for(; j<N; j++) {
		if(j!=k && j!=k+1)
			res2+=u[j]*u[j]/(d[j]-x);
	}
	return rho + res1 + res2 ;
}

double h(double *d, double *u, double rho, double x, int N, int k){
	return u[k]*u[k]/(d[k] - x) + u[k+1]*u[k+1]/(d[k+1] - x);
}

/**
 * Find the initial start guess
 */
double find_initial_guess(double *d, double *u, double rho, int N, int k) {
	double a, b, c, x, sum;
	double DELTA;
	DELTA=d[k+1] - d[k];
	x=(d[k]+d[k+1])/2;
	/* Case k on [1, N[ */
	if(k != N-1) {
		c=g(d, u, rho, x, N, k);
		if(compute_f(u, d, rho, x, N) >= 0) {
			a=c*DELTA+(u[k]*u[k] + u[k+1]*u[k+1]);
			b=u[k]*u[k]*DELTA;
			return a <=0 ? d[k] + (a - sqrt(fabs(a*a - 4*b*c))) / 2*c : d[k] +  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
		} else {
			a=-c*DELTA+(u[k]*u[k] + u[k+1]*u[k+1]);
			b=-u[k+1]*u[k+1]*DELTA;
			return a <= 0 ? d[k+1] + (a - sqrt(fabs(a*a - 4*b*c))) / 2*c : d[k+1] +  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
		}
	}
	/* Case k=N-1, last interval */
	else {
		c=g(d, u, rho, x, N, N-1);
		if(compute_f(u, d, rho, x, N) <= 0) {
			if(c + h(d, u, rho, d[N-1], N, N-1) <= 0) {
				sum=0;
				for(int i=0; i<N; i++) {
					sum+=u[i]*u[i];
				}
				return d[N] + sum/rho;
			} else {
				a= -c*DELTA+(u[N-2]*u[N-2] + u[N-1]*u[N-1]);
				b= -u[N-1]*u[N-1]*DELTA;
				return a >= 0 ? d[N] + (a + sqrt(fabs(a*a - 4*b*c))) / 2*c : d[N] +  2*b / (a - sqrt(fabs(a*a - 4*b*c)));
			}
		}
		else {
			a= -c*DELTA+(u[N-2]*u[N-2] + u[N-1]*u[N-1]);
			b= -u[N-1]*u[N-1]*DELTA;
			return a >= 0 ? d[N] + (a + sqrt(fabs(a*a - 4*b*c))) / 2*c : d[N] +  2*b / (a - sqrt(fabs(a*a - 4*b*c)));
		}
	}
}


/* COMPUTE OF THE INTERPOLATE EQUATION PARAMS */

/**
 * Compute gamma (i.e s)
 */
double compute_gamma(double *u, double *d, double lambda_j, int N, int k) {
	double res1=0, res2=0;
	int i=0;
	/* sum of positives terms */
	while(d[i]<lambda_j && i<N) {
		if(i!=k && i!=k+1)
			res1+=u[i]*u[i]*(d[i]-d[k+1])/pow(d[i]-lambda_j, 3);
		i++;
	}
	/* sum of negatives terms */
	for(; i<N; i++) {
		if(i!=k && i!=k+1)
			res2+=u[i]*u[i]*(d[i]-d[k+1])/pow(d[i]-lambda_j, 3);
	}
	return u[k]*u[k] + (pow(d[k] - lambda_j, 3)/(d[k] - d[k+1])) * (res1 + res2);
}

/**
 * Compute omega0 (i.e S)
 */
double compute_omega0(double *u, double *d, double lambda_j, int N, int k) {
	double res1=0, res2=0;
	int i=0;
	/* sum of positives terms */
	while(d[i]<lambda_j && i<N) {
		if(i!=k && i!=k+1)
			res1+=u[i]*u[i]*(d[i]-d[k])/pow(d[i]-lambda_j, 3);
		i++;
	}
	/* sum of negatives terms */
	for(; i<N; i++) {
		if(i!=k && i!=k+1) {
			res2+=u[i]*u[i]*(d[i]-d[k])/pow(d[i]-lambda_j, 3);
		}
	}
	return u[k+1]*u[k+1] + (pow(d[k+1] - lambda_j, 3)/(d[k+1] - d[k])) * (res1 + res2);
}

/**
 * Compute omega1 (i.e c)
 */
double compute_omega1(double *u, double *d, double rho, double lambda_j, int N, int k) {
	/* DELTA_k = delta_k - lambda_j
	 * DELTA_k+1 = delta_k+1 - lambda_j
	 */
	double DELTA_k = d[k] - lambda_j, DELTA_k_plus_1 = d[k+1] - lambda_j;
	return compute_f(u, d, rho, lambda_j, N) - (DELTA_k + DELTA_k_plus_1) * compute_fprim(u, d, rho, lambda_j, N) + DELTA_k * DELTA_k_plus_1 * (compute_fsnd(u, d, rho, lambda_j, N)/2);
}


/* SOLVING THE EQUATION : gamma + omega0/(delta_k - lambda_j) + omega1/(delta_k+1 - lambda_j) = 0
 * with lambda the unknown value
 * 	(E) is equivalent to :
 * 	(E') : gamma*lambda^2 - (gamma*d_k + gamma*d_k+1 + omega0 + omega1)*lambda + (gamma*d_k*d_k+1 + omega0*d_k+1 + omega1*d_k) = 0
 * 	which is a second order equation
 * 	This function compute The discriminant of (E') and find lambda with :
 * 	discr = (gamma*d_k + gamma*d_k+1 + omega0 + omega1)^2 - 4*gamma*(gamma*d_k*d_k+1 + omega0*d_k+1 + omega1*d_k)
 * 	This discr must be 0 because f has only one zero.
 * 	So lambda = sqrt(discr) / (2*gamma)
 */

/**
 * Compute parameter a and b for the solution
 */
double compute_param_a(double *u, double *d, double rho, double lambda_j, int N, int k) {
	double DELTA_k = d[k] - lambda_j, DELTA_k_plus_1 = d[k+1] - lambda_j;
	return (DELTA_k + DELTA_k_plus_1) * compute_f(u, d, rho, lambda_j, N) - (DELTA_k * DELTA_k_plus_1 * compute_fprim(u, d, rho, lambda_j, N));
}

double compute_param_b(double *u, double *d, double rho, double lambda_j, int N, int k) {
	double DELTA_k = d[k] - lambda_j, DELTA_k_plus_1 = d[k+1] - lambda_j;
	return DELTA_k * DELTA_k_plus_1 * compute_f(u, d, rho, lambda_j, N);
}

/**
 * Find the correction nano
 */
double find_correction(double a, double b, double c) {
	return a > 0 ? 2*b / (a + sqrt(fabs(a*a - 4*b*c))) : (a - sqrt(fabs(a*a - 4*b*c))) / 2*c;
}

/**
 * Find the solution of the equation
 */
double find_next_lambda(double a, double b, double c, double lambda_j) {
	return lambda_j + find_correction(a, b, c);
}

bool stopcriteria_non_monotically(double *u, double *d, double rho, double lambda_j, int N, int k) {
	double e, sumk=0, sumk_plus1=0, tho;
	double epsilon=pow(10, -14); /* Machine roundoff */
	int K;
	if(fabs(lambda_j-d[k]) < fabs(lambda_j-d[k+1])) K=k;
	else K=k+1;
	tho=lambda_j - d[K];
	for(int j=0; j<k; j++) {
		sumk+=(k-j+6) * fabs(u[j]*u[j]/(d[j] - lambda_j));
	}
	for(int j=k; j<N; j++) {
		sumk_plus1+=(j-k+5) * fabs(u[j]*u[j]/(d[j] - lambda_j));
	}
	e=2*rho + sumk + sumk_plus1 + fabs(compute_f(u, d, rho, d[K] + tho, N));
	return fabs(compute_f(u, d, rho, d[K] + tho, N)) <= e*epsilon + epsilon * fabs(tho) * fabs(compute_fprim(u, d, rho, d[K] - tho, N));
}

/**
 * Find the lamda_k such as f(lambda_j) ~= 0
 * on the interval ]k, k+1[
 */
results_t zero_finder_gragg(double *u, double *d, double rho, int N, int k) {
	results_t res;
	double lambda_j_moins_1, lambda_j, lambda_j_moins_2;
	double a, b, c;
	res.iter=0;
	bool stopcriteria=false;
	lambda_j = find_initial_guess(d, u, rho, N, k);

	do {
		/* Compute of the params */
		a=compute_param_a(u, d, rho, lambda_j, N, k);
		b=compute_param_b(u, d, rho, lambda_j, N, k);
		c=compute_omega1(u, d, rho, lambda_j, N, k);

		/* Stop criteria */
		/* Case of intermediates intervals */
		if((rho > 0 && k != N-1) || (rho < 0 && k != 0)) {
			/* updating lambdas and stop criteria */
			if(res.iter==0) lambda_j_moins_2=lambda_j;
			else if(res.iter==1) lambda_j_moins_1=lambda_j;
			else {
				lambda_j_moins_2=lambda_j_moins_1;
				lambda_j_moins_1=lambda_j;
				stopcriteria=(lambda_j_moins_1 - lambda_j_moins_2) * (lambda_j - lambda_j_moins_1) <= 0;
//				stopcriteria=stopcriteria_non_monotically(u, d, rho, lambda_j, N, k);
			}
		}
		/* Case of added interval */
		else stopcriteria=stopcriteria_non_monotically(u, d, rho, lambda_j, N, k);
		lambda_j=find_next_lambda(a, b, c, lambda_j);
		res.iter++;
	} while(!stopcriteria && res.iter<7);
	res.lambda=lambda_j;
	return res;
}

/**
 * Gragg method for finding eigenvalues of a matrix
 */
results_t gragg(double *u, double *d, double rho, int N) {
	results_t tmp, res;
	res.iter=0;
	double* lambdas = (double*)(calloc(sizeof(double), N));
#pragma omp parallel for private(tmp) schedule(runtime) ordered
	for(int k=0; k<N; k++) {
		tmp=zero_finder_gragg(u, d, fabs(rho), N, k);
		res.iter+=tmp.iter;
		lambdas[k]=tmp.lambda;
	}
	res.lambda=res.iter;
	res.iter/=N;
	res.lambdas=lambdas;
	return res;
}
