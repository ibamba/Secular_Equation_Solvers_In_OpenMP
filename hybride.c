/*
 * hybride.c
 * openmp algorithm for hybrid scheme
 * Authors: Ibrahim BAMBA
 * 			Paul NITCHEU
 */

#include "gragg.h"

double psi(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double res1=0, res2=0;
	int j=0;
	/* sum of negatives terms */
	while(d[j] < y && (fm ? j<k : j<=k)) {
		res1+=u[j]*u[j] / (d[j] - y);
		j++;
	}
	/* sum of positives terms */
	for(; (fm ? j<k : j<=k); j++) {
		res2+=u[j]*u[j] / (d[j] - y);
	}
	return rho * (res1+res2);

}

double psiprim(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double res1=0, res2=0;
	int j=0;
	/* sum of negatives terms */
	while(d[j] < y && (fm ? j<k : j<=k)) {
		res1+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
		j++;
	}
	/* sum of positives terms */
	for(; (!fm && j <= k) || (fm && j < k); j++) {
		res2+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
	}
	return rho * (res1+res2);
}

double fi(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double res1=0, res2=0;
	int j=(fm ? k+2 : k+1);
	/* sum of negatives terms */
	while(d[j] < y && j < N) {
		res1+=u[j]*u[j] / (d[j] - y);
		j++;
	}
	/* sum of positives terms */
	for(; j < N; j++) {
		res2+=u[j]*u[j] / (d[j] - y);
	}
	return rho * (res1+res2);
}

double fiprim(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double res1=0, res2=0;
	int j=(fm ? k+2 : k+1);
	/* sum of negatives terms */
	while(d[j] < y && j < N) {
		res1+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
		j++;
	}
	/* sum of positives terms */
	for(; j < N; j++) {
		res2+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
	}
	return rho * (res1+res2);
}

/**
 * Compute the secular function f(y) with
 * the mth term in the summation is removed
 */
double compute_fm(double *u, double *d, double rho, double y, int N, int m) {
	double res1=0, res2=0;
	int j=0;
	/* sum of negatives terms */
	while(d[j] < y && j<N) {
		if(j != m) res1+=u[j]*u[j] / (d[j] - y);
		j++;
	}
	/* sum of positive terms */
	for(; j<N; j++) {
		if(j != m) res2+=u[j]*u[j] / (d[j] - y);
	}
	return 1 + rho * (res1+res2);
}

/**
 * Compute fm first derivative
 */
double compute_fmprim(double *u, double *d, double rho, double y, int N, int m) {
	double res1=0, res2=0;
	int j=0;
	/* sum of negatives terms */
	while(d[j] < y && j<N) {
		if(j != m) res1+=u[j]*u[j] / (d[j] - y)*(d[j] - y);
		j++;
	}
	/* sum of positive terms */
	for(; j<N; j++) {
		if(j != m) res2+=u[j]*u[j] / (d[j] - y)*(d[j] - y);
	}
	return rho * (res1+res2);
}

/* sum of negatives terms of f */
/*double psi(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double sum=0;
	for(int j=0; (fm ? j<k : j<=k); j++) {
		sum+=u[j]*u[j] / (d[j] - y);
	}
	return rho * sum;
}
 */
/* sum of negatives terms of f' */
/*double psiprim(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double sum=0;
	for(int j=0; (fm ? j<k : j<=k); j++) {
		sum+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
	}
	return rho * sum;
}
 */
/* sum of positives terms of f */
/*double fi(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double sum=0;
	for(int j=(fm ? k+2 : k+1); j<N; j++) {
		sum+=u[j]*u[j] / (d[j] - y);
	}
	return rho * sum;
}
 */
/* sum of positives terms of f' */
/*double fiprim(double *u, double *d, double rho, double y, int N, int k, bool fm) {
	double sum=0;
	for(int j=(fm ? k+2 : k+1); j<N; j++) {
		sum+=u[j]*u[j] / ((d[j] - y)*(d[j] - y));
	}
	return rho * sum;
}
 */
/**
 * Compute the secular function f(y) with
 * the mth term in the summation is removed
 */
/*double compute_fm(double *u, double *d, double rho, double y, int N, int m) {
	return 1 + psi(u, d, rho, y, N, m, true) + fi(u, d, rho, y, N, m, true);
}
 */
/**
 * Compute fm first derivative
 */
/*double compute_fmprim(double *u, double *d, double rho, double y, int N, int m) {
	return 1 + psiprim(u, d, rho, y, N, m, true) + fiprim(u, d, rho, y, N, m, true);
}
 */
/**
 * Middle way method on f
 * @return c, s, and S
 */
double middle_way_f(double *u, double *d, double rho, double y, int N, int k) {
	double a, b, c, nano;
	double DELTA_k=d[k] - y;
	double DELTA_k_plus1=d[k+1] - y;

	if(k != N-1) {
		a=(DELTA_k + DELTA_k_plus1) * compute_f(u, d, rho, y, N) - DELTA_k * DELTA_k_plus1 * compute_fprim(u, d, rho, y, N);
		b=DELTA_k * DELTA_k_plus1 * compute_f(u, d, rho, y, N);
		c=compute_f(u, d, rho, y, N) - DELTA_k * psiprim(u, d, rho, y, N, k, false) - DELTA_k_plus1 * fiprim(u, d, rho, y, N, k, false);
		nano = a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
		return y + nano;
	}
	else { /* case k=n-1 : last interval */
		double DELTA_N = d[N-1] - y; /* N - 1 ???! */
		double DELTA_N_moins1 = d[N-2] - y;
		a=(DELTA_N_moins1 + DELTA_N) * compute_f(u, d, rho, y, N)  - DELTA_N_moins1 * DELTA_N * compute_fprim(u, d, rho, y, N);
		b=DELTA_N_moins1 * DELTA_N * compute_f(u, d, rho, y, N);
		c=compute_f(u, d, rho, y, N) - DELTA_N_moins1 * psiprim(u, d, rho, y, N, N-2, false) - u[N-1]*u[N-1]/DELTA_N;
		nano = a >= 0 ? (a + sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a - sqrt(fabs(a*a - 4*b*c)));
		return y + nano;
	}
}
/*
double* middle_way_f(double *u, double *d, double rho, double y, int N, int k) {
	double *res=(double*)malloc(sizeof(double)*3);
	double DELTA_k=d[k] - y;
	double DELTA_k_plus1=d[k+1] - y;

	res[1]=DELTA_k*DELTA_k*psiprim(u, d, rho, y, N, k, false);

	res[1]=DELTA_k_plus1*DELTA_k_plus1*fiprim(u, d, rho, y, N, k, false);

	if(k != N-1)
		res[0]=compute_f(u, d, rho, y, N) - DELTA_k * psiprim(u, d, rho, y, N, k, false) - DELTA_k_plus1 * fiprim(u, d, rho, y, N, k, false);
	else {
		double DELTA_N = d[N-1] - y;
		double DELTA_N_moins1 = d[N-2] - y;
		res[0]=compute_f(u, d, rho, y, N) - DELTA_N_moins1 * psiprim(u, d, rho, y, N, N-2, false) - u[N-1]*u[N-1]/DELTA_N;
	}

	return res;
}
 */

/**
 * Middle way method on fm
 * @return c, s and S
 */
double middle_way_fm(double *u, double *d, double rho, double y, int N, int k) {
	double a, b, c, nano;
	double DELTA_k=d[k] - y;
	double DELTA_k_plus1=d[k+1] - y;

	if(k != N-1) {
		a=(DELTA_k + DELTA_k_plus1) * compute_fm(u, d, rho, y, N, k) - DELTA_k * DELTA_k_plus1 * compute_fmprim(u, d, rho, y, N, k);
		b=DELTA_k * DELTA_k_plus1 * compute_fm(u, d, rho, y, N, k);
		c=compute_fm(u, d, rho, y, N, k) - DELTA_k * psiprim(u, d, rho, y, N, k, true) - DELTA_k_plus1 * fiprim(u, d, rho, y, N, k, true);
		nano = a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
		return y + nano;
	}
	else { /* case k=n-1 : last interval */
		double DELTA_N = d[N-1] - y; /* N - 1 ???! */
		double DELTA_N_moins1 = d[N-2] - y;
		a=(DELTA_N_moins1 + DELTA_N) * compute_fm(u, d, rho, y, N, k)  - DELTA_N_moins1 * DELTA_N * compute_fmprim(u, d, rho, y, N, k);
		b=DELTA_N_moins1 * DELTA_N * compute_fm(u, d, rho, y, N, k);
		c=compute_fm(u, d, rho, y, N, k) - DELTA_N_moins1 * psiprim(u, d, rho, y, N, N-2, true) - u[N-1]*u[N-1]/DELTA_N;
		nano = a >= 0 ? (a + sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a - sqrt(fabs(a*a - 4*b*c)));
		return y + nano;
	}
}
/*double* _middle_way_fm(double *u, double *d, double rho, double y, int N, int k) {
	double *res=(double*)malloc(sizeof(double)*3);
	double DELTA_k=d[k-1] - y;
	double DELTA_k_plus1=d[k+1] - y;
*/
	/* computing s */
//	res[1]=DELTA_k*DELTA_k*psiprim(u, d, rho, y, N, k, false);

	/* computing S */
//	res[1]=DELTA_k_plus1*DELTA_k_plus1*fiprim(u, d, rho, y, N, k, false);

	/* computing C */
//	if(k != N-1)
//		res[0]=compute_f(u, d, rho, y, N) - DELTA_k * psiprim(u, d, rho, y, N, k, false) - DELTA_k_plus1 * fiprim(u, d, rho, y, N, k, false);
//	else {/* case k=n-1 : last interval */
//		double DELTA_N = d[N-1] - y; /* N - 1 ???! */
//		double DELTA_N_moins1 = d[N-2] - y;
//		res[0]=compute_f(u, d, rho, y, N) - DELTA_N_moins1 * psiprim(u, d, rho, y, N, N-2, false) - u[N-1]*u[N-1]/DELTA_N;
//	}	
//	return res;
//}

/**
 * Fixed weight on f
 * @return c, s and S 
 */
double fixed_weight_f(double *u, double *d, double rho, double y, int N, int k) {
	double DELTA_k = d[k] - y;
	double DELTA_k_plus1 = d[k+1] - y;
	double a=compute_param_a(u, d, rho, y, N, k);
	double b=compute_param_b(u, d, rho, y, N, k);
	double c, nano;
	if(fabs(y - d[k]) < fabs(d[k+1] - y)) /* y is closer to d[k] */
		c=compute_f(u, d, rho, y, N) - DELTA_k_plus1 * compute_fprim(u, d, rho, y, N) - (u[k]*u[k] / (DELTA_k*DELTA_k)) * (d[k] - d[k+1]);
	else /* y is closer to d[k+1] */
		c=compute_f(u, d, rho, y, N) - DELTA_k * compute_fprim(u, d, rho, y, N) - (u[k+1]*u[k+1] / (DELTA_k_plus1*DELTA_k_plus1)) * (d[k+1] - d[k]);
	nano = a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
	return y+nano;
}
/*
double* fixed_weight_f(double *u, double *d, double rho, double y, int N, int k) {
	double *res=(double*)malloc(sizeof(double)*3);
	double DELTA_k = d[k] - y;
	double DELTA_k_plus1 = d[k+1] - y;
	double sum=0;

	if(fabs(y - d[k]) < fabs(d[k+1] - y)) {
		res[1]=u[k]*u[k];

		for(int j=0; j<N; j++) {
			if(j!=k && j!=k+1) sum+=(d[k+1]-y)*(d[k+1]-y)*u[j]*u[j]/((d[j]-y)*(d[j]-y));
		}
		res[2]=u[k+1]*u[k+1]+sum;

		res[0]=compute_f(u, d, rho, y, N) - DELTA_k_plus1 * compute_fprim(u, d, rho, y, N) - (u[k]*u[k] / (DELTA_k*DELTA_k)) * (d[k] - d[k+1]);
	}
	else {
		for(int j=0; j<N; j++) {
			if(j!=k && j!=k+1) sum+=(d[k]-y)*(d[k]-y)*u[j]*u[j]/((d[j]-y)*(d[j]-y));
		}
		res[2]=u[k]*u[k]+sum;

		res[1]=u[k+1]*u[k+1];

		res[0]=compute_f(u, d, rho, y, N) - DELTA_k * compute_fprim(u, d, rho, y, N) - (u[k+1]*u[k+1] / (DELTA_k_plus1*DELTA_k_plus1)) * (d[k+1] - d[k]);
	}

	return res;
}
 */

/**
 * Fixed weight on fm
 * @return c, s and S 
 */
double fixed_weight_fm(double *u, double *d, double rho, double y, int N, int k) {
	double DELTA_k = d[k] - y;
	double DELTA_k_plus1 = d[k+1] - y;
	double a=(DELTA_k + DELTA_k_plus1) * compute_fm(u, d, rho, y, N, k) - DELTA_k * DELTA_k_plus1 * compute_fmprim(u, d, rho, y, N, k);
	double b=DELTA_k * DELTA_k_plus1 * compute_fm(u, d, rho, y, N, k);
	double c, nano;
	if(fabs(y - d[k]) < fabs(d[k+1] - y)) /* y is closer to d[k-1] */
		c=compute_fm(u, d, rho, y, N, k) - DELTA_k_plus1 * compute_fmprim(u, d, rho, y, N, k) - (u[k]*u[k] / (DELTA_k*DELTA_k)) * (d[k] - d[k+1]);
	else /* y is closer to d[k+1] */
		c=compute_fm(u, d, rho, y, N, k) - DELTA_k * compute_fmprim(u, d, rho, y, N, k) - (u[k+1]*u[k+1] / (DELTA_k_plus1*DELTA_k_plus1)) * (d[k+1] - d[k]);
	nano = a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
	return y + nano;
}

/*double* _fixed_weight_fm(double *u, double *d, double rho, double y, int N, int k) {
	double *res=(double*)malloc(sizeof(double)*3);
	double DELTA_k = d[k-1] - y;
	double DELTA_k_plus1 = d[k+1] - y;
	double sum=0;

	if(fabs(y - d[k]) < fabs(d[k+1] - y)) { *//* y is closer to d[k] */
		/* computing s */
//		res[1]=u[k]*u[k];

		/* computing S */
/*		for(int j=0; j<N; j++) {
			if(j!=k && j!=k+1) sum+=(d[k+1]-y)*(d[k+1]-y)*u[j]*u[j]/((d[j]-y)*(d[j]-y));
		}
		res[2]=u[k+1]*u[k+1]+sum;
*/
		/* computing c */
//		res[0]=compute_f(u, d, rho, y, N) - DELTA_k_plus1 * compute_fprim(u, d, rho, y, N) - (u[k]*u[k] / (DELTA_k*DELTA_k)) * (d[k] - d[k+1]);
//	}
//	else { /* y is closer to d[k+1] */
		/* computing s */
/*		for(int j=0; j<N; j++) {
			if(j!=k && j!=k+1) sum+=(d[k]-y)*(d[k]-y)*u[j]*u[j]/((d[j]-y)*(d[j]-y));
		}
		res[2]=u[k]*u[k]+sum;
*/
		/* computing S */		
//		res[1]=u[k+1]*u[k+1];

		/* computing c */
//		res[0]=compute_f(u, d, rho, y, N) - DELTA_k * compute_fprim(u, d, rho, y, N) - (u[k+1]*u[k+1] / (DELTA_k_plus1*DELTA_k_plus1)) * (d[k+1] - d[k]);
//	}

//	return res;
//}

/**
 * Dichotomie to find a root of f beetwen d_k and d_k+1
 * with s -> param[1] and S -> param[2]
 */
double dich_Q(double *u, double *d, double rho, int N, int k, double* param, double x, double y, double fy, double deb, double fin) {
	double fx=1;
	int i=0;
	//	printf("x before %.6g ; ", x);
	while((fin-deb) > pow(10, -12)) {
		fx=fy + (x-y)*(param[1] / ((d[k-1] - y) - (x - y)) + u[k]*u[k]/((d[k] - y) - (x - y)) + param[2]/((d[k+1] - y) - (x - y)));
		if(fx*fy <= 0) fin=x;
		else deb=x;
		x=(fin+deb)/2.0;
		fy=fx;
		i++;
	}
	//	printf("x after %d tour %.6g & fx=%.6g\n", i, x, fx);
	return x;
}

results_t zero_finder_hybrid(double *u, double *d, double rho, int N, int k) {
	results_t res;
	int K=0;
	bool stopcriteria=false, fixed=true;
	double y=find_initial_guess(d, u, rho, N, k);
	double tho=0;
	double fpre, fnew, fK;
	res.iter=0;

	if(d[k] < y && y < (d[k]+d[k+1])/2) {
		K=k;
		tho=y-d[k];
	}
	else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
		K=k+1;
		tho=d[k+1]-y;
	}

	fpre=compute_f(u, d, rho, y, N); /* y=d[k]+tho -> fpre=f(d[k] + tho) */
	fK=fpre - u[K]*u[K] / (-tho);

	if(fK > 0)
		y=fixed_weight_f(u, d, rho, y, N, k); /* We interpolate f(x) */
	else {
		y=fixed_weight_fm(u, d, rho, y, N, k); /* We interpolate fk(x) */
/*		fst_y=y;
		fy=compute_f(u, d, rho, fst_y, N);
		s=compute_s(u, d, rho, y, N, k, true);
		S=compute_S(u, d, rho, y, N, k, true);*/
//		y=find_root_Q(u, d, rho, N, k, s, S, y, fst_y, fy, d[k], d[k+1]); /* Compute s and S */
	}

	fnew=compute_f(u, d, rho, y, N); /* y=d[k]+tho1 -> f1=f(d[k] + tho1) */
	if(fnew < 0 && fabs(fnew) > 0.1 * fabs(fpre)) /* First switch */
		fixed=false;

//	stopcriteria=stopcriteria_non_monotically(u, d, rho, y, N, k);

	while(!stopcriteria && res.iter<3) {
		if(d[k] < y && y < (d[k]+d[k+1])/2) {
			K=k;
			tho=y-d[k];
		}
		else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
			K=k+1;
			tho=d[k+1]-y;
		}
		fK=fnew - u[K]*u[K] / (-tho);
		if(fK > 0) { /* We interpolate f(x) */
			if(fixed)/* we use fixed weight */
				y=fixed_weight_f(u, d, rho, y, N, k);
			else /* we use middle way */
				y=middle_way_f(u, d, rho, y, N, k);
		} else { /* We interpolate fk(x) */
			if(fixed)/* we use fixed weight */
				y=fixed_weight_fm(u, d, rho, y, N, k);
			else /* we use middle way */
				y=middle_way_fm(u, d, rho, y, N, k);
/*			if(fst_y==-1) {
				fst_y=y;
				fy=compute_f(u, d, rho, fst_y, N);
			}
			s=compute_s(u, d, rho, y, N, k, fixed);
			S=compute_S(u, d, rho, y, N, k, fixed);
			y=find_root_Q(u, d, rho, N, k, s, S, y, fst_y, fy, d[k], d[k+1]); *//* Compute s and S */
		}
		fpre=fnew;
		fnew=compute_f(u, d, rho, y, N);
		if(fnew * fpre > 0 && fabs(fnew) > 0.1 * fabs(fpre)) /* New switch */
			fixed=!fixed;
		//		printf("%d fixed=%s\n", iter, fixed ? "Oui" : "False");
//		stopcriteria=stopcriteria_non_monotically(u, d, rho, y, N, k);
		res.iter++;
	}

	res.lambda=y;
	return res;
}
/*
results_t _zero_finder_hybrid(double *u, double *d, double rho, int N, int k) {
	results_t res;
	bool fixed=true;
	double tho=0, y=0, fst_y, fy;
	double *param; *//* param[0] -> c, param[1] -> s, param[2] -> S */
/*	double fpre, fnew, fk;
	int K;
	res.iter=0;

	tho=y-d[k];

	if(d[k] < y && y < (d[k]+d[k+1])/2) {
		K=k;
		tho=y-d[k];
	}
	else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
		K=k+1;
		tho=d[k+1]-y;
	}
	fst_y=find_initial_guess(d, u, rho, N, k);	
*/
	/* first iteration */
//	fy=compute_f(u, d, rho, fst_y, N); /* fst_y=d[k]+tho -> fpre=f(d[k] + tho) */
//	fk=fy - u[k]*u[k]/(-tho);

//	if(fk < 0) { /* Normal case */
//		y=fixed_weight_f(u, d, rho, fst_y, N, k); /* We interpolate f(x) */
//	}
//	else { /* dificult cases */
//		if(d[k] < fst_y && fst_y < (d[k]+d[k+1])/2) {
//			param=fixed_weight_fm(u, d, rho, fst_y, N, k); /* We interpolate fk(x) */
/*			y=dich_Q(u, d, rho, N, k, param, fst_y, fst_y, fy, d[k], (d[k]+d[k+1])/2);
		}
		else if((d[k]+d[k+1])/2 < fst_y && fst_y < d[k+1]) {*/
//			param=fixed_weight_fm(u, d, rho, fst_y, N, k+1); /* We interpolate fk+1(x) */
/*			y=dich_Q(u, d, rho, N, k, param, fst_y, fst_y, fy, (d[k]+d[k+1])/2, d[k+1]);
		}
	}
*/
	/* second iteration */
//	fpre=fy;
//	fnew=compute_f(u, d, rho, y, N); /* y=d[k]+tho1 -> f1=f(d[k] + tho1) */
//	if(fnew < 0 && fabs(fnew) > 0.1 * fabs(fpre)) /* First switch */
/*		fixed=false;
	tho=y-d[k];
	if(d[k] < y && y < (d[k]+d[k+1])/2) {
		K=k;
		tho=y-d[k];
	}
	else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
		K=k+1;
		tho=d[k+1]-y;
	}
	fk=fnew - u[k]*u[k]/(-tho);
	if(fk > 0) { *//* Normal case : We interpolate f(x) */
//		if(fixed)/* we use fixed weight */
//			y=fixed_weight_f(u, d, rho, y, N, k);
//		else /* we use middle way */
//			y=middle_way_f(u, d, rho, y, N, k);
//	}
//	else { /* dificult cases : We interpolate fk(x) */
//		if(d[k] < y && y < (d[k]+d[k+1])/2) {
//			if(fixed)/* we use fixed weight */
//				param=fixed_weight_fm(u, d, rho, y, N, k);
//			else /* we use middle way */
//				param=middle_way_fm(u, d, rho, y, N, k);
//			y=dich_Q(u, d, rho, N, k, param, y, y, fnew, d[k], (d[k]+d[k+1])/2);
//		}
//		else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
//			if(fixed)/* we use fixed weight */
//				param=fixed_weight_fm(u, d, rho, y, N, k);
//			else /* we use middle way */
/*				param=middle_way_fm(u, d, rho, y, N, k);
			y=dich_Q(u, d, rho, N, k, param, y, y, fnew, (d[k]+d[k+1])/2, d[k+1]);
		}
	}
*/
	/* Next iterations */
/*	while(!stopcriteria_non_monotically(u, d, rho, y, N, k) && res.iter<7 && !isnan(y)) {
		fpre=fnew;
		fnew=compute_f(u, d, rho, y, N);*/ /* y=d[k]+tho1 -> f1=f(d[k] + tho1) */
//		if(fnew * fpre > 0 && fabs(fnew) > 0.1 * fabs(fpre)) /* New switch */
//			fixed=!fixed;
//		tho=y-d[k];
		/*		if(d[k] < y && y < (d[k]+d[k+1])/2) {
				K=k;
				tho=y-d[k];
			}
			else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {
				K=k+1;
				tho=d[k+1]-y;
			}
		 */		//fk=fnew - u[k]*u[k] /(-tho);
//		 if(fk > 0) { /* Normal case : We interpolate f(x) */
//			 if(fixed)/* we use fixed weight */
//				 y=fixed_weight_f(u, d, rho, y, N, k);
//			 else /* we use middle way */
//				 y=middle_way_f(u, d, rho, y, N, k);
//		 }
//		 else { /* dificult cases : We interpolate fk(x) */
//			 if(d[k] < y && y < (d[k]+d[k+1])/2) {
//				 if(fixed)/* we use fixed weight */
//					 param=fixed_weight_fm(u, d, rho, y, N, k);
//				 else /* we use middle way */
/*					 param=middle_way_fm(u, d, rho, y, N, k);
				 y=dich_Q(u, d, rho, N, k, param, y, y, fnew, d[k], (d[k]+d[k+1])/2);
			 }
			 else if((d[k]+d[k+1])/2 < y && y < d[k+1]) {*/
//				 if(fixed)/* we use fixed weight */
//					 param=fixed_weight_fm(u, d, rho, y, N, k);
//				 else /* we use middle way */
//					 param=middle_way_fm(u, d, rho, y, N, k);
/*				 y=dich_Q(u, d, rho, N, k, param, y, y, fnew, (d[k]+d[k+1])/2, d[k+1]);
			 }
		 }
		 res.iter++;
	}

	res.lambda=y;
	return res;
}
*/
/**
 * Hybrid scheme for finding eigenvalues of a matrix
 */
results_t hybrid(double *u, double *d, double rho, int N) {
	results_t tmp, res;
	res.iter=0;
	double* lambdas = (double*)(calloc(sizeof(double), N));
#pragma omp parallel for private(tmp) schedule(runtime) ordered
	for(int k=0; k<N; k++) {
		tmp=zero_finder_hybrid(u, d, fabs(rho), N, k);
		res.iter+=tmp.iter;
		lambdas[k]=tmp.lambda;
	}
//	res.lambda=res.iter;
	res.iter/=N;
	res.lambdas=lambdas;
	return res;
}

int main(int argc, char **argv) {
	int ntests;
	int method;
	int N;
	if(argc == 4) {
		method=atoi(argv[1]);
		N=atoi(argv[2]);
		ntests=atoi(argv[3]);
	}
	else {
		fprintf(stderr, "Too few arguments. method, matrice length or number of test miss\n"
				"example : ./secular <0 | 1 : int> <matrice length : int> <ntest : int>\n"
				"          ./secular 0 1000 5\n"
				"Where 0 -> Gragg method\n"
				"      1 -> Hybrid scheme\n");
		exit(1);
	}

	if(N < 4 || ntests<1) exit(1);

	results_t results;
	double begin, end;
	double *d, *u;
	double rho;
	double betas[3]={pow(10, -3), pow(10, -6), pow(10, -10)};
	double w[4];
	char *s;
	if(!method) s="Gragg Method";
	else s="Hybrid Scheme";

	for(int i=1; i<=ntests; i++) {
		w[0]=w[3]=2;
		w[1]=w[2]=betas[i%4];
		rho=pow(-1, i%2+1)*pow(norm_2(w, 4), -2);
		u=(double*)malloc(sizeof(double)*N);
		d=(double*)malloc(sizeof(double)*(N+1));
		generate_U(N, u);
		generate_deltas(N+1, u, d, betas[i%4], rho);

		printf("u={");
		for(int j=0; j<N; j++) printf("%g ", u[j]);
		printf("}");

		printf("d={");
		for(int j=0; j<N; j++) printf("%g ", d[j]);
		printf("}");
		
		/* Beginning of chrono */
		begin = my_gettimeofday();
		if(!method) results=hybrid(u, d, rho, N);
		else results=gragg(u, d, rho, N);
		end = my_gettimeofday();
		print_test_result(results, s, i, u, d, rho, N, true);
		/* End of chrono */
		printf("Temps total de calcul : %g seconde(s) \n\n", end - begin);

		free(u);
		free(d);
	}
}
