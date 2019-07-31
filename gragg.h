#ifndef _GRAGG_H_
#define _GRAGG_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* les fonctions math */
#include <time.h>   /* chronometrage */
#include <sys/time.h>
#include <stdbool.h> /* Pour les booléens */
#include <omp.h>

typedef struct results_s {
  int iter;
  double lambda;
  double *lambdas;
} results_t;

double my_gettimeofday();
double min(double a, double b);
double norm_2(double *v, int n);
void normalize(int n, double *v);
void generate_U(int n, double* u);
double compute_complete_delta(double delta, double *u, double rho, int N);
void generate_deltas(int n, double* d, double *u, double beta, double rho);
void print_test_result(results_t res, char* s, int numtest, double *u, double *d, double rho, int N, bool pr);
double compute_f(double *u, double *d, double rho, double lambda_j, int N);
double compute_fprim(double *u, double *d, double rho, double lambda_j, int N);
double compute_fsnd(double *u, double *d, double rho, double lambda_j, int N);
double g(double *d, double *u, double rho, double x, int N, int k);
double h(double *d, double *u, double rho, double x, int N, int k);
double find_initial_guess(double *d, double *u, double rho, int N, int k);
double compute_gamma(double *u, double *d, double lambda_j, int N, int k);
double compute_omega0(double *u, double *d, double lambda_j, int N, int k);
double compute_omega1(double *u, double *d, double rho, double lambda_j, int N, int k);
double compute_param_a(double *u, double *d, double rho, double lambda_j, int N, int k);
double compute_param_b(double *u, double *d, double rho, double lambda_j, int N, int k);
double find_correction(double a, double b, double c);
double find_next_lambda(double a, double b, double c, double lambda_j);
bool stopcriteria_non_monotically(double *u, double *d, double rho, double lambda_j, int N, int k);
results_t zero_finder_gragg(double *u, double *d, double rho, int N, int k);
results_t gragg(double *u, double *d, double rho, int N);

#endif
