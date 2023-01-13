/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <bits/time.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
	i32 ierr;
	i32 jj;
	i32 nbpoints, la;
	i32 ku, kl, lab, kv;
	i32 *ipiv;
	i32 info;
	i32 NRHS;
	f64 T0, T1;
	f64 *RHS, *SOL, *EX_SOL, *X;
	f64 *AB;
	f64 *MB;
	f64 *tmp_vec;
	f64 temp, relres;

	f64 opt_alpha, forward_error;

	/* Size of the problem */
	NRHS = 1;
	nbpoints = atoi(argv[1]);
	la = nbpoints - 2;

	/* Dirichlet Boundary conditions */
	T0 = 5.0;
	T1 = 20.0;

	// printf("--------- Poisson 1D ---------\n\n");
	RHS = (f64 *)malloc(sizeof(f64) * la);
	SOL = (f64 *)calloc(la, sizeof(f64));
	EX_SOL = (f64 *)malloc(sizeof(f64) * la);
	X = (f64 *)malloc(sizeof(f64) * la);
	tmp_vec = malloc(sizeof(f64) * la);
	/* Setup the Poisson 1D problem */
	/* General Band Storage */
	set_grid_points_1D(X, &la);
	set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
	set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

	write_vec(RHS, &la, "RHS.dat");
	write_vec(EX_SOL, &la, "EX_SOL.dat");
	write_vec(X, &la, "X_grid.dat");

	kv = 0;
	ku = 1;
	kl = 1;
	lab = kv + kl + ku + 1;

	AB = (f64 *)malloc(sizeof(f64) * lab * la);
	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

	/* uncomment the following to check matrix A */
	// write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

	/********************************************/
	/* Solution (Richardson with optimal alpha) */

	/* Computation of optimum alpha */
	opt_alpha = atof(argv[2]); // richardson_alpha_opt(&la);
	// printf("Optimal alpha for simple Richardson iteration is : %lf",
	//        opt_alpha);

	/* Solve */
	f64 tol = 1e-3;
	i32 maxit = 1000;
	f64 *resvec;
	i32 nbite = 0;
	struct timespec t1, t2;
	f64 elapsed_alpha, elapsed_jacobi, elapsed_gauss;
	i32 value = 1;
	resvec = (f64 *)calloc(maxit, sizeof(f64));

	/* Solve with Richardson alpha */
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol,
			 &maxit, resvec, &nbite);
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	elapsed_alpha = (((f64)(t2.tv_sec - t1.tv_sec)) +
			 (t2.tv_nsec - t1.tv_nsec) * 1e-9L);

	/* Richardson General Tridiag */

	/* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
	kv = 1;
	ku = 1;
	kl = 1;
	MB = (f64 *)malloc(sizeof(f64) * (lab)*la);
	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

	extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
	write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "JACOBI.dat");

	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit,
		      resvec, &nbite);
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	elapsed_jacobi = (((f64)(t2.tv_sec - t1.tv_sec)) +
			  (t2.tv_nsec - t1.tv_nsec) * 1e-9L);

	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
	extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
	write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "GAUSS-SEIDEL.dat");

	/* Solve with General Richardson */
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit,
		      resvec, &nbite);
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	elapsed_gauss = (((f64)(t2.tv_sec - t1.tv_sec)) +
			 (t2.tv_nsec - t1.tv_nsec) * 1e-9L);

	/* Forward error :    ||x-xø|| / ||x|| */

	memcpy(tmp_vec, RHS, la * sizeof(f64));

	cblas_daxpy(la, -1, EX_SOL, 1, tmp_vec, 1);

	f64 norm_x_xø = cblas_dnrm2(la, tmp_vec, 1);

	f64 norm_x = cblas_dnrm2(la, EX_SOL, 1);

	forward_error = norm_x_xø / norm_x;

	// printf("\n%lf\n", forward_error);
	/* Write solution */
	write_vec(SOL, &la, "SOL.dat");

	/* Write convergence history */
	write_vec(resvec, &nbite, "RESVEC.dat");
	printf("%lf : %lf\n", opt_alpha, elapsed_alpha);
	// printf("%d : alpha : %lf\n \
	// 				%d : jacobi: %lf\n \
	// 				%d : gauss : %lf\n", nbpoints, elapsed_alpha, nbpoints, elapsed_jacobi, nbpoints, elapsed_gauss);

	free(resvec);
	free(RHS);
	free(SOL);
	free(EX_SOL);
	free(X);
	free(AB);
	free(MB);
} // printf("\n\n--------- End -----------\n");
