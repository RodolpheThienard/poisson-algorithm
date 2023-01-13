/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc, char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
	i32 ierr;
	i32 jj;
	i32 nbpoints, la;
	i32 ku, kl, kv, lab;
	i32 *ipiv;
	i32 info;
	i32 NRHS;
	f64 T0, T1;
	f64 *RHS, *EX_SOL, *X;
	f64 **AAB;
	f64 *AB;

	f64 *tmp_vec, *sol;

	f64 temp, forward_error, backward_error;

	struct timespec t1, t2;
	f64 elapsed_dgbtrs, elapsed_dgbtrf, elapsed_dgbsv,
		elapsed_dgbtrftridiag;

	elapsed_dgbsv = 0.0;
	elapsed_dgbtrs = 0.0;
	elapsed_dgbtrf = 0.0;
	elapsed_dgbtrftridiag = 0.0;

	u32 array_size = 10;
	if (argc > 1) {
		array_size = atoi(argv[1]);
	}

	NRHS = 1;
	nbpoints = array_size;
	//    printf("nb ptn : %d", nbpoints);
	la = nbpoints - 2;
	T0 = -5.0;
	T1 = 5.0;

	//    printf("--------- Poisson 1D ---------\n\n");
	RHS = (f64 *)malloc(sizeof(f64) * la);
	EX_SOL = (f64 *)malloc(sizeof(f64) * la);
	X = (f64 *)malloc(sizeof(f64) * la);

	sol = malloc(sizeof(double) * la);
	tmp_vec = malloc(sizeof(double) * la);

	// TODO : you have to implement those functions
	set_grid_points_1D(X, &la);
	set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
	set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

	write_vec(RHS, &la, "RHS.dat");
	write_vec(EX_SOL, &la, "EX_SOL.dat");
	write_vec(X, &la, "X_grid.dat");

	kv = 1;
	ku = 1;
	kl = 1;
	lab = kv + kl + ku + 1;

	AB = (f64 *)malloc(sizeof(f64) * lab * la);

	set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
	cblas_dgbmv(CblasColMajor, CblasTrans, la, la, kl, ku, 1, AB + 1, lab,
		    EX_SOL, 1, 0, RHS, 1);
	write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

	//    printf("Solution with LAPACK\n");
	/* LU Factorization */

	// f64 *ABB = (f64 *)malloc(sizeof(f64) * lab * la);
	// set_GB_operator_colMajor_poisson1D(ABB, &lab, &la, &kv);

	info = 0;
	ipiv = (i32 *)calloc(la, sizeof(i32));
	// clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	// dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
	// clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	// elapsed_dgbtrf = (((f64)(t2.tv_sec - t1.tv_sec)) +
	//  (t2.tv_nsec - t1.tv_nsec) * 10e-9L);

	// set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
	// cblas_dgbmv(CblasColMajor, CblasTrans, la, la, kl, ku, 1, AB + 1, lab, EX_SOL, 1, 0, RHS, 1);
	/* LU for tridiagonal matrix  (can replace dgbtrf_) */
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
	clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	elapsed_dgbtrftridiag = (((f64)(t2.tv_sec - t1.tv_sec)) +
				 (t2.tv_nsec - t1.tv_nsec) * 10e-9L);

	write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

	/* Solution (Triangular) */
	if (info == 0) {
		clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
		dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la,
			&info, 1);
		clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
		elapsed_dgbtrs = (((f64)(t2.tv_sec - t1.tv_sec)) +
				  (t2.tv_nsec - t1.tv_nsec) * 10e-9L);
		if (info != 0) {
			printf("\n INFO DGBTRS = %d\n", info);
		}
	} else {
		printf("\n INFO = %d\n", info);
	}

	// set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
	// cblas_dgbmv(CblasColMajor, CblasTrans, la, la, kl, ku, 1, AB + 1, lab, EX_SOL, 1, 0, RHS, 1);
	/* It can also be solved with dgbsv */
	// TODO : use dgbsv
	// clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	// dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
	// clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
	// elapsed_dgbsv = (((f64)(t2.tv_sec - t1.tv_sec)) +
	// 		 (t2.tv_nsec - t1.tv_nsec) * 10e-9L);
	// write_xy(RHS, X, &la, "SOL.dat");

	/* Relative forward and backward error */
	// TODO : Compute relative norm of the residual

	/* Forward error :    ||x-xø|| / ||x|| */

	f64 norm_x = cblas_dnrm2(la, EX_SOL, 1);
	cblas_daxpy(la, -1, EX_SOL, 1, RHS, 1);

	f64 norm_x_xø = cblas_dnrm2(la, RHS, 1);

	forward_error = norm_x_xø / norm_x;
	i32 err_i = 1;
	// write_value(&forward_error, &err_i, "dgbsv_err.dat");
	//    printf("\nThe relative forward error is relres = %e\n", forward_error );
	printf("%d : dgbtrf : %lf \n \
            %d : dgbtrs : %lf \n \
            %d : dgbsv : %lf\n \
            %d : dgbtrftridiag : %lf \n \
        %d : forward_error : %e\n",
	       nbpoints, elapsed_dgbtrf, nbpoints, elapsed_dgbtrs, nbpoints,
	       elapsed_dgbsv, nbpoints, elapsed_dgbtrftridiag, nbpoints,
	       forward_error);
	free(RHS);
	free(EX_SOL);
	free(X);
	free(AB);
	//    printf("\n\n--------- End -----------\n");
}
