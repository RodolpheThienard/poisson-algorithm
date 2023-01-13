/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <cblas.h>
#include <lapack.h>
#include <math.h>
#include <string.h>

//
void set_GB_operator_colMajor_poisson1D(f64 *AB, i32 *lab, i32 *la, i32 *kv)
{
	i32 ii, jj, kk;
	for (jj = 0; jj < (*la); jj++) {
		kk = jj * (*lab);
		if (*kv >= 0) {
			for (ii = 0; ii < *kv; ii++) {
				AB[kk + ii] = 0.0;
			}
		}
		AB[kk + *kv] = -1.0;
		AB[kk + *kv + 1] = 2.0;
		AB[kk + *kv + 2] = -1.0;
	}
	AB[0] = 0.0;
	if (*kv == 1) {
		AB[1] = 0;
	}

	AB[(*lab) * (*la) - 1] = 0.0;
}

//
void set_GB_operator_colMajor_poisson1D_Id(f64 *AB, i32 *lab, i32 *la, i32 *kv)
{
	i32 ii, jj, kk;
	for (jj = 0; jj < (*la); jj++) {
		kk = jj * (*lab);
		if (*kv >= 0) {
			for (ii = 0; ii < *kv; ii++) {
				AB[kk + ii] = 0.0;
			}
		}
		AB[kk + *kv] = 0.0;
		AB[kk + *kv + 1] = 1.0;
		AB[kk + *kv + 2] = 0.0;
	}
	AB[1] = 0.0;
	AB[(*lab) * (*la) - 1] = 0.0;
}

//
void set_dense_RHS_DBC_1D(f64 *RHS, i32 *la, f64 *BC0, f64 *BC1)
{
	RHS[0] = *BC0;
	RHS[(*la) - 1] = *BC1;
	for (i32 jj = 1; jj < (*la) - 1; jj++) {
		RHS[jj] = 0.0;
	}
}

//
void set_analytical_solution_DBC_1D(f64 *EX_SOL, f64 *X, i32 *la, f64 *BC0,
				    f64 *BC1)
{
	for (i32 i = 0; i < *la; i++, EX_SOL++, X++) {
		*EX_SOL = *BC0 + *X * (*BC1 - *BC0);
	}
}

//
void set_grid_points_1D(f64 *x, i32 *la)
{
	i32 i;
	f64 h = 1.0 / (1.0 * (*la + 1.0));
	for (i = 0; i < (*la); i++)
		x[i] = (i + 1) * h;
}

void write_GB_operator_rowMajor_poisson1D(f64 *AB, i32 *lab, i32 *la,
					  ascii *filename)
{
	FILE *file;
	i32 ii, jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (ii = 0; ii < (*lab); ii++) {
			for (jj = 0; jj < (*la); jj++) {
				fprintf(file, "%lf\t", AB[ii * (*la) + jj]);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

void write_GB_operator_colMajor_poisson1D(f64 *AB, i32 *lab, i32 *la,
					  ascii *filename)
{
	FILE *file;
	i32 ii, jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (ii = 0; ii < (*la); ii++) {
			for (jj = 0; jj < (*lab); jj++) {
				fprintf(file, "%lf\t", AB[ii * (*lab) + jj]);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

void write_GB2AIJ_operator_poisson1D(f64 *AB, i32 *la, ascii *filename)
{
	FILE *file;
	i32 jj;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 1; jj < (*la); jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj, jj + 1,
				AB[(*la) + jj]);
		}
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj + 1, jj + 1,
				AB[2 * (*la) + jj]);
		}
		for (jj = 0; jj < (*la) - 1; jj++) {
			fprintf(file, "%d\t%d\t%lf\n", jj + 2, jj + 1,
				AB[3 * (*la) + jj]);
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

void write_vec(f64 *vec, i32 *la, ascii *filename)
{
	i32 jj;
	FILE *file;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%lf\n", vec[jj]);
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

void write_value(f64 *vec, i32 *la, ascii *filename)
{
	i32 jj;
	FILE *file;
	file = fopen(filename, "a");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%e\n", vec[jj]);
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

void write_xy(f64 *vec, f64 *x, i32 *la, ascii *filename)
{
	i32 jj;
	FILE *file;
	file = fopen(filename, "w");
	// Numbering from 1 to la
	if (file != NULL) {
		for (jj = 0; jj < (*la); jj++) {
			fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
		}
		fclose(file);
	} else {
		perror(filename);
	}
}

i32 indexABCol(i32 i, i32 j, i32 *lab)
{
	return j * (*lab) + i + 1;
}
i32 dgbtrftridiag(i32 *la, i32 *n, i32 *kl, i32 *ku, f64 *AB, i32 *lab,
		  i32 *ipiv, i32 *info)
{
	for (u32 i = 1; i < *la; i++) {
		u32 a_1 = (*lab) * (i - 1) + 1;
		u32 a = (*lab) * (i) + 1;
		u32 b_1 = (*lab) * (i - 1) + 2;
		u32 c_1 = (*lab) * (i);

		ipiv[i - 1] = i;

		/* b(i-1) = b(i-1)/a(i-1) */
		AB[(*kl) + b_1] /= AB[(*kl) + a_1];

		/* a(i) = a(i) - b(i-1) / a(i-1) * c(i-1) */
		AB[(*kl) + a] -= AB[(*kl) + b_1] * AB[(*kl) + c_1];
	}

	return 0; // Value never used
}

void eig_poisson1D(f64 *eigval, i32 *la)
{
}

f64 eigmax_poisson1D(i32 *la)
{
	f64 h = 1 / (*la) + 1;
	f64 sin = (M_PI * h * (*la)) / 2;
	return 4 * sin * sin;
}

f64 eigmin_poisson1D(i32 *la)
{
	f64 h = 1 / (*la) + 1;
	f64 sin = (M_PI * h) / 2;
	return 4 * sin * sin;
}

f64 richardson_alpha_opt(i32 *la)
{
	return 0.5;
	return 2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(f64 *AB, f64 *RHS, f64 *X, f64 *alpha_rich, i32 *lab,
		      i32 *la, i32 *ku, i32 *kl, f64 *tol, i32 *maxit,
		      f64 *resvec, i32 *nbite)
{
	/*
    r^0 = b - Ax^0
    while ||r^k+1|| > epsilon do:
      x^k+1 = x^k + alpha(b-Ax^k)
  */

	f64 norm = 0;
	f64 *b = (f64 *)calloc(*la, sizeof(f64));
	f64 norm_b = cblas_dnrm2(*la, RHS, 1);
	cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB,
		    *lab, X, 1, 1, b, 1);

	resvec[*nbite] = cblas_dnrm2(*la, b, 1) / norm_b;

	// x^k+1 = x^k +alpha(b-Ax^k)
	cblas_daxpy(*la, *alpha_rich, b, 1, X, 1);

	do {
		memcpy(b, RHS, *la * sizeof(f64));

		// b - Ax^k
		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku,
			    -1.0, AB, *lab, X, 1, 1, b, 1);
		resvec[*nbite] = cblas_dnrm2(*la, b, 1) / norm_b;

		// x^k+1 = x^k +alpha(b-Ax^k)
		cblas_daxpy(*la, *alpha_rich, b, 1, X, 1);
		(*nbite)++;
	} while (resvec[*nbite - 1] > *tol && *nbite - 1 < *maxit);
	(*nbite)--;
}

void extract_MB_jacobi_tridiag(f64 *AB, f64 *MB, i32 *lab, i32 *la, i32 *ku,
			       i32 *kl, i32 *kv)
{
	/* M = D */
	for (i32 i = 0; i < *la; i++) {
		MB[*lab * i + 1] = AB[*lab * i + 1];
	}
}

void extract_MB_gauss_seidel_tridiag(f64 *AB, f64 *MB, i32 *lab, i32 *la,
				     i32 *ku, i32 *kl, i32 *kv)
{
	/* M = D - E */

	for (i32 i = 0; i < *la; i++) {
		MB[*lab * i + 1] = AB[*lab * i + 1];
		MB[*lab * i + 2] = AB[*lab * i + 2];
	}
}

void richardson_MB(f64 *AB, f64 *RHS, f64 *X, f64 *MB, i32 *lab, i32 *la,
		   i32 *ku, i32 *kl, f64 *tol, i32 *maxit, f64 *resvec,
		   i32 *nbite)
{
	f64 *B = malloc(*la * sizeof(f64));
	f64 norme_B = cblas_dnrm2(*la, RHS, 1);
	i32 *ipiv = malloc(*la * sizeof(i32));
	i32 info = 0;
	i32 NRHS = 1;
	i32 ku_ = *ku - 1;

	dgbtrf_(la, la, kl, &ku_, MB, lab, ipiv, &info);
	for ((*nbite) = 0; (*nbite) < *maxit; (*nbite)++) {
		if (resvec[*nbite] <= *tol)
			break;
		memcpy(B, RHS, *la * sizeof(f64));

		cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku,
			    -1.0, AB, *lab, X, 1, 1.0, B, 1);

		resvec[*nbite] = cblas_dnrm2(*la, B, 1) / norme_B;

		dgbtrs_("N", la, kl, &ku_, &NRHS, MB, lab, ipiv, B, la, &info,
			1);

		cblas_daxpy(*la, 1, B, 1, X, 1);
	}
	free(B);
	free(ipiv);
}