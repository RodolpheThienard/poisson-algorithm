/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <lapack.h>

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

    f64 * tmp_vec, *sol;

    f64 temp, forward_error, backward_error;

    NRHS = 1;
    nbpoints = 10;
    la = nbpoints - 2;
    T0 = -5.0;
    T1 = 5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS = (f64 *)malloc(sizeof(f64) * la);
    EX_SOL = (f64 *)malloc(sizeof(f64) * la);
    X = (f64 *)malloc(sizeof(f64) * la);

    sol = malloc (sizeof (double) * la);
    tmp_vec = malloc (sizeof (double) * la);

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

    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    printf("Solution with LAPACK\n");
    /* LU Factorization */
    info = 0;
    ipiv = (i32 *)calloc(la, sizeof(i32));
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* LU for tridiagonal matrix  (can replace dgbtrf_) */
    // ierr= dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");

    /* Solution (Triangular) */
    if (info == 0)
    {
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info != 0)
        {
            printf("\n INFO DGBTRS = %d\n", info);
        }
    }
    else
    {
        printf("\n INFO = %d\n", info);
    }

    /* It can also be solved with dgbsv */
    // TODO : use dgbsv
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

    write_xy(RHS, X, &la, "SOL.dat");

    /* Relative forward and backward error */
    // TODO : Compute relative norm of the residual

    /* Forward error :    ||x-xø|| / ||x|| */

    memcpy(tmp_vec, RHS, la * sizeof (f64));

    cblas_daxpy(la, -1, EX_SOL, 1, tmp_vec, 1);

    f64 norm_x_xø = cblas_dnrm2(la, tmp_vec, 1);
    
    f64 norm_x = cblas_dnrm2(la, EX_SOL, 1);

    forward_error = norm_x_xø / norm_x;


    /* Backward error :    ||b-Axø|| / ||b|| */

    memcpy(tmp_vec, RHS, la * sizeof (f64));

    cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, -1.0, AB, lab, sol, 1, 1.0, tmp_vec, 1);

    f64 tmp = cblas_dnrm2(la, tmp_vec, 1);

    backward_error = tmp / cblas_dnrm2(la, RHS, 1);


    printf("\nThe relative forward error is relres = %e\n", forward_error );
    printf("\nThe relative backward error is relres = %e\n", backward_error );

    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    printf("\n\n--------- End -----------\n");
}
