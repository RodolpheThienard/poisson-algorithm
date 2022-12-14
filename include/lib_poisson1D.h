/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */ 
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "types.h"
#include <string.h>

void set_GB_operator_colMajor_poisson1D(f64* AB, i32* lab, i32 *la, i32 *kv);
void set_GB_operator_colMajor_poisson1D_Id(f64* AB, i32* lab, i32 *la, i32 *kv);
void set_dense_RHS_DBC_1D(f64* RHS, i32* la, f64* BC0, f64* BC1);
void set_analytical_solution_DBC_1D(f64* EX_SOL, f64* X, i32* la, f64* BC0, f64* BC1);
void set_grid_points_1D(f64* x, i32* la);
void write_GB2AIJ_operator_poisson1D(f64* AB, i32 *la, ascii* filename);
void write_GB_operator_rowMajor_poisson1D(f64* AB, i32* lab, i32 *la, ascii* filename);
void write_GB_operator_colMajor_poisson1D(f64* AB, i32* lab, i32* la, ascii* filename);
void write_vec(f64* vec, i32* la, ascii* filename);
void write_xy(f64* vec, f64* x, i32* la, ascii* filename);
i32 indexABCol(i32 i, i32 j, i32 *lab);
i32 dgbtrftridiag(i32 *la, i32 *n, i32 *kl, i32 *ku, f64 *AB, i32 *lab, i32 *ipiv, i32 *info);
