/*
 *  snoextras.c
 *  
 *
 *  Created by W. Ross Morrow on 3/19/12.
 *  Copyright 2012 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "snoextras.h"

void lintest();
void quadtest();
void nonlintest();
void mixtest();

int main(int argc, char *argv[])
{
	lintest();
	quadtest();
	nonlintest();
	mixtest();
	
	return 0;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Linear function test
 * 
 * Has F(x) = Ax for a sparse matrix A
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int lintest_fun_(integer		*Status,	// SNOPT status code
				 integer		*N,			// number of variables
				 doublereal		*x,			// current variable values
				 integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
				 integer		*FN,		// length of the vector of objective and constraint values
				 doublereal		*F,			// values (to be calculated) for objective and constraint values
				 integer		*needG,     // 0 if G(x) not needed, > 0 if it is
				 integer		*Gnnz,		// length of arrays iGvar and jGfun
				 doublereal		*Gdata,		// derivative values (MMF format)
				 char			*cu,		// character workspace
				 integer		*lencu,		// length of character workspace
				 integer		*iu,		// integer workspace
				 integer		*leniu,		// length of integer workspace
				 doublereal		*ru,		// double workspace
				 integer		*lenru )	// length of double workspace
{
	// F[0] =	x[0] +  x[1]
	// F[1] =			x[1]
	// F[2] =					x[2]
	// F[3] =	x[0]
	
	if( needF[0] > 0 ) {
		F[0] = x[0] + x[1];
		F[1] = x[1];
		F[2] = x[2];
		F[3] = x[0];
	}
	
	if( needG[0] > 0 ) {
		Gdata[0] = 1.0;
		Gdata[1] = 1.0;
		Gdata[2] = 1.0;
		Gdata[3] = 1.0;
		Gdata[4] = 1.0;
	}
	
	return 0;
}

void lintest()
{
	integer N, M, FN;
	double *x, *xLoBnds, *xUpBnds;
	integer Gnnz, *Grows, *Gcols;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * PROBLEM DATA  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	N = 3;
	M = 3;
	FN = N*M + 1;
	
	Gnnz = 5; // 5 out of twelve nonzeros in Jacobian
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * ALLOCATE PROBLEM VARIABLES  * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
	
	// Jacobian of the nonlinear part of objective and constraints
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// x point
	x[0] = 1.0;
	x[1] = 2.0;
	x[2] = 3.0;
	
	// x bound(s)
	xLoBnds[0] = -1.0e20; xUpBnds[0] =  1.0e20;
	xLoBnds[1] = -1.0e20; xUpBnds[1] =  1.0e20;
	xLoBnds[2] = -1.0e20; xUpBnds[2] =  1.0e20;
	
	// Constraint Jacobian sparsity pattern using FORTRAN-style indices, 
	// and a non-standard pattern (i.e., not row- or column-major)
	// SNOPT shouldn't care about the order of the appearance of derivatives
	Grows[0] = 1; Gcols[0] = 1;
	Grows[1] = 1; Gcols[1] = 2;
	Grows[2] = 2; Gcols[2] = 2;
	Grows[3] = 3; Gcols[3] = 3;
	Grows[4] = 4; Gcols[4] = 1;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// check jacobian 
	snopta_eval_G_check(N, FN,
						lintest_fun_,
						Gnnz, Grows, Gcols, 
						x, xLoBnds, xUpBnds,
						NULL, NULL, 
						NULL, NULL,
						NULL, NULL);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
	free(xLoBnds);
	free(xUpBnds);
	free(Grows);
	free(Gcols);
	
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Quadratic function test
 * 
 *		F(1) = x(1)^2
 *		F(2) = x(2)^2
 *		F(3) = x(3)^2
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int quadtest_fun_(integer		*Status,	// SNOPT status code
				 integer		*N,			// number of variables
				 doublereal		*x,			// current variable values
				 integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
				 integer		*FN,		// length of the vector of objective and constraint values
				 doublereal		*F,			// values (to be calculated) for objective and constraint values
				 integer		*needG,     // 0 if G(x) not needed, > 0 if it is
				 integer		*Gnnz,		// length of arrays iGvar and jGfun
				 doublereal		*Gdata,		// derivative values (MMF format)
				 char			*cu,		// character workspace
				 integer		*lencu,		// length of character workspace
				 integer		*iu,		// integer workspace
				 integer		*leniu,		// length of integer workspace
				 doublereal		*ru,		// double workspace
				 integer		*lenru )	// length of double workspace
{
	
	if( needF[0] > 0 ) {
		F[0] = x[0] * x[0];
		F[1] = x[1] * x[1];
		F[2] = x[2] * x[2];
	}
	
	if( needG[0] > 0 ) {
		Gdata[0] = 2.0 * x[0];
		Gdata[1] = 2.0 * x[1];
		Gdata[2] = 2.0 * x[2];
	}
	
	return 0;
}

void quadtest()
{
	integer N, M, FN;
	double *x, *xLoBnds, *xUpBnds;
	integer Gnnz, *Grows, *Gcols;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * PROBLEM DATA  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	N  = 3;
	M  = 2;
	FN = 3;
	
	Gnnz = 3; // 3 out of 9 nonzeros in Jacobian
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * ALLOCATE PROBLEM VARIABLES  * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
	
	// Jacobian of the nonlinear part of objective and constraints
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// x point
	x[0] = 1.0;
	x[1] = 2.0;
	x[2] = 3.0;
	
	// x bound(s)
	xLoBnds[0] = -1.0e20; xUpBnds[0] =  1.0e20;
	xLoBnds[1] = -1.0e20; xUpBnds[1] =  1.0e20;
	xLoBnds[2] = -1.0e20; xUpBnds[2] =  1.0e20;
	
	// Constraint Jacobian sparsity pattern using FORTRAN-style indices, 
	// and a non-standard pattern (i.e., not row- or column-major)
	// SNOPT shouldn't care about the order of the appearance of derivatives
	Grows[0] = 1; Gcols[0] = 1;
	Grows[1] = 2; Gcols[1] = 2;
	Grows[2] = 3; Gcols[2] = 3;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// check jacobian 
	snopta_eval_G_check(N, FN,
						quadtest_fun_,
						Gnnz, Grows, Gcols, 
						x, xLoBnds, xUpBnds,
						NULL, NULL, 
						NULL, NULL,
						NULL, NULL);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
	free(xLoBnds);
	free(xUpBnds);
	free(Grows);
	free(Gcols);
	
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Nonlinear function test
 * 
 *		F(1) = x(1)^3 + exp( x(2) / 4.0 )
 *		F(2) = log( x(1) ) * x(2)
 *		F(3) = sin( x(3) )
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int nonlintest_fun_(integer		*Status,	// SNOPT status code
				  integer		*N,			// number of variables
				  doublereal		*x,			// current variable values
				  integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
				  integer		*FN,		// length of the vector of objective and constraint values
				  doublereal		*F,			// values (to be calculated) for objective and constraint values
				  integer		*needG,     // 0 if G(x) not needed, > 0 if it is
				  integer		*Gnnz,		// length of arrays iGvar and jGfun
				  doublereal		*Gdata,		// derivative values (MMF format)
				  char			*cu,		// character workspace
				  integer		*lencu,		// length of character workspace
				  integer		*iu,		// integer workspace
				  integer		*leniu,		// length of integer workspace
				  doublereal		*ru,		// double workspace
				  integer		*lenru )	// length of double workspace
{
	
	if( needF[0] > 0 ) {
		F[0] = x[0] * x[0] * x[0] + exp( x[1] / 4.0 );
		F[1] = log( x[0] ) * x[1];
		F[2] = sin( x[2] );
	}
	
	if( needG[0] > 0 ) {
		Gdata[0] = 3.0 * x[0] * x[0];
		Gdata[1] = 0.25 * exp( x[1] / 4.0 );
		Gdata[2] = x[1] / x[0];
		Gdata[3] = log( x[0] );
		Gdata[4] = cos( x[2] );
	}
	
	return 0;
}

void nonlintest()
{
	integer N, M, FN;
	double *x, *xLoBnds, *xUpBnds;
	integer Gnnz, *Grows, *Gcols;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * PROBLEM DATA  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	N  = 3;
	M  = 2;
	FN = 3;
	
	Gnnz = 5; // 5 out of 9 nonzeros in Jacobian
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * ALLOCATE PROBLEM VARIABLES  * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
	
	// Jacobian of the nonlinear part of objective and constraints
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// x point (x[0] cannot be zero for differentiability)
	x[0] = 1.0;
	x[1] = 2.0;
	x[2] = 3.0;
	
	// x bound(s)
	xLoBnds[0] = -1.0e20; xUpBnds[0] =  1.0e20;
	xLoBnds[1] = -1.0e20; xUpBnds[1] =  1.0e20;
	xLoBnds[2] = -1.0e20; xUpBnds[2] =  1.0e20;
	
	// Constraint Jacobian sparsity pattern using FORTRAN-style indices, 
	// and a non-standard pattern (i.e., not row- or column-major)
	// SNOPT shouldn't care about the order of the appearance of derivatives
	Grows[0] = 1; Gcols[0] = 1;
	Grows[1] = 1; Gcols[1] = 2;
	Grows[2] = 2; Gcols[2] = 1;
	Grows[3] = 2; Gcols[3] = 2;
	Grows[4] = 3; Gcols[4] = 3;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// check jacobian 
	snopta_eval_G_check(N, FN,
						nonlintest_fun_,
						Gnnz, Grows, Gcols, 
						x, xLoBnds, xUpBnds,
						NULL, NULL, 
						NULL, NULL,
						NULL, NULL);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
	free(xLoBnds);
	free(xUpBnds);
	free(Grows);
	free(Gcols);
	
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Nonlinear function test
 * 
 *		F(1) = x(1)^3 + exp( x(2) / 4.0 )
 *		F(2) = log( x(1) ) * x(2)
 *		F(3) = sin( x(3) )
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int mixtest_fun_(integer		*Status,	// SNOPT status code
					integer		*N,			// number of variables
					doublereal		*x,			// current variable values
					integer		*needF,		// 0 if f(x) is not needed, > 0 if it is
					integer		*FN,		// length of the vector of objective and constraint values
					doublereal		*F,			// values (to be calculated) for objective and constraint values
					integer		*needG,     // 0 if G(x) not needed, > 0 if it is
					integer		*Gnnz,		// length of arrays iGvar and jGfun
					doublereal		*Gdata,		// derivative values (MMF format)
					char			*cu,		// character workspace
					integer		*lencu,		// length of character workspace
					integer		*iu,		// integer workspace
					integer		*leniu,		// length of integer workspace
					doublereal		*ru,		// double workspace
					integer		*lenru )	// length of double workspace
{
	
	if( needF[0] > 0 ) {
		F[0] = x[0] * x[0] * x[0] + exp( x[1] / 4.0 );
		F[1] = log( x[0] ) * x[1];
		F[2] = sin( x[2] );
	}
	
	if( needG[0] > 0 ) {
		Gdata[0] = cos( x[2] );					// (3,3)
		Gdata[1] = 0.25 * exp( x[1] / 4.0 );	// (1,2)
		Gdata[2] = x[1] / x[0];					// (2,1)
		Gdata[3] = 3.0 * x[0] * x[0];			// (1,1)
		Gdata[4] = log( x[0] );					// (2,2)
	}
	
	return 0;
}

void mixtest()
{
	integer N, M, FN;
	double *x, *xLoBnds, *xUpBnds;
	integer Gnnz, *Grows, *Gcols;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * PROBLEM DATA  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	N  = 3;
	M  = 2;
	FN = 3;
	
	Gnnz = 5; // 5 out of 9 nonzeros in Jacobian
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * ALLOCATE PROBLEM VARIABLES  * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// problem variables
    x			 = (double *) calloc (N, sizeof(double));
    xLoBnds      = (double *) calloc (N, sizeof(double));
    xUpBnds      = (double *) calloc (N, sizeof(double));
	
	// Jacobian of the nonlinear part of objective and constraints
    Grows		 = (integer *)calloc (Gnnz, sizeof(integer));
    Gcols		 = (integer *)calloc (Gnnz, sizeof(integer));
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * INITIALIZE SNOPT PROBLEM DATA * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// x point (x[0] cannot be zero for differentiability)
	x[0] = 1.0;
	x[1] = 2.0;
	x[2] = 3.0;
	
	// x bound(s)
	xLoBnds[0] = -1.0e20; xUpBnds[0] =  1.0e20;
	xLoBnds[1] = -1.0e20; xUpBnds[1] =  1.0e20;
	xLoBnds[2] = -1.0e20; xUpBnds[2] =  1.0e20;
	
	// Constraint Jacobian sparsity pattern using FORTRAN-style indices, 
	// and a non-standard pattern (i.e., not row- or column-major)
	// SNOPT shouldn't care about the order of the appearance of derivatives
	Grows[0] = 3; Gcols[0] = 3;
	Grows[1] = 1; Gcols[1] = 2;
	Grows[2] = 2; Gcols[2] = 1;
	Grows[3] = 1; Gcols[3] = 1;
	Grows[4] = 2; Gcols[4] = 2;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// check jacobian 
	snopta_eval_G_check(N, FN,
						mixtest_fun_,
						Gnnz, Grows, Gcols, 
						x, xLoBnds, xUpBnds,
						NULL, NULL, 
						NULL, NULL,
						NULL, NULL);
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	free(x);
	free(xLoBnds);
	free(xUpBnds);
	free(Grows);
	free(Gcols);
	
}