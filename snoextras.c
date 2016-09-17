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

#include <vecLib/cblas.h>

#include "snoextras.h"

#define ABS(a)     ( (a < 0) ? -(a) : (a) )
#define MAX(a,b)   ( (a < b) ?  (b) : (a) )
#define MIN(a,b)   ( (a > b) ?  (b) : (a) )

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Check first derivatives for a SNOPT callback routine
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void snopta_eval_G_check(integer		 N,			// number of variables
						 integer		 M,			// number of nonlinear functions
						 U_fp			 usrfun,	// callback to evaluate nonlinear function values and derivatives
						 integer		 Gnnz,		// (constant) constraint Jacobian number of nonzeros
						 integer		*Grows,		// (constant) constraint Jacobian structure - row indices
						 integer		*Gcols,		// (constant) constraint Jacobian structure - column indices
						 double			*x,			// point at which we want to check first derivatives
						 double			*xLoBnds,	// lower bounds on variables
						 double			*xUpBnds,	// upper bounds on variables
						 char			*cu,		// character workspace
						 integer		*lencu,		// length of character workspace
						 integer		*iu,		// integer workspace
						 integer		*leniu,		// length of integer workspace
						 doublereal		*ru,		// double workspace
						 integer		*lenru )	// length of double workspace
{
	int S = 16;
	
	integer INFO;
	integer needF = 1;
	integer needG = 1;
	
	double * curr_F;
	double * curr_G;
	double * new_x;
	double * new_F;
	double * fdG;
	
	double h, hinv, H, Hinv;
	
	double abstmp, reltmp;
	double absGdiffnorm[S];
	double relGdiffnorm[S];
	
	int s, m, n, e;
	
	int flag;
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// need several more argument checks to be robust. 
	if( usrfun == NULL || x == NULL ) { return; }
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// allocate memory needed
	curr_F  = (double *)calloc( M    , sizeof(double) );
	curr_G  = (double *)calloc( Gnnz , sizeof(double) );
	new_x   = (double *)calloc( N    , sizeof(double) );
	new_F   = (double *)calloc( M    , sizeof(double) );
	fdG     = (double *)calloc( M*N  , sizeof(double) );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	printf("\n\n");
	printf("SNO(A)_EVALGA_CHECK:: F : R(%i) -> R(%i)\n",(int)N,(int)M);
	printf("SNO(A)_EVALGA_CHECK:: DF is %i x %i, with %i nonzeros (%0.4f %% dense)\n",(int)M,(int)N,(int)Gnnz,100.0*(double)Gnnz/((double)(M*N)));
	// for( e = 0 ; e < Gnnz ; e++ ) {
	//	printf("                      element %i gives (%i,%i)\n",e,Grows[e],Gcols[e]);
	//}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// ensure x is within bounds (using Euclidean projection)
	for( n = 0 ; n < N ; n++ ) {
		if( x[n] < xLoBnds[n] ) { x[n] = xLoBnds[n]; }
		if( x[n] > xUpBnds[n] ) { x[n] = xUpBnds[n]; }
	}
	
	printf("SNO(A)_EVALGA_CHECK:: Current Point (projected to bounds) is\n");
	printf("  x = [ %0.6f ",x[0]);
	for( n = 1 ; n < N ; n++ ) { printf(", %0.6f ",x[n]); }
	printf("]\n");
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// get both function and Jacobian
	needF = 1;
	needG = 1;
	usrfun(&INFO,  &N,	  x, 
		   &needF, &M,	  curr_F, 
		   &needG, &Gnnz, curr_G, 
		   cu, lencu, 
		   iu, leniu, 
		   ru, lenru);
	
	// don't compute Jacobian any more
	needG = 0;
	
	printf("SNO(A)_EVALGA_CHECK:: Jacobian Differences are...\n");
	
	// finite differences
	for( s = 0 ; s < S ; s++ ) {
		
		// stepsize
		h = pow( 10.0 , -1.0 * s ); hinv = 1.0 / h;
		
		printf("  h = %0.16f, ",h);
		
		// for each direction...
		for( n = 0 ; n < N ; n++ ) {
			
			// new_x <- x
			cblas_dcopy( N , x , 1 , new_x , 1 );
			
			// relative step size
			H = ( 1.0 + fabs(x[n]) ) * h; Hinv = 1.0 / H;
			
			// respect bounds, with new point
			if( x[n] + H <= xUpBnds[n] ) {
				
				// use a * forward * difference, safeguarded to lie within bounds
				
				// new_x <- min{ x + H * e_n , xUpBnds[n] }
				new_x[n] += H;
				
				// get *new* nonlinear function values
				usrfun(&INFO,  &N, new_x, 
					   &needF, &M, new_F, 
					   &needG, NULL, NULL, // don't needG, 
					   cu, lencu, 
					   iu, leniu, 
					   ru, lenru);
				
				// form * forward * finite differences
				
				// fdG[:,n] = ( new_F - curr_F ) / h
				cblas_dcopy( M , new_F , 1 , fdG+M*n , 1 );
				cblas_daxpy( M , -1.0 , curr_F , 1 , fdG+M*n , 1 );
				cblas_dscal( M , Hinv , fdG+M*n , 1 );
				
			} else {
				
				// we assume (or assert) that 
				//
				//		xLoBnds[n] <= x[n] <= xUpBnds[n]
				//
				// thus, if the if statement above has been violated, 
				// we must have that x[n] + H > xUpBnds[n]
				
				// here use a * backward * difference, safeguarded to lie within bounds
				
				// new_x <- max{ x - H * e_n , xLoBnds[n] }
				new_x[n] -= H;
				if( new_x[n] < xLoBnds[n] ) { new_x[n] = xLoBnds[n]; }
				
				// get *new* nonlinear function values
				usrfun(&INFO,  &N, new_x, 
					   &needF, &M, new_F, 
					   &needG, NULL, NULL, 
					   cu, lencu, 
					   iu, leniu, 
					   ru, lenru);
				
				// form * backward * finite differences
				
				// fdG[:,n] = ( curr_F - new_F ) / h
				cblas_dcopy( M , curr_F , 1 , fdG+M*n , 1 );
				cblas_daxpy( M , -1.0 , new_F , 1 , fdG+M*n , 1 );
				cblas_dscal( M , Hinv , fdG+M*n , 1 );
				
			}
			
		}
		
		// finite difference Jacobian has been formed
		// compare to computed Jacobian
		
		// first assess sparsity pattern
		for( m = 0 ; m < M ; m++ ) {
			for( n = 0 ; n < N ; n++ ) {
				
				if( fdG[ m + M * n ] != 0.0 ) {
					
					// make sure this entry is included in the sparse structure
					// (using FORTRAN-style indexing)
					flag = 0;
					for( e = 0 ; e < Gnnz ; e++ ) {
						if( Grows[e] - 1 == m && Gcols[e] - 1 == n ) { flag = 1; break; }
					}
					if( flag == 0 ) {
						printf("SNO(A)_EVALGA_CHECK WARNING:: element (%i,%i) in the Jacobian has a nonzero finite difference, \n",m+1,n+1);
						printf("                              but is not included in the given sparsity structure.\n");
					}
					
				}
			}
		}
		
		// then compute differences (we can overwrite elements of finite difference)
		absGdiffnorm[s] = 0.0;
		relGdiffnorm[s] = 0.0;
		for( e = 0 ; e < Gnnz ; e++ ) {
			
			m = Grows[e]-1;
			n = Gcols[e]-1;
			
			abstmp = ABS( fdG[ m + M * n ] - curr_G[e] );
			absGdiffnorm[s] = MAX( abstmp , absGdiffnorm[s] );
			
			// take denominator to be max{ fdG[.] , 1 }
			reltmp = abstmp / MAX( fabs( fdG[ m + M * n ] ) , 1.0 );
			relGdiffnorm[s] = MAX( relGdiffnorm[s] , reltmp );
			
		}
		
		printf("|| abs diff ||_inf = %0.16f, || rel diff ||_inf = %0.16f\n",absGdiffnorm[s],relGdiffnorm[s]);
		
	}
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
	// free memory allocated above
	free( curr_F );
	free( curr_G );
	free(  new_x );
	free(  new_F );
	free(    fdG );
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
}
