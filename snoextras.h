/*
 *  snoextras.h
 *  
 *
 *  Created by W. Ross Morrow on 3/19/12.
 *  Copyright 2012 Iowa State University. All rights reserved.
 *
 */

#include "f2c.h"
#include "snopt.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Check first derivatives for a SNOPT(A) callback routine
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
					   double		*x,			// point at which we want to check first derivatives
					   double		*xLoBnds,	// lower bounds on variables
					   double		*xUpBnds,	// upper bounds on variables
					   char			*cu,		// character workspace
					   integer		*lencu,		// length of character workspace
					   integer		*iu,		// integer workspace
					   integer		*leniu,		// length of integer workspace
					   doublereal	*ru,		// double workspace
					   integer		*lenru );	// length of double workspace

/*
 
int usrfun(
 integer		*Status,	// SNOPT status code
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
 
*/