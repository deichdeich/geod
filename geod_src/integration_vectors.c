#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "ad_rkf78.h"
#include "eoms.h"
#include "integration_vectors.h"
#include "single_step.h"
#include "steppers.h"
#include "fileio.h"

/* declaring the integration vector struct*/


///////////////////
// integration_vector_init initializes the integration_vectors struct
// integration_vector_free frees the struct
// No arguments
////////////////////
int init_count = 0;
void integration_vector_init(int dof)
{
    init_count += 1;
    int i, Nelt;

    gsl_vector **p;		// pointer to a vector

    Nelt = sizeof(integration_vectors) / sizeof(integration_vectors.k1_in);	// number of elements in x.

    p = &integration_vectors.k1_in;	// pointing to the first vector in the struct

    for (i = 0; i < Nelt; i++, p++) {	// p[0] points to integration_vectors.k1_in,
	    *p = gsl_vector_calloc(dof);	// p[1] points to integration_vectors.k2_in, etc...
    }
}

void integration_vector_zero(int dof)
{
    int i, Nelt, j;
    gsl_vector **p;

    Nelt = sizeof(integration_vectors) / sizeof(integration_vectors.k1_in);	// number of elements in x.

    p = &integration_vectors.k1_in;
    for (i = 0; i < Nelt; i++, p++) {
        for (j = 0; j < dof; j++){
	        gsl_vector_set(*p, j, 0);
	    }
    }
}

void integration_vector_free()
{
    int i, Nelt;
    gsl_vector **p;

    Nelt = sizeof(integration_vectors) / sizeof(integration_vectors.k1_in);	// number of elements in x.

    p = &integration_vectors.k1_in;
    for (i = 0; i < Nelt; i++, p++) {
	    gsl_vector_free(*p);
    }
}
