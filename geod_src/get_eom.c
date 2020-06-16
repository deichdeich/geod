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

void get_eom(char eom_name[]){
    if(strcmp("kerr_eom_exact", eom_name) == 0){
        f = kerr_eom_exact;
    }
    else if(strcmp("kerr_eom_SR_o2", eom_name) == 0){
        f = kerr_eom_SR_o2;
    }
    else if(strcmp("edgb_a2_z1", eom_name) == 0){
        f = edgb_a2_z1;
    }
}
