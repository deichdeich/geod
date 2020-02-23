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
#include "definitions.h"


int check_err(double tol){
    //printf("check_err: ");
    //print_vec(integration_vectors.err_vec);
    //printf("%f\n", tol);
    int good_err = 1;
    for(int i = 0; i < 6; i++){
        if(gsl_vector_get(integration_vectors.err_vec, i) > tol){
            //printf("%f\n", gsl_vector_get(integration_vectors.err_vec, i));
            good_err = 0;
        }
    }
    return good_err;
}

int stepper2(int (*f) (double, gsl_vector *, gsl_vector *),
	     int dof,
	     gsl_vector * in_state_vec,
	     gsl_vector * out_state_vec,
	     double x,
	     double h,
	     double xmax,
	     double tolerance){
    double err;
    double scale;
    int step = 0;

    cur_state = gsl_vector_calloc(dof);
    gsl_vector_memcpy(integration_vectors.temp_in_vec, in_state_vec);
    gsl_vector_memcpy(integration_vectors.temp_out_vec, out_state_vec);
    gsl_vector_memcpy(cur_state, in_state_vec);
    while (x < xmax+h){
        int attempts = 0;
        int err_val = 0;

        err = fabs(int_func(f,
                        dof,
                        cur_state,
                        integration_vectors.temp_out_vec,
                        x, h));
        err_val = check_err(tolerance);
        
        while (err_val == 0){
            scale = pow(0.8 * h * (1/tolerance),(1./7.));
            scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
            h *= scale;
            err = fabs(int_func(f,
                        dof,
                        cur_state,
                        integration_vectors.temp_out_vec,
                        x, h));
            err_val = check_err(tolerance);
            attempts += 1;
        
            if(attempts > 100){
                printf(">100 attempts to get the error below tolerance");
                return -2;
                }
            
        }

        if (make_section){
            poincare_check(cur_state);
            if (poincare_yes){
                add_to_history = 1;
                }
            }


        if(add_to_history){
            populate_history(dof, lines_written, x, cur_state);
            lines_written++;
            if(make_section){
                add_to_history = 0;
                poincare_yes = 0;
            }
        }
        
        // clean up integration_vector struct, copy out state to file
        gsl_vector_memcpy(cur_state,
                          integration_vectors.temp_out_vec);
        integration_vector_zero(dof);
        
        step += 1;
        
        // advance timestep
        x += h;
        
    }	
    write_history(lines_written%hist_len);
    integration_vector_free();
    //free(history);

    return 0;     

}    
 
/////////////////////
// stepper steps the state from one timestamp to another.  It makes sure that the error
// in the final state is less than a given tolerance.
//
// Arguments:
//      *f: pointer to the function which gives the derivative of a system at a given time
//          and given initial condition.  It's an equation of motion.
//      in_state_vec: pointer to the initial condition, gsl_vector
//      out_state_vec: pointer to the state, to be filled, at the time x0+h, gsl_vector
//      x: the clock value corresponding to in_vec, double
//      h: the initial timestep guess, double
//      xmax: the final clock value to integrate to, double
//      *h_next: pointer to where the next timestep is stored, double
//      tolerance:  The error in the integration is guaranteed to be less than this, double
//
//  Returns:
//      0 if successful
//      No failure case yet
//
/////////////////////


int stepper(int (*f) (double, gsl_vector *, gsl_vector *),
	     int dof,
	     gsl_vector * in_state_vec,
	     gsl_vector * out_state_vec,
	     double x,
	     double h,
	     double xmax,
	     double *h_next,
	     double tolerance)
{


    integration_vector_init(dof);
    const double err_exponent = 1.0 / 7.0;

    double scale;
    double err;
    double yy;
    int i;
    int last_interval = 0;

    // Verify that the step size is positive and that the upper endpoint //
    // of integration is greater than the initial enpoint.               //

    if (xmax < x || h <= 0.0)
	return -2;

    // If the upper endpoint of the independent variable agrees with the //
    // initial value of the independent variable.  Set the value of the  //
    // dependent variable and return success.                            //

    *h_next = h;
    gsl_vector_memcpy(out_state_vec, in_state_vec);
    if (xmax == x)
	return 0;

    // ensure that the step size h is not larger than the length of the //
    // integration interval.                                            //

    if (h > (xmax - x)) {
	h = xmax - x;
	last_interval = 1;
    }
    // Redefine the error tolerance to an error tolerance per unit    //
    // length of the integration interval.                            //

    tolerance /= (xmax - x);

    // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
    // maintain an error less than tolerance * (xmax-x) using an     //
    // initial step size of h and initial value: y = y[0]            //

    gsl_vector_memcpy(integration_vectors.temp_in_vec, in_state_vec);
    int foo = 0;
    while (x < xmax) {
        //printf("%f, %f\n",x, xmax);
	    scale = 1.0;
	    i = 0;
        for (i = 0; i < ATTEMPTS; i++) {
            err = fabs(single_stepRKF78(f,
                        dof,
                        integration_vectors.temp_in_vec,
                        integration_vectors.temp_out_vec,
                        x, h));
            printf("%0.10f, %0.10f, %0.10f, %0.10f, %0.10f\n\n", err, tolerance, h, x, xmax);
        
            if (err == 0.0) {
                scale = MAX_SCALE_FACTOR;
                break;
            }

            double statenorm = gsl_blas_dnrm2(integration_vectors.temp_in_vec);
            yy = (statenorm == 0.0) ? tolerance : statenorm;
            scale = 0.8 * pow(tolerance * yy / err, err_exponent);
            scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);

            if (err < (tolerance * yy)){
                break;
                }
            h *= scale;

            if (x + h > xmax){
                h = xmax - x;
                }
            
            else if (x + h + 0.5 * h > xmax){
                h = 0.5 * h;
                }
            }
        
        if (make_section){
            poincare_check(integration_vectors.temp_out_vec);
            if (poincare_yes){
                add_to_history = 1;
                }
            }

        //add_to_file(dof, x, integration_vectors.temp_out_vec);

        if (add_to_history){
	        
	        if(make_section){
	            add_to_history = 0;
	            poincare_yes = 0;
	        }
	    }

        if (i >= ATTEMPTS) {
            *h_next = h * scale;
            return -1;
        }
        
        gsl_vector_memcpy(integration_vectors.temp_in_vec,
                  integration_vectors.temp_out_vec);
        x += h;
        printf("%f\n", x);
        h *= scale;
        *h_next = h;

        if (last_interval) {
            printf("last interval\n");
            break;
        }

        if (x + h > xmax) {
            last_interval = 1;
            h = xmax - x;
        }

        else if (x + h + 0.5 * h > xmax)
            h = 0.5 * h;
    foo++;
    }
    //printf("blah");

    gsl_vector_memcpy(out_state_vec, integration_vectors.temp_out_vec);
    integration_vector_free();
    return 0;
}
