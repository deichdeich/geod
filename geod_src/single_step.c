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


double single_stepRK4(int (*f) (double, gsl_vector *, gsl_vector *),
                          int dof,
                          gsl_vector * in_vec,
                          gsl_vector * out_vec,
                          double x0,
                          double h){
    //print_vec(out_vec);
    
    gsl_vector_memcpy(integration_vectors.k1_in, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k2_in, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k3_in, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k4_in, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k5_in, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k1_out, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k2_out, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k3_out, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k4_out, integration_vectors.zero_vec);
    gsl_vector_memcpy(integration_vectors.k5_out, integration_vectors.zero_vec);
    
    (*f) (x0, in_vec, integration_vectors.k1_in);
    gsl_blas_dscal(h, integration_vectors.k1_in);
    /*
    k1 = dt * derivative_function(state)
    k2 = dt * derivative_function(state + (k1 / 2))
    k3 = dt * derivative_function(state + (k2 / 2))
    k4 = dt * derivative_function(state + k3)
    k = (1 / 6) * (k1 + (2 * k2) + (2 * k3) + k4)
    */
    // divide k1 by 2
    gsl_blas_dscal(1./2., integration_vectors.k1_in);
    
    // add it to the state
    vect_add(integration_vectors.k1_out, 6, 2, in_vec, integration_vectors.k1_in);
    // get k2
    (*f) (x0, integration_vectors.k1_out, integration_vectors.k2_in);
    // multiply k2 by dt
    gsl_blas_dscal(h, integration_vectors.k2_in);

    
    // divide k2 by 2
    gsl_blas_dscal(1./2., integration_vectors.k2_in);
    // add it to the state
    vect_add(integration_vectors.k2_out, 6, 2, in_vec, integration_vectors.k2_in);
    // get k3
    (*f) (x0, integration_vectors.k2_out, integration_vectors.k3_in);
    // multiply k3 by dt
    gsl_blas_dscal(h, integration_vectors.k3_in);

    

    // add it to the state
    vect_add(integration_vectors.k3_out, 6, 2, in_vec, integration_vectors.k3_in);
    (*f) (x0, integration_vectors.k3_out, integration_vectors.k4_in);
    gsl_blas_dscal(h, integration_vectors.k4_in);

    // scale the k's
    gsl_blas_dscal(2., integration_vectors.k1_in);
    gsl_blas_dscal(4., integration_vectors.k2_in);
    gsl_blas_dscal(2., integration_vectors.k3_in);
    
    vect_add(integration_vectors.k5_in, 6, 4, integration_vectors.k1_in,
                                              integration_vectors.k2_in,
                                              integration_vectors.k3_in,
                                              integration_vectors.k4_in);
    
    gsl_blas_dscal(1./6.,integration_vectors.k5_in);
    
    vect_add(out_vec, 6, 2, in_vec, integration_vectors.k5_in);

    
    return 0.;
    }


////////////////
// single_step computes the next value of the system for a single timestep.
// It's a shitty, awful function that just encodes the RKF7(8) Butcher tableau.
// I hate editing this function.  Don't touch it if you don't think you have to.
// It basically makes sense if you look at any pseudocode implementation of RKF7(8); not sure which source I used. 
//
// It relies on the integration_vectors struct to hold each interstitial integration vector.
// 
// Arguments:
//      *f: pointer to the function which gives the derivative of a system at a given time
//          and given initial condition.  It's an equation of motion.
//      in_vec: pointer to the initial condition, gsl_vector
//      out_vec: pointer to the state, to be filled, at the time x0+h, gsl_vector
//      x0: the clock value corresponding to in_vec
//      h: the timestep
//
// Returns:
//      err: The error in the 7th-order timestep as evaluated by the 8th-order. This is used
//           by stepper() to _determine the size of the next timestep (if err is too big).  Double.
//
////////////////
 double single_stepRKF78(int (*f) (double, gsl_vector *, gsl_vector *),
                          int dof,
                          gsl_vector * in_vec,
                          gsl_vector * out_vec,
                          double x0,
                          double h,
                          int check){
     if(check ==1){
        printf("input values:\n\t in_vec: ");
        print_vec(in_vec);
        printf("\t out_vec: ");
        print_vec(out_vec);
     }

    const double c_1_11 = 41.0 / 840.0;
    const double c6 = 34.0 / 105.0;
    const double c_7_8 = 9.0 / 35.0;
    const double c_9_10 = 9.0 / 280.0;

    const double a2 = 2.0 / 27.0;
    const double a3 = 1.0 / 9.0;
    const double a4 = 1.0 / 6.0;
    const double a5 = 5.0 / 12.0;
    const double a6 = 1.0 / 2.0;
    const double a7 = 5.0 / 6.0;
    const double a8 = 1.0 / 6.0;
    const double a9 = 2.0 / 3.0;
    const double a10 = 1.0 / 3.0;

    const double b31 = 1.0 / 36.0;
    const double b32 = 3.0 / 36.0;
    const double b41 = 1.0 / 24.0;
    const double b43 = 3.0 / 24.0;
    const double b51 = 20.0 / 48.0;
    const double b53 = -75.0 / 48.0;
    const double b54 = 75.0 / 48.0;
    const double b61 = 1.0 / 20.0;
    const double b64 = 5.0 / 20.0;
    const double b65 = 4.0 / 20.0;
    const double b71 = -25.0 / 108.0;
    const double b74 = 125.0 / 108.0;
    const double b75 = -260.0 / 108.0;
    const double b76 = 250.0 / 108.0;
    const double b81 = 31.0 / 300.0;
    const double b85 = 61.0 / 225.0;
    const double b86 = -2.0 / 9.0;
    const double b87 = 13.0 / 900.0;
    const double b91 = 2.0;
    const double b94 = -53.0 / 6.0;
    const double b95 = 704.0 / 45.0;
    const double b96 = -107.0 / 9.0;
    const double b97 = 67.0 / 90.0;
    const double b98 = 3.0;
    const double b10_1 = -91.0 / 108.0;
    const double b10_4 = 23.0 / 108.0;
    const double b10_5 = -976.0 / 135.0;
    const double b10_6 = 311.0 / 54.0;
    const double b10_7 = -19.0 / 60.0;
    const double b10_8 = 17.0 / 6.0;
    const double b10_9 = -1.0 / 12.0;
    const double b11_1 = 2383.0 / 4100.0;
    const double b11_4 = -341.0 / 164.0;
    const double b11_5 = 4496.0 / 1025.0;
    const double b11_6 = -301.0 / 82.0;
    const double b11_7 = 2133.0 / 4100.0;
    const double b11_8 = 45.0 / 82.0;
    const double b11_9 = 45.0 / 164.0;
    const double b11_10 = 18.0 / 41.0;
    const double b12_1 = 3.0 / 205.0;
    const double b12_6 = -6.0 / 41.0;
    const double b12_7 = -3.0 / 205.0;
    const double b12_8 = -3.0 / 41.0;
    const double b12_9 = 3.0 / 41.0;
    const double b12_10 = 6.0 / 41.0;
    const double b13_1 = -1777.0 / 4100.0;
    const double b13_4 = -341.0 / 164.0;
    const double b13_5 = 4496.0 / 1025.0;
    const double b13_6 = -289.0 / 82.0;
    const double b13_7 = 2193.0 / 4100.0;
    const double b13_8 = 51.0 / 82.0;
    const double b13_9 = 33.0 / 164.0;
    const double b13_10 = 12.0 / 41.0;
    const double b13_12 = 1.0;

    const double err_factor = -41.0 / 840.0;
    double h2_7 = a2 * h;

    //////////////
    // k1
    /////////////
    if(check ==1){
        printf("\t integration_vectors.k1_out: ");
        print_vec(integration_vectors.k1_out);
 }
    gsl_vector_memcpy(integration_vectors.k1_in, in_vec);
    (*f) (x0, integration_vectors.k1_in, integration_vectors.k1_out);
    if(check==1) printf("\t integration_vectors.k1_out: "); //print_vec(integration_vectors.k1_out);

    //////////////
    // k2
    ///////////// 
    gsl_vector_memcpy(integration_vectors.k2_in,
		      integration_vectors.k1_out);
    gsl_blas_dscal(h * a2, integration_vectors.k2_in);
    gsl_vector_add(integration_vectors.k2_in, in_vec);

    (*f) (x0 + h2_7, integration_vectors.k2_in,
	  integration_vectors.k2_out);

    //////////////
    // k3
    /////////////
    gsl_vector_memcpy(integration_vectors.k3k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k3k2_out,
		      integration_vectors.k2_out);
    gsl_blas_dscal(b31, integration_vectors.k3k1_out);
    gsl_blas_dscal(b32, integration_vectors.k3k2_out);
    vect_add(integration_vectors.k3_in, dof, 2,
	      integration_vectors.k3k1_out, integration_vectors.k3k2_out);
    gsl_blas_dscal(h, integration_vectors.k3_in);
    gsl_vector_add(integration_vectors.k3_in, in_vec);

    (*f) (x0 + a3 * h, integration_vectors.k3_in,
	  integration_vectors.k3_out);
	  
    //////////////
    // k4
    /////////////
    gsl_vector_memcpy(integration_vectors.k4k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k4k3_out,
		      integration_vectors.k3_out);
    gsl_blas_dscal(b41, integration_vectors.k4k1_out);
    gsl_blas_dscal(b43, integration_vectors.k4k3_out);
    vect_add(integration_vectors.k4_in, dof, 2,
	      integration_vectors.k4k1_out, integration_vectors.k4k3_out);
    gsl_blas_dscal(h, integration_vectors.k4_in);
    gsl_vector_add(integration_vectors.k4_in, in_vec);

    (*f) (x0 + a4 * h, integration_vectors.k4_in,
	  integration_vectors.k4_out);
	  
    //////////////
    // k5
    /////////////
    gsl_vector_memcpy(integration_vectors.k5k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k5k3_out,
		      integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k5k4_out,
		      integration_vectors.k4_out);
    gsl_blas_dscal(b51, integration_vectors.k5k1_out);
    gsl_blas_dscal(b53, integration_vectors.k5k3_out);
    gsl_blas_dscal(b54, integration_vectors.k5k4_out);
    vect_add(integration_vectors.k5_in, dof, 3,
	      integration_vectors.k5k1_out, integration_vectors.k5k3_out,
	      integration_vectors.k5k4_out);
    gsl_blas_dscal(h, integration_vectors.k5_in);
    gsl_vector_add(integration_vectors.k5_in, in_vec);

    (*f) (x0 + a5 * h, integration_vectors.k5_in,
	  integration_vectors.k5_out);

    /////////////
    // k6
    /////////////
    gsl_vector_memcpy(integration_vectors.k6k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k6k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k6k5_out,
		      integration_vectors.k5_out);
    gsl_blas_dscal(b61, integration_vectors.k6k1_out);
    gsl_blas_dscal(b64, integration_vectors.k6k4_out);
    gsl_blas_dscal(b65, integration_vectors.k6k5_out);
    vect_add(integration_vectors.k6_in, dof, 3,
	      integration_vectors.k6k1_out, integration_vectors.k6k4_out,
	      integration_vectors.k6k5_out);
    gsl_blas_dscal(h, integration_vectors.k6_in);
    gsl_vector_add(integration_vectors.k6_in, in_vec);

    (*f) (x0 + a6 * h, integration_vectors.k6_in,
	  integration_vectors.k6_out);

    //////////////
    // k7
    /////////////
    gsl_vector_memcpy(integration_vectors.k7k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k7k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k7k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k7k6_out,
		      integration_vectors.k6_out);
    gsl_blas_dscal(b71, integration_vectors.k7k1_out);
    gsl_blas_dscal(b74, integration_vectors.k7k4_out);
    gsl_blas_dscal(b75, integration_vectors.k7k5_out);
    gsl_blas_dscal(b76, integration_vectors.k7k6_out);
    vect_add(integration_vectors.k7_in, dof, 4,
	      integration_vectors.k7k1_out, integration_vectors.k7k4_out,
	      integration_vectors.k7k5_out, integration_vectors.k7k6_out);
    gsl_blas_dscal(h, integration_vectors.k7_in);
    gsl_vector_add(integration_vectors.k7_in, in_vec);

    (*f) (x0 + a7 * h, integration_vectors.k7_in,
	  integration_vectors.k7_out);

    //////////////
    // k8
    /////////////
    gsl_vector_memcpy(integration_vectors.k8k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k8k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k8k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k8k7_out,
		      integration_vectors.k7_out);
    gsl_blas_dscal(b81, integration_vectors.k8k1_out);
    gsl_blas_dscal(b85, integration_vectors.k8k5_out);
    gsl_blas_dscal(b86, integration_vectors.k8k6_out);
    gsl_blas_dscal(b87, integration_vectors.k8k7_out);
    vect_add(integration_vectors.k8_in, dof, 4,
	      integration_vectors.k8k1_out, integration_vectors.k8k5_out,
	      integration_vectors.k8k6_out, integration_vectors.k8k7_out);
    gsl_blas_dscal(h, integration_vectors.k8_in);
    gsl_vector_add(integration_vectors.k8_in, in_vec);

    (*f) (x0 + a8 * h, integration_vectors.k8_in,
	  integration_vectors.k8_out);

    //////////////
    // k9
    /////////////
    gsl_vector_memcpy(integration_vectors.k9k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k9k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k9k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k9k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k9k7_out,
		      integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k9k8_out,
		      integration_vectors.k8_out);
    gsl_blas_dscal(b91, integration_vectors.k9k1_out);
    gsl_blas_dscal(b94, integration_vectors.k9k4_out);
    gsl_blas_dscal(b95, integration_vectors.k9k5_out);
    gsl_blas_dscal(b96, integration_vectors.k9k6_out);
    gsl_blas_dscal(b97, integration_vectors.k9k7_out);
    gsl_blas_dscal(b98, integration_vectors.k9k8_out);
    vect_add(integration_vectors.k9_in, dof, 6,
	      integration_vectors.k9k1_out, integration_vectors.k9k4_out,
	      integration_vectors.k9k5_out, integration_vectors.k9k6_out,
	      integration_vectors.k9k7_out, integration_vectors.k9k8_out);

    gsl_blas_dscal(h, integration_vectors.k9_in);
    gsl_vector_add(integration_vectors.k9_in, in_vec);
    
    (*f) (x0 + a9 * h, integration_vectors.k9_in,
	  integration_vectors.k9_out);

    //////////////
    // k10
    /////////////
    gsl_vector_memcpy(integration_vectors.k10k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k10k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k10k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k10k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k10k7_out,
		      integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k10k8_out,
		      integration_vectors.k8_out);
    gsl_vector_memcpy(integration_vectors.k10k9_out,
		      integration_vectors.k9_out);
    gsl_blas_dscal(b10_1, integration_vectors.k10k1_out);
    gsl_blas_dscal(b10_4, integration_vectors.k10k4_out);
    gsl_blas_dscal(b10_5, integration_vectors.k10k5_out);
    gsl_blas_dscal(b10_6, integration_vectors.k10k6_out);
    gsl_blas_dscal(b10_7, integration_vectors.k10k7_out);
    gsl_blas_dscal(b10_8, integration_vectors.k10k8_out);
    gsl_blas_dscal(b10_9, integration_vectors.k10k9_out);
    vect_add(integration_vectors.k10_in, dof, 7,
	      integration_vectors.k10k1_out, integration_vectors.k10k4_out,
	      integration_vectors.k10k5_out, integration_vectors.k10k6_out,
	      integration_vectors.k10k7_out, integration_vectors.k10k8_out,
	      integration_vectors.k10k9_out);
    gsl_blas_dscal(h, integration_vectors.k10_in);
    gsl_vector_add(integration_vectors.k10_in, in_vec);
    
    (*f) (x0 + a10 * h, integration_vectors.k10_in,
	  integration_vectors.k10_out);

    //////////////
    // k11
    /////////////
    gsl_vector_memcpy(integration_vectors.k11k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k11k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k11k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k11k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k11k7_out,
		      integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k11k8_out,
		      integration_vectors.k8_out);
    gsl_vector_memcpy(integration_vectors.k11k9_out,
		      integration_vectors.k9_out);
    gsl_vector_memcpy(integration_vectors.k11k10_out,
		      integration_vectors.k10_out);
    gsl_blas_dscal(b11_1, integration_vectors.k11k1_out);
    gsl_blas_dscal(b11_4, integration_vectors.k11k4_out);
    gsl_blas_dscal(b11_5, integration_vectors.k11k5_out);
    gsl_blas_dscal(b11_6, integration_vectors.k11k6_out);
    gsl_blas_dscal(b11_7, integration_vectors.k11k7_out);
    gsl_blas_dscal(b11_8, integration_vectors.k11k8_out);
    gsl_blas_dscal(b11_9, integration_vectors.k11k9_out);
    gsl_blas_dscal(b11_10, integration_vectors.k11k10_out);
    vect_add(integration_vectors.k11_in, dof, 8,
	      integration_vectors.k11k1_out, integration_vectors.k11k4_out,
	      integration_vectors.k11k5_out, integration_vectors.k11k6_out,
	      integration_vectors.k11k7_out, integration_vectors.k11k8_out,
	      integration_vectors.k11k9_out,
	      integration_vectors.k11k10_out);
    gsl_blas_dscal(h, integration_vectors.k11_in);
    gsl_vector_add(integration_vectors.k11_in, in_vec);
    
    (*f) (x0 + h, integration_vectors.k11_in, integration_vectors.k11_out);

    //////////////
    // k12
    /////////////
    gsl_vector_memcpy(integration_vectors.k12k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k12k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k12k7_out,
		      integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k12k8_out,
		      integration_vectors.k8_out);
    gsl_vector_memcpy(integration_vectors.k12k9_out,
		      integration_vectors.k9_out);
    gsl_vector_memcpy(integration_vectors.k12k10_out,
		      integration_vectors.k10_out);
    gsl_blas_dscal(b12_1, integration_vectors.k12k1_out);
    gsl_blas_dscal(b12_6, integration_vectors.k12k6_out);
    gsl_blas_dscal(b12_7, integration_vectors.k12k7_out);
    gsl_blas_dscal(b12_8, integration_vectors.k12k8_out);
    gsl_blas_dscal(b12_9, integration_vectors.k12k9_out);
    gsl_blas_dscal(b12_10, integration_vectors.k12k10_out);
    vect_add(integration_vectors.k12_in, dof, 6,
	      integration_vectors.k12k1_out, integration_vectors.k12k6_out,
	      integration_vectors.k12k7_out, integration_vectors.k12k8_out,
	      integration_vectors.k12k9_out,
	      integration_vectors.k12k10_out);
    gsl_blas_dscal(h, integration_vectors.k12_in);
    gsl_vector_add(integration_vectors.k12_in, in_vec);

    (*f) (x0, integration_vectors.k12_in, integration_vectors.k12_out);

    //////////////
    // k13
    /////////////
    gsl_vector_memcpy(integration_vectors.k13k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k13k4_out,
		      integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k13k5_out,
		      integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k13k6_out,
		      integration_vectors.k6_out);
    gsl_vector_memcpy(integration_vectors.k13k7_out,
		      integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k13k8_out,
		      integration_vectors.k8_out);
    gsl_vector_memcpy(integration_vectors.k13k9_out,
		      integration_vectors.k9_out);
    gsl_vector_memcpy(integration_vectors.k13k10_out,
		      integration_vectors.k10_out);
    gsl_vector_memcpy(integration_vectors.k13k12_out,
		      integration_vectors.k12_out);
    gsl_blas_dscal(b13_1, integration_vectors.k13k1_out);
    gsl_blas_dscal(b13_4, integration_vectors.k13k4_out);
    gsl_blas_dscal(b13_5, integration_vectors.k13k5_out);
    gsl_blas_dscal(b13_6, integration_vectors.k13k6_out);
    gsl_blas_dscal(b13_7, integration_vectors.k13k7_out);
    gsl_blas_dscal(b13_8, integration_vectors.k13k8_out);
    gsl_blas_dscal(b13_9, integration_vectors.k13k9_out);
    gsl_blas_dscal(b13_10, integration_vectors.k13k10_out);
    vect_add(integration_vectors.k13_in, dof, 9,
	      integration_vectors.k13k1_out, integration_vectors.k13k4_out,
	      integration_vectors.k13k5_out, integration_vectors.k13k6_out,
	      integration_vectors.k13k7_out, integration_vectors.k13k8_out,
	      integration_vectors.k13k9_out,
	      integration_vectors.k13k10_out,
	      integration_vectors.k13k12_out);
    gsl_blas_dscal(h, integration_vectors.k13_in);
    gsl_vector_add(integration_vectors.k13_in, in_vec);

    (*f) (x0, integration_vectors.k13_in, integration_vectors.k13_out);

    //////////////
    // out_vec
    /////////////
    vect_add(integration_vectors.c_1_11_vec, dof, 2,
	      integration_vectors.k1_out, integration_vectors.k11_out);
    gsl_vector_memcpy(integration_vectors.c_6_vec,
		      integration_vectors.k6_out);
    vect_add(integration_vectors.c_7_8_vec, dof, 2,
	      integration_vectors.k7_out, integration_vectors.k8_out);
    vect_add(integration_vectors.c_9_10_vec, dof, 2,
	      integration_vectors.k9_out, integration_vectors.k10_out);
    gsl_blas_dscal(c_1_11, integration_vectors.c_1_11_vec);
    gsl_blas_dscal(c6, integration_vectors.c_6_vec);
    gsl_blas_dscal(c_7_8, integration_vectors.c_7_8_vec);
    gsl_blas_dscal(c_9_10, integration_vectors.c_9_10_vec);
    vect_add(integration_vectors.c_tot_vec, dof, 4,
	      integration_vectors.c_1_11_vec, integration_vectors.c_6_vec,
	      integration_vectors.c_7_8_vec,
	      integration_vectors.c_9_10_vec);

    gsl_blas_dscal(h, integration_vectors.c_tot_vec);
    gsl_vector_add(integration_vectors.c_tot_vec, in_vec);
    gsl_vector_memcpy(out_vec, integration_vectors.c_tot_vec);
    if(check==1){
        printf("\t output values: ");
        print_vec(out_vec);
    }
    //////////////
    // err_factor
    /////////////
    gsl_vector_memcpy(integration_vectors.ek12,
		      integration_vectors.k12_out);
    gsl_vector_memcpy(integration_vectors.ek13,
		      integration_vectors.k13_out);
    gsl_blas_dscal(-1, integration_vectors.ek12);
    gsl_blas_dscal(-1, integration_vectors.ek13);

    vect_add(integration_vectors.err_vec, dof, 4,
	      integration_vectors.k1_out, integration_vectors.k11_out,
	      integration_vectors.ek12, integration_vectors.ek13);
    //printf("error vector:\n");
    //print_vec(integration_vectors.err_vec);
    gsl_blas_dscal(err_factor, integration_vectors.err_vec);
    /*
    int i = 0;
    double maxbad = INFINITY;
    for(i; i < 1; i++){
        double val = gsl_vector_get(integration_vectors.err_vec, i);
        if(val < maxbad) maxbad = val;
    }
    printf("maxbad: %.2e \n", maxbad);
    */
    return gsl_blas_dnrm2(integration_vectors.err_vec);
}


 double single_stepRKV65(int (*f) (double, gsl_vector *, gsl_vector *),
                          int dof,
                          gsl_vector * in_vec,
                          gsl_vector * out_vec,
                          double x0,
                          double h,
                          int check){
    const double c_2 = 1./6.;
    const double c_3 = 4./15.;
    const double c_4 = 2./3.;
    const double c_5 = 5./6.;
    const double c_6 = 1.;
    const double c_7 = 1./15.;
    const double c_8 = 1.;


    const double a_21 = 1./6.;

    const double a_31 = 4./75.;
    const double a_32 = 16./75.;
    
    const double a_41 = 5./6.;
    const double a_42 = -8./3.;
    const double a_43 = 5./2.;
    
    const double a_51 = -165./64.;
    const double a_52 = 55./6.;
    const double a_53 = -425./64.;
    const double a_54 = 85./96.;

    const double a_61 = 12./5.;
    const double a_62 = -8.;
    const double a_63 = 4015./612;
    const double a_64 = -11./36.;
    const double a_65 = 88./255.;

    const double a_71 = -8263./15000.;
    const double a_72 = 124./75.;
    const double a_73 = -643./680.;
    const double a_74 = -81./250.;
    const double a_75 = 2484./10625.;
    const double a_76 = 0.;
    
    const double a_81 = 3501./1720.;
    const double a_82 = -300./43.;
    const double a_83 = 297275./52632.;
    const double a_84 = -319./2322.;
    const double a_85 = 24068./84065.;
    const double a_86 = 0.;
    const double a_87 = 3850./26703.;


    const double b_11 = 3./40.;
    const double b_12 = 0.;
    const double b_13 = 875./2244.;
    const double b_14 = 23./72.;
    const double b_15 = 264./1955.;
    const double b_16 = 0.;
    const double b_17 = 125./11592.;
    const double b_18 = 43./616.;
    
    const double b_21 = 13./160.;
    const double b_22 = 0.;
    const double b_23 = 2375./5984.;
    const double b_24 = 5./16.;
    const double b_25 = 12./85.;
    const double b_26 = 3./44.;
    const double b_27 = 0.;
    const double b_28 = 0.;

    ///////////
    // K1 /////
    ///////////
    gsl_vector_memcpy(integration_vectors.k1_in, in_vec);
    (*f) (x0, integration_vectors.k1_in, integration_vectors.k1_out);
    if(check==1) printf("\t integration_vectors.k1_out: "); //print_vec(integration_vectors.k1_out);

    //////////
    // k2
    /////////

    gsl_vector_memcpy(integration_vectors.k2_in,
		      integration_vectors.k1_out);
    gsl_blas_dscal(h * a_21, integration_vectors.k2_in);
    gsl_vector_add(integration_vectors.k2_in, in_vec);

    (*f) (x0 + h, integration_vectors.k2_in,
	  integration_vectors.k2_out);

    ///////
    // k3
    //////
    
    gsl_vector_memcpy(integration_vectors.k3k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k3k2_out,
		      integration_vectors.k2_out);
    gsl_blas_dscal(a_31, integration_vectors.k3k1_out);
    gsl_blas_dscal(a_32, integration_vectors.k3k2_out);
    vect_add(integration_vectors.k3_in, dof, 2,
	      integration_vectors.k3k1_out, integration_vectors.k3k2_out);
    gsl_blas_dscal(h, integration_vectors.k3_in);
    gsl_vector_add(integration_vectors.k3_in, in_vec);

    (*f) (x0 + h, integration_vectors.k3_in,
	  integration_vectors.k3_out);


    //////
    // k4
    /////
    gsl_vector_memcpy(integration_vectors.k4k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k4k2_out,
		      integration_vectors.k2_out);
    gsl_vector_memcpy(integration_vectors.k4k3_out,
              integration_vectors.k3_out);
    gsl_blas_dscal(a_41, integration_vectors.k4k1_out);
    gsl_blas_dscal(a_42, integration_vectors.k4k2_out);
    gsl_blas_dscal(a_43, integration_vectors.k4k3_out);
    vect_add(integration_vectors.k4_in, dof, 3,
	      integration_vectors.k4k1_out, integration_vectors.k4k2_out,
          integration_vectors.k4k3_out);
    gsl_blas_dscal(h, integration_vectors.k4_in);
    gsl_vector_add(integration_vectors.k4_in, in_vec);

    (*f) (x0 + h, integration_vectors.k4_in,
	  integration_vectors.k4_out);



    //////
    // k5
    //////
    gsl_vector_memcpy(integration_vectors.k5k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k5k2_out,
		      integration_vectors.k2_out);
    gsl_vector_memcpy(integration_vectors.k5k3_out,
              integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k5k4_out,
              integration_vectors.k4_out);
    gsl_blas_dscal(a_51, integration_vectors.k5k1_out);
    gsl_blas_dscal(a_52, integration_vectors.k5k2_out);
    gsl_blas_dscal(a_53, integration_vectors.k5k3_out);
    gsl_blas_dscal(a_54, integration_vectors.k5k4_out);
    vect_add(integration_vectors.k5_in, dof, 4,
	      integration_vectors.k5k1_out, integration_vectors.k5k2_out,
          integration_vectors.k5k3_out, integration_vectors.k5k4_out);
    gsl_blas_dscal(h, integration_vectors.k5_in);
    gsl_vector_add(integration_vectors.k5_in, in_vec);
    (*f) (x0 + h, integration_vectors.k5_in,
	  integration_vectors.k5_out);

    //////
    // k6
    ////////
    gsl_vector_memcpy(integration_vectors.k6k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k6k2_out,
		      integration_vectors.k2_out);
    gsl_vector_memcpy(integration_vectors.k6k3_out,
              integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k6k4_out,
              integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k6k5_out,
              integration_vectors.k5_out);
    gsl_blas_dscal(a_61, integration_vectors.k6k1_out);
    gsl_blas_dscal(a_62, integration_vectors.k6k2_out);
    gsl_blas_dscal(a_63, integration_vectors.k6k3_out);
    gsl_blas_dscal(a_64, integration_vectors.k6k4_out);
    gsl_blas_dscal(a_65, integration_vectors.k6k5_out);
    vect_add(integration_vectors.k6_in, dof, 5,
	      integration_vectors.k6k1_out, integration_vectors.k6k2_out,
          integration_vectors.k6k3_out, integration_vectors.k6k4_out,
          integration_vectors.k6k5_out);
    gsl_blas_dscal(h, integration_vectors.k6_in);
    gsl_vector_add(integration_vectors.k6_in, in_vec);
    (*f) (x0 + h, integration_vectors.k6_in,
          integration_vectors.k6_out);

    //////
    // k7
    ////// 
    gsl_vector_memcpy(integration_vectors.k7k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k7k2_out,
		      integration_vectors.k2_out);
    gsl_vector_memcpy(integration_vectors.k7k3_out,
              integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k7k4_out,
              integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k7k5_out,
              integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k7k6_out,
              integration_vectors.k6_out);
    gsl_blas_dscal(a_71, integration_vectors.k7k1_out);
    gsl_blas_dscal(a_72, integration_vectors.k7k2_out);
    gsl_blas_dscal(a_73, integration_vectors.k7k3_out);
    gsl_blas_dscal(a_74, integration_vectors.k7k4_out);
    gsl_blas_dscal(a_75, integration_vectors.k7k5_out);
    gsl_blas_dscal(a_76, integration_vectors.k7k6_out);
    vect_add(integration_vectors.k7_in, dof, 6,
	      integration_vectors.k7k1_out, integration_vectors.k7k2_out,
          integration_vectors.k7k3_out, integration_vectors.k7k4_out,
          integration_vectors.k7k5_out, integration_vectors.k7k6_out);
    gsl_blas_dscal(h, integration_vectors.k7_in);
    gsl_vector_add(integration_vectors.k7_in, in_vec);
    (*f) (x0 + h, integration_vectors.k7_in,
          integration_vectors.k7_out);

    //////
    // k8
    //////
    gsl_vector_memcpy(integration_vectors.k8k1_out,
		      integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k8k2_out,
		      integration_vectors.k2_out);
    gsl_vector_memcpy(integration_vectors.k8k3_out,
              integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k8k4_out,
              integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k8k5_out,
              integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k8k6_out,
              integration_vectors.k6_out); 
    gsl_vector_memcpy(integration_vectors.k8k7_out,
              integration_vectors.k7_out);
    gsl_blas_dscal(a_81, integration_vectors.k8k1_out);
    gsl_blas_dscal(a_82, integration_vectors.k8k2_out);
    gsl_blas_dscal(a_83, integration_vectors.k8k3_out);
    gsl_blas_dscal(a_84, integration_vectors.k8k4_out);
    gsl_blas_dscal(a_85, integration_vectors.k8k5_out);
    gsl_blas_dscal(a_86, integration_vectors.k8k6_out);
    gsl_blas_dscal(a_87, integration_vectors.k8k7_out);
    vect_add(integration_vectors.k8_in, dof, 7,
	      integration_vectors.k8k1_out, integration_vectors.k8k2_out,
          integration_vectors.k8k3_out, integration_vectors.k8k4_out,
          integration_vectors.k8k5_out, integration_vectors.k8k6_out, integration_vectors.k8k7_out);
    gsl_blas_dscal(h, integration_vectors.k8_in);
    gsl_vector_add(integration_vectors.k8_in, in_vec);
    (*f) (x0 + h, integration_vectors.k8_in,
          integration_vectors.k8_out);
    
    // the 5th order step will be stored in the k13k5 vector
    // the 6th order step will be stored in the k13k6 vector
    // the k12 vectors will be used for temporary storage
    
    gsl_vector_memcpy(integration_vectors.k12k1_out, integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k12k3_out, integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k12k4_out, integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k12k5_out, integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k12k7_out, integration_vectors.k7_out);
    gsl_vector_memcpy(integration_vectors.k12k8_out, integration_vectors.k8_out);
    gsl_blas_dscal(b_11, integration_vectors.k12k1_out);
    gsl_blas_dscal(b_13, integration_vectors.k12k3_out);
    gsl_blas_dscal(b_14, integration_vectors.k12k4_out);
    gsl_blas_dscal(b_15, integration_vectors.k12k5_out);
    gsl_blas_dscal(b_17, integration_vectors.k12k7_out);
    gsl_blas_dscal(b_18, integration_vectors.k12k8_out);
    vect_add(integration_vectors.k13k3_out, dof, 6,
            integration_vectors.k12k1_out, integration_vectors.k12k3_out,
            integration_vectors.k12k4_out, integration_vectors.k12k5_out,
            integration_vectors.k12k7_out, integration_vectors.k12k8_out);
    gsl_blas_dscal(h, integration_vectors.k13k3_out);
    vect_add(integration_vectors.k13k5_out, dof, 2,
            integration_vectors.k13k3_out, in_vec);

    gsl_vector_memcpy(out_vec, integration_vectors.k13k5_out);
    
    gsl_vector_memcpy(integration_vectors.k12k1_out, integration_vectors.k1_out);
    gsl_vector_memcpy(integration_vectors.k12k3_out, integration_vectors.k3_out);
    gsl_vector_memcpy(integration_vectors.k12k4_out, integration_vectors.k4_out);
    gsl_vector_memcpy(integration_vectors.k12k5_out, integration_vectors.k5_out);
    gsl_vector_memcpy(integration_vectors.k12k6_out, integration_vectors.k6_out);
    gsl_blas_dscal(b_21, integration_vectors.k12k1_out);
    gsl_blas_dscal(b_23, integration_vectors.k12k3_out);
    gsl_blas_dscal(b_24, integration_vectors.k12k4_out);
    gsl_blas_dscal(b_25, integration_vectors.k12k5_out);
    gsl_blas_dscal(b_26, integration_vectors.k12k6_out);
    vect_add(integration_vectors.k13k4_out, dof, 5,
            integration_vectors.k12k1_out, integration_vectors.k12k3_out,
            integration_vectors.k12k4_out, integration_vectors.k12k5_out,
            integration_vectors.k12k6_out);
    gsl_blas_dscal(h, integration_vectors.k13k4_out);
    vect_add(integration_vectors.k13k6_out, dof, 2,
            integration_vectors.k13k4_out, in_vec);
    //print_vec(integration_vectors.k13k6_out);
    //print_vec(integration_vectors.k13k5_out);
    gsl_blas_dscal(-1, integration_vectors.k13k6_out);
    vect_add(integration_vectors.k13k7_out, dof, 2,
            integration_vectors.k13k6_out, integration_vectors.k13k5_out);
    return gsl_blas_dnrm2(integration_vectors.k13k7_out);
 }

