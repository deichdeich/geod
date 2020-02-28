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

double oldh;
_Bool is_near = 0;
double dist1 = 100;
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
int has_crossed = 0;
int warn=0;
double prev_x;
double last_time_rec;
double hinit;
double h_rec = 999.;
int stepper2(int (*f) (double, gsl_vector *, gsl_vector *),
	     int dof,
	     gsl_vector * in_state_vec,
	     gsl_vector * out_state_vec,
	     double x,
	     double h,
	     double xmax,
	     double tolerance){
    hinit = h;
    double err;
    double scale;
    int step = 0;
    cur_state = gsl_vector_calloc(dof);
    gsl_vector_memcpy(integration_vectors.temp_in_vec, in_state_vec);
    gsl_vector_memcpy(integration_vectors.temp_out_vec, out_state_vec);
    gsl_vector_memcpy(cur_state, in_state_vec);
    gsl_vector_memcpy(integration_vectors.prev_state, in_state_vec);
    h_rec = h;
    while (x < xmax+h){
        int attempts = 0;
        int err_val = 0;

        err = fabs(int_func(f,
                        dof,
                        cur_state,
                        integration_vectors.temp_out_vec,
                        x, h,0));
        err_val = check_err(tolerance);
        
        while (err_val == 0){
            scale = pow(0.8 * h * (1/tolerance),(1./7.));
            scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
            h *= scale;
            err = fabs(int_func(f,
                        dof,
                        cur_state,
                        integration_vectors.temp_out_vec,
                        x, h, 0));
            err_val = check_err(tolerance);
            attempts += 1;
        
            if(attempts > 100){
                printf(">100 attempts to get the error below tolerance");
                return -2;
                }
            
        }
        if(lines_written > 1 && warn!=2){
            //print_vec(integration_vectors.temp_out_vec);
            warn+=1;
        }
        oldh = h;
        if (make_section){
            //h = poincare_check2(cur_state, h);
            poincare_check_bisect(f, cur_state, h, dof);
            gsl_vector_memcpy(cur_state, integration_vectors.prev_state);
            //if (h!=h_rec) printf("%f\n",h);h_rec = h;
            if (poincare_yes){
                add_to_history = 1;
                has_crossed = 0;
                }
            }


        if(add_to_history){
           //populate_history(dof, lines_written, x, cur_state);
           //printf("%f\n",x);
           //print_vec(cur_state);
           // this is fucking up the order of things in the history array
            if (x == 0. || !make_section){
                populate_history(dof, lines_written, x, cur_state);
            }
            else if (x > (last_time_rec + 10.)){
                populate_history(dof, lines_written, prev_x, integration_vectors.prev_state);
                poincare_yes = 0;
                last_time_rec = x;
            }
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
        prev_x = x;
        x += h;
        
    }	
    write_history((lines_written-1)%hist_len);
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

/* if the particle is near (10x tol), then:
 *      -slow down timestep by some factor (start with 10)
 *      if the particle is close (within tol), then:
 *          -start keeping track of how close it's getting
 *          -as soon as the distance from the surface increases,
 *          take the previous state vector.
 *
 *
 */
int orbit_num = 0;
_Bool is_it_closer = 0;

void poincare_check_bisect(int (*f) (double, gsl_vector *, gsl_vector *), gsl_vector * state, double h, int dof){
    
    double r, th, Pth;
    r = gsl_vector_get(state, 0);
    th = gsl_vector_get(state, 2);
    Pth = gsl_vector_get(state, 3);

    double prev_r, prev_th, prev_Pth;
    prev_r = gsl_vector_get(integration_vectors.prev_state, 0);
    prev_th = gsl_vector_get(integration_vectors.prev_state, 2);
    prev_Pth = gsl_vector_get(integration_vectors.prev_state, 3);
    
    gsl_vector *temp_vec = gsl_vector_calloc(dof);
    
    th = fmod(th, TPI);
    prev_th = fmod(prev_th, TPI);

    double thdist, prev_thdist;
    thdist = th - poincare_condition[0];
    prev_thdist = prev_th - poincare_condition[0];

    
    if ((Pth / poincare_condition[1]) > 0){
        if ((thdist * prev_thdist) < 0 && fabs(prev_thdist) > poincare_tolerance){
            //printf("pc sees: ");
            //print_vec(state);
            bisect(f, prev_thdist, h, temp_vec, dof);
            poincare_yes = 1;
            has_crossed = 1;
            gsl_vector_memcpy(integration_vectors.prev_state, temp_vec);
            //printf("pc out:");
            //print_vec(integration_vectors.prev_state);
        }
        else gsl_vector_memcpy(integration_vectors.prev_state, state);
    }
    else gsl_vector_memcpy(integration_vectors.prev_state, state);
    gsl_vector_free(temp_vec);
}

void bisect(int (*f) (double, gsl_vector *, gsl_vector *), double prev_dist, double h, gsl_vector *bisect_out, int dof){
    gsl_vector *bisect_temp_vec1 = gsl_vector_calloc(dof);
    gsl_vector *bisect_temp_vec2 = gsl_vector_calloc(dof);
    gsl_vector *bisect_temp_vec3 = gsl_vector_calloc(dof);
    gsl_vector_memcpy(bisect_temp_vec3, integration_vectors.prev_state);
    double prev_th = gsl_vector_get(integration_vectors.prev_state, 2);
    double h1 = -h;
    double h2 = h;
    double err1, err2, err3;
    double dist1, dist2, dist3;
    double th1, th2, th3;
    double c;
    dist3 = prev_dist;
    //printf("Surface crossing detected, doing a bisection\n");
    //printf("\tcurrent theta: %0.10f (%e) \n", prev_th, fabs(dist3));
    int n = 0;
    while(fabs(dist3) > poincare_tolerance){  
        c = (h1 + h2) / 2.;
       // printf("\tc: %f\n",c);
       // printf("\th1: %0.10f\n",h1);
        
        integration_vector_zero(dof);
        err1 = fabs(int_func(f,
                        dof,
                        integration_vectors.prev_state,
                        integration_vectors.temp_out_vec,
                        0, h1, 0));
        gsl_vector_memcpy(bisect_temp_vec1, integration_vectors.temp_out_vec); 
        
        integration_vector_zero(dof);
        err3 = fabs(int_func(f,
                        dof,
                        integration_vectors.prev_state,
                        integration_vectors.temp_out_vec,
                        0, c, 0));
        gsl_vector_memcpy(bisect_temp_vec3, integration_vectors.temp_out_vec);

        th1 = fmod(gsl_vector_get(bisect_temp_vec1, 2),TPI);
        th3 = fmod(gsl_vector_get(bisect_temp_vec3, 2),TPI);
        
        dist1 = th1 - poincare_condition[0];
        dist3 = th3 - poincare_condition[0];
        //printf("\tdist1, dist3 = %.10e, %.10e\n", dist1, dist3);

        if (dist3 / dist1 > 0){
            h1 = c;
        }
        
        else if(dist3 / dist1 < 0){
            h2 = c;
        }
        
        //printf("\tcurrent theta: %0.10f (%e) \n", th3, fabs(dist3));
        n+=1;
        if(n == 100){
            //printf("\tunable to converge\n");
            break;
        }
    }

   //printf("xxx");
    //print_vec(integration_vectors.temp_out_vec);
    //print_vec(bisect_temp_vec3);
    //printf("\taccepted theta: %0.10f (%e) \n", th3, fabs(dist3));
    gsl_vector_memcpy(bisect_out, integration_vectors.temp_out_vec);
    gsl_vector_free(bisect_temp_vec3);
    gsl_vector_free(bisect_temp_vec2);
    gsl_vector_free(bisect_temp_vec1);
}


double poincare_check2(gsl_vector * state, double h){
    double r, prev_r, th, Pth;
    
    prev_r = gsl_vector_get(integration_vectors.prev_state, 0);
    r = gsl_vector_get(state, 0);

    th = gsl_vector_get(state, 2);
    Pth = gsl_vector_get(state, 3);
    
    // the distance of the particle from the surface
    double dist2 = fmod(fabs(th - poincare_condition[0]),TPI);
    
    double rdist = fabs(r-prev_r);


    // it's getting close and is going the correct direction
    if (dist2  < (poincare_tolerance * SLOWDOWN_DIST) && (Pth / poincare_condition[1]) > 0 && is_near == 0){
        h /= HIRES_SCALE_FACTOR;
        is_near = 1;
    }
    else if (dist2  > (poincare_tolerance * SLOWDOWN_DIST)){
        is_near = 0;
        h = hinit;
        }
    // now it's within tolerance and going the correct direction
    if (dist2 < poincare_tolerance && (Pth / poincare_condition[1]) > 0){
        // while the particle is getting closer to the plane, update the prev_state vector
        if (dist2 < dist1){
            is_it_closer = 1;
            gsl_vector_memcpy(integration_vectors.prev_state, state);
            //printf("orbit %d:\n",orbit_num);
            //printf("  %f\n",dist1);
            //printf("  %f\n\n",dist2);
            dist1 = dist2;
        }
        // now that it's closest approach was the previous timestep, allow the state to be printed
        else if (dist2 > dist1 && poincare_yes == 0 && rdist > (r * h * 100)){
            //printf("%d\n",orbit_num);
            is_it_closer = 0;
            poincare_yes = 1;
            // return the timestep to the size it was before it approached the surface
            dist1 = 100;
            h = hinit;
            orbit_num += 1.;
     }
    }
    return h;
}


void poincare_check(gsl_vector * state){
    double r, th, Pth;
    r = gsl_vector_get(state, 0);
    th = gsl_vector_get(state, 2);
    Pth = gsl_vector_get(state, 3);

    /* the first if checks to see if Th1 is within the tolerance of the desired Th1 */
    if (fmod(fabs(th - poincare_condition[0]), TPI) < poincare_tolerance )  {
        /* the second checks to see if Th1_d and the desired Th1_d are on the 
        same side of 0 */
        if ((Pth / poincare_condition[1]) > 0){
            poincare_yes = 1;
        }
    }
}

