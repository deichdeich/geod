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

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )
int nunc = 0;
double oldh;
_Bool is_near = 0;
double dist1 = 100;
/*
 * steppers2.c:
 * written in May 2020 in order to try to fix
 * the rkf78 method
 *
 */
int check_err(double tol, int verbose){
    //printf("check_err: ");
    //print_vec(integration_vectors.err_vec);
    //printf("tol: %0.10e\n", tol);
    int good_err = 1;
    //printf("\t%d\n",good_err);
    int i;
    //printf("\terrors:\n");
    for(i = 0; i < 6; i++){
        if(fabs(gsl_vector_get(integration_vectors.err_vec, i)) > tol){
            if(verbose == 1){
                printf("\t bad error:  %0.10e\n", gsl_vector_get(integration_vectors.err_vec, i));
            }
            good_err = 0;
        }
    }
    //printf("\t%d\n",good_err);
    return good_err;
}

void get_scale(scale, range, tol){
    double newscale = 0;
    int i;
    for(i = 0; i < 6; i++){
        double err = fabs(gsl_vector_get(integration_vectors.err_vec, i));
        double yi = gsl_vector_get(integration_vectors.temp_out_vec, i);
        double tempscale = 0.8 * pow(fabs(err * yi / (err * range)), 1./7.);
        if (tempscale > newscale){
            newscale = tempscale;
        }
    }
    scale = newscale;
}


int has_crossed = 0;
double prev_r_acc_dir;
double cur_r_acc_dir;
double prev_x;
double theta;
int passpi = 1;
double prev_pr = 0;
double cur_pr;
double cur_pd, prev_pd;
double last_time_rec;
double hinit;
gsl_vector *acc_state;
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
    int step = 0;
    cur_state = gsl_vector_calloc(dof);
    acc_state = gsl_vector_calloc(dof);
    gsl_vector_memcpy(integration_vectors.temp_out_vec, out_state_vec);
    gsl_vector_memcpy(cur_state, in_state_vec);
    gsl_vector_memcpy(integration_vectors.prev_state, in_state_vec);
    double maxerr = -INFINITY;
    double minerr = INFINITY;
    double maxdt = -INFINITY;
    double mindt = INFINITY;
    double maxdttime = 0;
    double maxerrtime = 0;
    double minerrtime = 0;
    int attempts = 0;
    f(0., cur_state, acc_state);
    prev_r_acc_dir = gsl_vector_get(acc_state, 1);
    prev_pr = 0.;
    //print_vec(acc_state);
    //printf("ipr: %f\n",gsl_vector_get(cur_state, 1));
    //printf("prev_pr: %f\n",prev_pr);
    //printf("0: prev_r_acc_dir: %f\n",prev_r_acc_dir);
    while (x < xmax+h){
        //printf("1: prev_r_acc_dir: %f\n",prev_r_acc_dir);
        h = fabs(h);
        //print_vec(integration_vectors.temp_out_vec);
        //printf("2: prev_r_acc_dir: %f\n",prev_r_acc_dir);
        err = int_func(f,
                        dof,
                        cur_state,
                        integration_vectors.temp_out_vec,
                        x, h, 0);
        //printf("3:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        //print_vec(integration_vectors.temp_out_vec);
        //printf("4:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        if(err < minerr){
            minerr = err;
            minerrtime = x;
        }
        //printf("5:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        else if(err > maxerr){
            maxerr = err;
            maxerrtime = x;
        }
        //printf("6:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        double R = err/h;
        //printf("7:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        gsl_vector *temp_vec = gsl_vector_calloc(6);
        //printf("8:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        gsl_vector_memcpy(temp_vec, integration_vectors.err_vec);
        //printf("9:prev_r_acc_dir: %f\n",prev_r_acc_dir);
        gsl_blas_dscal(h, temp_vec);
        //printf("10: prev_r_acc_dir: %f\n",prev_r_acc_dir);break;
        double S = (err == 0) ? 3 : pow((tolerance/R), (1. / 5.));
        //printf("\n time: %.2e  dt: %.2e  err: %.2e  R: %.2e  delta: %.2e, next dt: %.2e\n\n", x, h, err, R, S, h*S);
        //int ii;
        //for(ii = 0; ii < 6; ii++){
        //    printf("%.3e, ",gsl_vector_get(temp_vec, ii)/gsl_vector_get(cur_state, ii));
        //}
        //printf("\n");
        if(S >=1){//(R < tolerance){
            // check for a section crossing & bisect to get it accurate 
            if (make_section){
                h = poincare_check_bisect(f, cur_state, h, dof, x);
                gsl_vector_memcpy(cur_state, integration_vectors.prev_state);
                theta = fmod(gsl_vector_get(cur_state, 4), TPI);
                //printf("%f\n",theta);
                if (cos(theta) < -1 * 0.9) {
                    passpi = 1;
                 //   printf("yes\n");
                }
                if (poincare_yes){
                   // printf("%f\n",x);
                    add_to_history = 1;
                    has_crossed = 0;
                    }
                }


            if(add_to_history){
                int ii = 1;
                if (x == 0. || !make_section){
                    populate_history(dof, lines_written, x, cur_state);
                    lines_written++;    
                }
                else if (x>0 && make_section){
                   // gsl_vector_memcpy(acc_vec, integration_vectors.prev_state);
                    //f(0, integration_vectors.prev_state, acc_state);
                    //cur_r_acc_dir = gsl_vector_get(acc_state, 1);
                    //printf("prev_r_acc_dir: %f\n",prev_r_acc_dir);
                    if (x > (prev_pr + 2)){//((cur_r_acc_dir / prev_r_acc_dir) < 0){
                        //printf("%f: ", prev_x);
                        //printf("stepper is giving: ");
                       // print_vec(integration_vectors.prev_state);
                        poincare_yes = 0;

                        populate_history(dof, lines_written, prev_x, integration_vectors.prev_state);
                        last_time_rec = x;
                        lines_written++;
                        //printf("%f: ", prev_x);
                        //printf("last hist line is: [");
                        //for(ii; ii < 7; ii++){
                        //   printf("%0.10e ",history[(lines_written-1)%hist_len][ii]);
                        //}
                        //printf("\b]\n\n");
                        prev_pr = x;
                        prev_r_acc_dir = cur_r_acc_dir;
                        passpi = 0;
                        }
                }
                if(make_section){
                    add_to_history = 0;
                    poincare_yes = 0;
                }
            }
            
            // clean up integration_vector struct, copy out state to file
            gsl_vector_memcpy(cur_state,
                              integration_vectors.temp_out_vec);
            
            step += 1;
            
            // advance timestep
            prev_x = x;
            x += h;
            h *= min(S, 2);   
       
        if(h < mindt){
            mindt = h;
            //printf("new min dt: %.10e\n", h);
        }
        if(h > maxdt){
            maxdt = h;
            maxdttime = x;
        }
        }
        if(x == 0.){
            //printf("backwards, dt = %0.05e\n", h);
        }
        else if(S < 1){//(R > tolerance){
            h *= max(S, 1/2);
        }
        //if(h > 1) h = 1/h;
        integration_vector_zero(dof);
        /*
        h *= S;
        if isnan(h){
            //printf("%.2e\n", h);
            h = 0.1;
        }
       */ 
        if(h > 1){
            h = 1/h;
        }
        /*
        attempts += 1;
        if(attempts > 50){
            break;
        }
        */
    }	
    write_history((lines_written-1)%hist_len);
    integration_vector_free();
    printf("minerr: %.2e at %.2f\n", minerr, minerrtime);
    printf("maxerr: %.2e at %.2f\n", maxerr, maxerrtime);
    printf("maxdt: %.2e at %.2f\n", maxdt, maxdttime);
    printf("final dt: %.2f \n", h);
    printf("unconverged bisections: %d\n", nunc);
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


double poincare_check_bisect(int (*f) (double, gsl_vector *, gsl_vector *), gsl_vector * state, double h, int dof, double x){
        
    gsl_vector *temp_bisect_vec = gsl_vector_calloc(dof);
    
    double r, th, Pth;
    r = gsl_vector_get(state, 0);
    th = gsl_vector_get(state, 2);
    Pth = gsl_vector_get(state, 3);

    double prev_r, prev_th, prev_Pth;
    prev_r = gsl_vector_get(integration_vectors.prev_state, 0);
    prev_th = gsl_vector_get(integration_vectors.prev_state, 2);
    prev_Pth = gsl_vector_get(integration_vectors.prev_state, 3);
    
    th = fmod(th, TPI);
    prev_th = fmod(prev_th, TPI);

    double thdist, prev_thdist;
    thdist = th - poincare_condition[0];
    prev_thdist = prev_th - poincare_condition[0];

    double newh = h;
    if ((Pth / poincare_condition[1]) > 0){
        if ((thdist * prev_thdist) < 0){
            //printf("%f: ", x);
            //print_vec(state);
            newh = bisect(f, prev_thdist, h, temp_bisect_vec, dof, x);
            poincare_yes = 1;
            has_crossed = 1;
            gsl_vector_memcpy(integration_vectors.prev_state, temp_bisect_vec);
            //printf("pc out:");
            //print_vec(integration_vectors.prev_state);
        }
        else gsl_vector_memcpy(integration_vectors.prev_state, state);
    }
    else gsl_vector_memcpy(integration_vectors.prev_state, state);

    gsl_vector_free(temp_bisect_vec);
    return(newh);
}

double bisect(int (*f) (double, gsl_vector *, gsl_vector *), double prev_dist, double h, gsl_vector *bisect_out, int dof, double x){
    gsl_vector *bisect_temp_vec1 = gsl_vector_calloc(dof);
    gsl_vector *bisect_temp_vec2 = gsl_vector_calloc(dof);
    gsl_vector *bisect_temp_vec3 = gsl_vector_calloc(dof);
    gsl_vector_memcpy(bisect_temp_vec3, integration_vectors.prev_state);
    double prev_th = gsl_vector_get(integration_vectors.prev_state, 2);
    //printf("%f: %f\n",x,prev_th - 1.5707963268);
    double h1 = 0;
    double h2 = 5 * h;
    //printf("h1: %0.5e, h2: %0.5e\n", h1, h2);
    double err1, err2, err3;
    double dist1 = INFINITY;
    double dist3 = INFINITY;
    double th1, th2, th3;
    double c;
    //dist3 = prev_dist;
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
        //if(fabs(dist3) < poincare_tolerance) printf("dist3 is within bounds: %.5e\n",dist3);
        //if(fabs(dist1) < poincare_tolerance) printf("dist1 is within bounds: %.5e\n",dist1);
       // printf("\tcurrent theta: %0.10f (%e) \n", th3, fabs(dist3));
        n+=1;
        if(n == 10000){
            printf("failed to converge\n");
            nunc += 1;
            break;
        }
    }
    th3 = fmod(gsl_vector_get(integration_vectors.temp_out_vec, 2),TPI);
    dist3 = th3 - poincare_condition[0];
    if (fabs(dist3) > poincare_tolerance){
        //printf("bad bisect: %.5e\n", dist3);
        nunc++;
    }
   //printf("xxx");
    //print_vec(integration_vectors.temp_out_vec);
    //print_vec(bisect_temp_vec3);
    //printf("\taccepted theta: %0.10f (%e), h1: %0.5e, h2: %0.5e, t = %0.5f\n\n", th3, fabs(dist3), h1, h2, x);
    gsl_vector_memcpy(bisect_out, integration_vectors.temp_out_vec);
    //printf("\t%f\n",gsl_vector_get(bisect_out, 3));
    gsl_vector_free(bisect_temp_vec3);
    gsl_vector_free(bisect_temp_vec2);
    gsl_vector_free(bisect_temp_vec1);
    return(c);
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
            if (r < 10){
                poincare_yes = 1;
        }
    }
}
}

