#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "single_step.h"
#include "integration_vectors.h"

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

// testing the rkf78 implementation in single_step.c
// need to write:
//      - exponential derivative function
//      - timestepping function
//      - main function

int exponential(double t, gsl_vector *in_vec, gsl_vector *out_vec);
void print_vec(gsl_vector * vec, int N);
int *vect_add(gsl_vector * sum, int vect_len, int count, ...);

int timestep(int (*f) (double, gsl_vector *, gsl_vector *),
             int dof,
             gsl_vector * in_vec,
             gsl_vector * out_vec,
             double t0,
             double tmax,
             double h_init,
             double tol,
             int step_check,
             int time_check);

int exponential(double t, gsl_vector *in_vec, gsl_vector *out_vec){
    gsl_vector_memcpy(out_vec, in_vec);

    return 0;
}

int sine(double t, gsl_vector *in_vec, gsl_vector *out_vec){
    double x1 = gsl_vector_get(in_vec, 0);
    double x2 = gsl_vector_get(in_vec, 1);

    x1 *= -1;

    gsl_vector_set(out_vec, 0, x2);
    gsl_vector_set(out_vec, 1, x1);

    return 0;
}

void print_vec(gsl_vector * vec, int N){
     int i;
     printf("[");
     for (i = 0; i < N; i++) {
         printf("%0.10e ",gsl_vector_get(vec, i));
     }
     printf("\b]\n");
}

int *vect_add(gsl_vector * sum, int vect_len, int count, ...){   
     va_list p; 
     int i;
     gsl_vector *tmp_vec;
     va_start(p, count);
         
     for (i = 0; i < count; i++) {
         tmp_vec = va_arg(p, gsl_vector *);
         gsl_blas_daxpy(1, tmp_vec, sum);
     }
     va_end(p);
     return 0;
}

int timestep(int (*f) (double, gsl_vector *, gsl_vector *),
             int dof,
             gsl_vector * in_vec,
             gsl_vector * out_vec,
             double t0,
             double tmax,
             double h_init,
             double tol,
             int step_check,
             int time_check){
    // timestep the function f from t0 to tmax

    double h = h_init;
    double t = t0;
    gsl_vector *temp_in_vec = gsl_vector_calloc(dof);
    gsl_vector_memcpy(temp_in_vec, in_vec);

    gsl_vector *temp_out_vec = gsl_vector_calloc(dof);

    double err_i = 999;
    double R = 999;
    double delta = 999;

    int iters = 0;

    while(t < tmax){
        h = (h < (tmax - t)) ? h : (tmax - t);
        iters += 1;
        //if(iters > 10){
        //    printf("iteration threshold reached\n");
        //    break;
        //}
         
        err_i = single_stepRKF78(sine,
                               dof,
                               temp_in_vec,
                               temp_out_vec,
                               t,
                               h,
                               step_check);
        integration_vector_zero(dof);
        R = err_i;
        delta = 0.95 * pow((tol/R), (1. / 8.));
        //delta = min(max(delta, 0.125), 4); 
        if (time_check == 1){
            printf("\t temporary in vec:\n");
            printf("\t\t\t ");
            print_vec(temp_in_vec, dof);
            printf("\t\t\t ");
            printf("current time: %.2e ", t);
            printf("dt: %.2e, ", h);
            printf("err: %.2e, ", err_i);
            printf("R: %.2e, ", R);
            printf("delta: %.2e,\n\n", delta);
            printf("\t temporary out vec:\n");
            printf("\t\t\t ");
            print_vec(temp_out_vec, dof);
            printf("\n");
        }

        if (R < tol){
            t = t + h;
            gsl_vector_memcpy(temp_in_vec, temp_out_vec);
            if(time_check == 1){
                printf("\t new starting point:\n\t\t\t");
                print_vec(temp_in_vec, dof);
            }
        }
        if(isinf(delta)) delta = 10;
        h = h * delta;

    }
    gsl_vector_memcpy(out_vec, temp_out_vec);

    return 0;
}

int main(){
    gsl_vector *x_0 = gsl_vector_calloc(1);
    gsl_vector_set(x_0, 0, 1.);
    gsl_vector *y_0 = gsl_vector_calloc(1);
    integration_vector_init(2);
    // double single_stepRKF78(int (*f) (double, gsl_vector *, gsl_vector *),
    //                         int dof,
    //                         gsl_vector * in_vec,
    //                         gsl_vector * out_vec,
    //                         double x0,
    //                         double h,
    //                         int check){

    int dof = 1;
    double h_init = 0.1;
    int check = 0;
    double t = 0.;
    gsl_vector *init_state = gsl_vector_calloc(2);
    gsl_vector_set(init_state, 0, 0.);
    gsl_vector_set(init_state, 1, 1.);
    gsl_vector *out_state = gsl_vector_calloc(2);
    double err = single_stepRKF78(sine,
                                  2,
                                  init_state,
                                  out_state,
                                  t,
                                  h_init,
                                  0);
    printf("\n\n\n");
    printf("############\n");
    printf("SINGLE STEP TEST\n");
    printf("############\n");
    printf("\t input vec:\n");
    printf("\t\t");
    print_vec(init_state, 2);
    printf("\t\t");
    printf("dt: %.2e\n", h_init);
    printf("\t output vec:\n");
    printf("\t\t");
    print_vec(out_state, 2);
    printf("\t error:\n");
    printf("\t\t");
    printf("%.4e\n", err);

    integration_vector_zero(2);

// int timestep(int (*f) (double, gsl_vector *, gsl_vector *),
//              int dof,
//              gsl_vector * in_vec,
//              gsl_vector * out_vec,
//              double t0,
//              double tmax,
//              double h_init,
//              double tol,
//              int step_check,
//              int time_check){

    printf("\n\n\n");
    printf("############\n");
    printf("MULTI STEP TEST\n");
    printf("############");
    printf("\n");
    double tm = 1;
    
    timestep(sine,
             2,
             init_state,
             out_state,
             t,
             tm,
             h_init,
             1e-15,
             0,
             1);

    double ans = gsl_vector_get(out_state, 0);
    double ans_0 = sin(tm);
    printf("ans: %.4e,  pct err: %.4e\n", ans, abs(ans - ans_0));
    
}
