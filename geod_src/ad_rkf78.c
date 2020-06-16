/*
ad_gr_rkf78.c:
Alex Deich's implementation of an RKF7(8) integrator.
Now with arbitrary degrees of freedom for generalization beyond pendula.
Date: March, 2020
adeich2@illinois.edu
*/
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

_Bool add_to_history = 1.;
_Bool poincare_yes = 0.;
int lines_written = 0;


/* poincare section flags and condition */



////////////////////////////////////////////////////////////////////
//  rkf78: an embedded Runge-Kutta-Felberg method
//
//  It's an adaptive time-step method.  My implementation will take a grid of timesteps
// which the integrator will adaptively step through to achieve the desired tolerance
// at each timestep.  So, even though the time is initially given in a rigid grid (and
// ultimately returned in the same rigid grid) getting between those grids takes a variable
// number of timesteps to ensure the desired precision.
//
// Arguments
//    *f            pointer to the function which gives the derivative of the system
//    xmin          The starting point along the integration coordinate
//    xmax          The ending point
//    init_state    The inish condish at xmin: [Th1, Th2, Th1_d, Th2_d]
//    hist_len      The total number of timesteps (and length of the history array)
//    history       An array of size (5,N) containing all integration steps, edited in-place
//    tolerance     The tolerance for each step
//    h             Initial step size

// current bottleneck is in file i/o which is currently done naively.  see
// this discussion: https://stackoverflow.com/questions/41210227/fastest-way-to-write-integer-to-file-in-c
/////////////////////////////////////////////////////////////////////////
int rkf78(double xmin,
	      double xmax,
    	  int dof,
	      double init_state[dof],
	      double tol,
    	  double h,
	      _Bool init_make_section,
    	  double init_poincare_condition[2],
	      double init_poincare_tolerance,
	      char init_filename[],
	      char integrator[])
{

    double hpt;
    if(strcmp("rkf78", integrator) == 0){
        int_func = single_stepRKF78;
    }
    else if(strcmp("rkv65", integrator) == 0){
        int_func = single_stepRKV65;
    }
    else if(strcmp("rk4", integrator) == 0){
        //int_func = single_stepRK4;
        printf("RK4 integration disabled at the moment\n");
    }
    else{
        printf("Bad integrator specification: %s\n  Options: rk4, rkf78\n", integrator);
        return -2;
    }
    
    
    filename = init_filename;
    clear_file();

    integration_vector_init(dof);

    make_section = init_make_section;
    poincare_condition[0] = init_poincare_condition[0];
    poincare_condition[1] = init_poincare_condition[1];
    poincare_tolerance = init_poincare_tolerance;
    
    if(make_section){
        add_to_history = 1;
    }
    /* initialize initial state */
    gsl_vector *init_state_vec = gsl_vector_calloc(dof);
    
    /* load the array from python into the gsl_vector */
    arr2vec(dof, init_state, init_state_vec);
    gsl_vector_memcpy(integration_vectors.temp_in_vec, init_state_vec);

    /* initialize the final state */
    gsl_vector *out_state_vec = gsl_vector_calloc(dof);
    
    /* step through all the times */
	stepper2(EOM, dof, init_state_vec, out_state_vec, xmin, h, xmax, tol);
	
    poincare_yes = 0;
    return 0;
}




/////////////////////////
//   Helper functions  //
/////////////////////////

/////////////////////////////////////////////////////
// vect_add adds arbitrary numbers of gsl_vectors 
// this is used to calculate the k1...k13
// Arguments:
//        sum: A pointer to the gsl_vector which will hold the sum of the vectors, gsl_vector
//        vect_len: the length of the vectors, int
//        count: the number of vectors, int
//        ...: the vectors to be added, gsl_vector
// Returns:
//        0 if successful
//        No failure case yet
/////////////////////////////////////////////////////
int *vect_add(gsl_vector * sum, int vect_len, int count, ...)
{
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




//////////////////////
// print_vec prints a gsl_vector
// used for troubleshooting
//
// Argument:
//      vec: a gsl_vector pointer
//
//////////////////////
void print_vec(gsl_vector * vec)
{
    int i;
    printf("[");
    for (i = 0; i < 6; i++) {
        printf("%0.10e ",gsl_vector_get(vec, i));
    }
    printf("\b]\n");
}

/////////////////////
// populate_history populates the history array
// which is the output of the integration.
// The history array has a shape nsteps * (1 + PhaseSpaceSize)
// for some arbitrary size of a phase space.
// the 0th column is used for the timestamp of that step.
//
// Arguments:
//      hist_len: The number of time steps, integer
//      history: The history array to be populated, double
//      step: The step being populated, int
//      clock: The value of the clock at the current step
//      state: The state at the clock time, gsl_vector
//
// Returns:
//      0 if successful
//      No failure case yet.
///////////////////

double lastclock = INFINITY;
int populate_history(int dof,
                     int step,
                     double clock,
                     gsl_vector * state){
    int line = step % hist_len;
    history[line][0] = clock;
    int i;
   // printf("%f: populate_history : ", clock);
   // print_vec(state);
    //printf("");
    lastclock = clock;
    for (i = 1; i < dof+1; i++) {
	    history[line][i] = gsl_vector_get(state, i - 1);
    }
    //printf("line %d\n", line);
    if(line == (hist_len-1)){
        //printf("calling it...\n");
        write_history(line);
    }
    return 0;
}

/////////////////
// arr2vec converts a regular ol' array to a gsl_vector
//
// Arguments:
//      len: The length of the array and vector, integer
//      in_arr:  The array being converted, double
//      out_vec:  A pointer to the gsl_vector being populated, gsl_vector
//
// Returns:
//      0 if successful
//      No failure case yet
////////////////
int arr2vec(int len, double in_arr[len], gsl_vector * out_vec)
{
    int i;
    for (i = 0; i < len; i++) {
	    gsl_vector_set(out_vec, i, in_arr[i]);
    }
    return 0;
}



// ./ad_gr_rkf78_wtf.e 0. 1e5 1e-10 1e-1 blah/out_ms.bin
//command line syntax: ./ad_gr_rkf78_wtf.c end_time mass spin energy r pr th pth ph Jz filename
// argc = 11
clear_history();
int main(int argc, char *argv[])
{
	double t0 = atof(argv[1]);
	double t1 = atof(argv[2]);
	double tol = 1e-15;
	double h = atof(argv[3]);
	double ipc[2] = {1.5707963268, 1.};
	double ipt = atof(argv[4]);
	_Bool ims = atoi(argv[5]);
    double test_vec[6] = {atof(argv[6]),
                          atof(argv[7]),
                          atof(argv[8]),
                          atof(argv[9]),
                          atof(argv[10]),
                          atof(argv[11])};
    char *ifn = argv[12];
	rkf78(t0,
          t1,
          6,
          test_vec,
          tol,
          h,
          ims,
          ipc,
          ipt,
          ifn,
          "rkv65");

return 0.;
}
