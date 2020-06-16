#include <stdio.h>
#include <stdlib.h>
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



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/////////
/////////  Naming conventions
/////////  Coordinates
/////////  o Coordinates with roman letters are written as the lower-case
/////////    version of that letter
/////////    Ex: the radial coordinate is written as "r"
/////////  o Coordinates with greek letters are written as the two-letter abbreviation
/////////    with the first letter capitalized.
/////////    Ex: theta becomes "Th", Ph becomes "Ph".
/////////
/////////  Momenta and time derivatives
/////////  o Momenta are written by prepending "P_" to the beginning of the coordinate.
/////////    Ex: The Ph momentum becomes "P_Ph"
/////////  o Time derivatives are written by appending "_d" to the end of the momentum
/////////    or coordinate
/////////    Ex: The time derivative of the r coordinate becomes "r_d",
/////////        the time derivative of the Ph momentum becomes "P_Ph_d"
///////// 
/////////  Coupling parameters
/////////  o Coupling/expansion parameters are written as the name of their full greek
/////////    letter, with appropriate capitalization
/////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int double_pendulum(double t, gsl_vector * in_state, gsl_vector * out_state){
    double Th1, Th1_d, Th2, Th2_d;
    Th1 = gsl_vector_get(in_state, 0);
    Th1_d = gsl_vector_get(in_state, 1);
    Th2 = gsl_vector_get(in_state, 2);
    Th2_d = gsl_vector_get(in_state, 3);

    double a1, a2, a3, a4, Th1_dd;
    double b1, b2, b3, b4, Th2_dd;

    a1 = GRAVITY * (sin(Th2) * cos(Th1 - Th2) - 2 * sin(Th1));
    a2 = -(Th2_d * Th2_d + Th1_d * Th1_d * cos(Th1 - Th2));
    a3 = sin(Th1 - Th2);
    a4 = (2 - cos(Th1 - Th2) * cos(Th1 - Th2));
    Th1_dd = (a1 + (a2 * a3)) / a4;

    b1 = 2 * GRAVITY * (sin(Th1) * cos(Th1 - Th2) - sin(Th2));
    b2 = 2 * Th1_d * Th1_d + Th2_d * Th2_d * cos(Th1 - Th2);
    b3 = sin(Th1 - Th2);
    b4 = (2 - cos(Th1 - Th2) * cos(Th1 - Th2));

    Th2_dd = (b1 + (b2 * b3)) / b4;

    gsl_vector_set(out_state, 0, Th1_d);
    gsl_vector_set(out_state, 1, Th1_dd);
    gsl_vector_set(out_state, 2, Th2_d);
    gsl_vector_set(out_state, 3, Th2_dd);
    return 0;
}
int test_sine(double t, gsl_vector * in_state, gsl_vector * out_state){
    double x1 = gsl_vector_get(in_state, 0);
    double x2 = gsl_vector_get(in_state, 1);
    
    x1 *= -1;

    gsl_vector_set(out_state, 0, x2);
    gsl_vector_set(out_state, 1, x1);
    return 0;
}

int test_high_pow(double t, gsl_vector * in_state, gsl_vector * out_state){
    double x1 = gsl_vector_get(in_state, 0);
    double x2 = gsl_vector_get(in_state, 1);
    
    x1 *= -1;
    gsl_vector_set(out_state, 0, pow(x2, 9));
    gsl_vector_set(out_state, 1, pow(x1, 9));
    return 0;
}

static double csc(double arg){
    return(1/sin(arg));
}
static double cot(double arg){
    return(1/tan(arg));
}

int polar2d_lagr(double t, gsl_vector * in_state, gsl_vector * out_state){
       
    double r, r_d, Th, Th_d;
    
    r = gsl_vector_get(in_state, 0);
    r_d = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    Th_d = gsl_vector_get(in_state, 3);
    
    double r_dd, Th_dd;
    
    r_dd = r*pow(Th_d, 2);
    Th_dd = (- (2 * r * r_d * Th_d)) / (pow(r, 2)); 

    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, r_dd);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, Th_dd);    

return 0;
}

int polar2d_ham(double t, gsl_vector * in_state, gsl_vector * out_state){
    
    double r, P_r, Th, PTh;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    PTh = gsl_vector_get(in_state, 3);
    
    double r_d, P_r_d, Th_d, PTh_d;
    
    r_d = 2 * P_r;
    P_r_d = 2 * pow(PTh, 2) / (pow(r, 3));
    Th_d = 2 * PTh / pow(r, 2);
    PTh_d = 0;

    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, PTh_d);    

return 0;
}

     ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////							KERR METRICS                          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////

// The equations of motion for the Kerr slow-rotation expansion to second order
// NOTE:  This expects the state vector to be in coordinates and momenta!
// NOT coordinates and time derivatives!!

int kerr_eom_SR_o2(double t, gsl_vector * in_state, gsl_vector * out_state){

    double r, P_r, Th, PTh, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    PTh = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, PTh_d, Ph_d, P_Ph_d;
    
    ///////////////////////////////////////////
    //         coordinate derivatives        //
    /////////////////////////////////////////// 
    r_d = (2*P_r*(2*pow(SPIN,2)*MASS*pow(cos(Th),2) - 2*MASS*pow(r,2) + pow(r,3) + pow(SPIN,2)*r*pow(sin(Th),2)))/pow(r,3);
    Th_d = (2*PTh*(-(pow(SPIN,2)*pow(cos(Th),2)) + pow(r,2)))/pow(r,4);
    Ph_d = (-4*pow(SPIN,2)*JZ*MASS*pow(cot(Th),2) + 2*pow(SPIN,2)*JZ*pow(csc(Th),2)*r + (-4*SPIN*ENERGY*MASS + 4*JZ*MASS*pow(csc(Th),2))*pow(r,2) - 2*JZ*pow(csc(Th),2)*pow(r,3))/((2*MASS - r)*pow(r,4));
    
    //////////////////////////////////////////////
    //          momentum derivatives            //
    //////////////////////////////////////////////
    P_r_d = -((32*pow(SPIN,2)*pow(JZ,2)*pow(MASS,3)*pow(cot(Th),2) - 12*pow(SPIN,2)*pow(MASS,2)*(pow(JZ,2)*(3*pow(cot(Th),2) + pow(csc(Th),2)) + 4*pow(MASS,2)*pow(cos(Th),2)*pow(P_r,2))*r + (2*(pow(JZ,2) + pow(ENERGY,2)*pow(MASS,2) - pow(ENERGY,2)*pow(MASS,2)*cos(2*Th))*pow(csc(Th),2) + (pow(SPIN,2) + 12*pow(MASS,2) - pow(SPIN,2)*cos(2*Th))*pow(P_r,2))*pow(r,5) - 2*MASS*(pow(ENERGY,2) + pow(P_r,2))*pow(r,6) - 2*pow(PTh,2)*pow(2*MASS - r,3)*(-2*pow(SPIN,2)*pow(cos(Th),2) + pow(r,2)) + 3*MASS*pow(r,4)*(3*pow(SPIN,2)*pow(ENERGY,2) + 4*SPIN*ENERGY*JZ - pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2) - 4*pow(JZ,2)*pow(csc(Th),2) + (-pow(SPIN,2) - 8*pow(MASS,2) + 3*pow(SPIN,2)*cos(2*Th))*pow(P_r,2) - 3*pow(SPIN,2)*pow(ENERGY,2)*pow(sin(Th),2)) + 2*MASS*pow(r,2)*(16*SPIN*ENERGY*JZ*pow(MASS,2) + 5*pow(SPIN,2)*pow(JZ,2)*pow(cot(Th),2) + 7*pow(SPIN,2)*pow(JZ,2)*pow(csc(Th),2) - 8*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),2) + 2*pow(SPIN,2)*pow(MASS,2)*(7 + 11*cos(2*Th))*pow(P_r,2) - 8*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(sin(Th),2)) + pow(r,3)*(-6*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2) - 40*SPIN*ENERGY*JZ*pow(MASS,2) + 2*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(cos(Th),2) - 4*pow(SPIN,2)*pow(JZ,2)*pow(csc(Th),2) + 24*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),2) + 2*pow(MASS,2)*(-3*pow(SPIN,2) + 8*pow(MASS,2) - 15*pow(SPIN,2)*cos(2*Th))*pow(P_r,2) + 22*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(sin(Th),2)))/(pow(2*MASS - r,3)*pow(r,5)));           
    PTh_d = (2*(2*pow(SPIN,2)*pow(JZ,2)*MASS*cot(Th)*pow(csc(Th),3) + pow(SPIN,2)*cos(Th)*pow(PTh,2)*(2*MASS - r) - pow(SPIN,2)*(pow(JZ,2)*cot(Th)*pow(csc(Th),3) + 4*pow(MASS,2)*cos(Th)*pow(P_r,2))*r + 2*MASS*(pow(SPIN,2)*pow(ENERGY,2)*cos(Th) - pow(JZ,2)*cot(Th)*pow(csc(Th),3) + 2*pow(SPIN,2)*cos(Th)*pow(P_r,2))*pow(r,2) + (pow(JZ,2)*cot(Th)*pow(csc(Th),3) - pow(SPIN,2)*cos(Th)*pow(P_r,2))*pow(r,3))*sin(Th))/(pow(r,4)*(-2*MASS + r));
    P_Ph_d = 0;
    
    gsl_vector_set(out_state, 0, r_d/2);
    gsl_vector_set(out_state, 1, P_r_d/2);
    gsl_vector_set(out_state, 2, Th_d/2);
    gsl_vector_set(out_state, 3, PTh_d/2);
    gsl_vector_set(out_state, 4, Ph_d/2);
    gsl_vector_set(out_state, 5, P_Ph_d/2); 

return 0;
}


int kerr_eom_exact(double t, gsl_vector * in_state, gsl_vector * out_state){

    double r, P_r, Th, PTh, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    PTh = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, PTh_d, Ph_d, P_Ph_d;
    
    r_d = (2*P_r*(pow(SPIN,2) - 2*MASS*r + pow(r,2)))/(pow(SPIN,2)*
             pow(cos(Th),2) + pow(r,2));
    Th_d = (2*PTh)/(pow(SPIN,2)*pow(cos(Th),2) + pow(r,2));
    Ph_d = (2*(pow(SPIN,2)*JZ*pow(cot(Th),2) + 2*MASS*(SPIN*ENERGY - JZ*
               pow(csc(Th),2))*r + JZ*pow(csc(Th),2)*pow(r,2)))/
               ((pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*(pow(SPIN,2) - 2*MASS*
               r + pow(r,2)));
    
    
    P_r_d = -((pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*(pow(SPIN,2) - 2*MASS*
               r + pow(r,2))*(2*(-(pow(SPIN,2)*pow(ENERGY,2)*(3 + cos(2*
               Th))) + 2*pow(JZ,2)*pow(csc(Th),2) + 2*pow(PTh,2))*
               r - 8*pow(ENERGY,2)*pow(r,3) + 4*pow(P_r,2)*(-2*MASS + 2*
               r)*(pow(SPIN,2) - 2*MASS*r + pow(r,2)) + 4*MASS*
               (-pow(PTh,2) - pow(JZ*csc(Th) - SPIN*ENERGY*sin(Th),2))) -
               (-2*MASS + 2*r)*(pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*
               (2*pow(SPIN,2)*(-(pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2)) +
               pow(JZ,2)*pow(cot(Th),2) + pow(PTh,2)) + (-(pow(SPIN,2)*
               pow(ENERGY,2)*(3 + cos(2*Th))) + 2*pow(JZ,2)*pow(csc(Th),2) +
               2*pow(PTh,2))*pow(r,2) - 2*pow(ENERGY,2)*pow(r,4) +
               2*pow(P_r,2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2) +
               4*MASS*r*(-pow(PTh,2) - pow(JZ*csc(Th) - SPIN*ENERGY*
               sin(Th),2))) - 2*r*(pow(SPIN,2) - 2*MASS*r + pow(r,2))*
               (2*pow(SPIN,2)*(-(pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2)) +
               pow(JZ,2)*pow(cot(Th),2) + pow(PTh,2)) + (-(pow(SPIN,2)*
               pow(ENERGY,2)*(3 + cos(2*Th))) + 2*pow(JZ,2)*pow(csc(Th),2) +
               2*pow(PTh,2))*pow(r,2) - 2*pow(ENERGY,2)*pow(r,4) + 2*
               pow(P_r,2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2) + 4*MASS*
               r*(-pow(PTh,2) - pow(JZ*csc(Th) - SPIN*ENERGY*
               sin(Th),2))))/(2.*pow(pow(SPIN,2)*pow(cos(Th),2) +
               pow(r,2),2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2));
    PTh_d = -((2*MASS*pow(r,3)*(2*pow(JZ,2)*cot(Th)*pow(csc(Th),2) -
                pow(SPIN,2)*pow(ENERGY,2)*sin(2*Th) - 2*pow(SPIN,2)*pow(P_r,2)*
                sin(2*Th)) + pow(r,4)*(-2*pow(JZ,2)*cot(Th)*
                pow(csc(Th),2) + pow(SPIN,2)*pow(P_r,2)*sin(2*Th)) +
                pow(SPIN,4)*(-2*pow(JZ,2)*pow(cos(Th),2)*pow(cot(Th),3) +
                pow(SPIN,2)*pow(P_r,2)*sin(2*Th) + pow(PTh,2)*
                sin(2*Th)) - 2*pow(SPIN,2)*MASS*r*(2*pow(JZ,2)*cot(Th) -
                2*pow(JZ,2)*pow(cot(Th),3) + 2*pow(SPIN,2)*pow(ENERGY,2)*
                pow(cos(Th),3)*sin(Th) + 2*pow(SPIN,2)*pow(ENERGY,2)*
                cos(Th)*pow(sin(Th),3) - 2*SPIN*ENERGY*JZ*sin(2*Th) +
                2*pow(SPIN,2)*pow(P_r,2)*sin(2*Th) + pow(PTh,2)*
                sin(2*Th)) + pow(SPIN,2)*pow(r,2)*(-4*pow(JZ,2)*
                pow(cot(Th),3) + 2*(pow(SPIN,2) + 2*pow(MASS,2))*pow(P_r,2)*
                sin(2*Th) + pow(PTh,2)*sin(2*Th)))/(pow(pow(SPIN,2)*
                pow(cos(Th),2) + pow(r,2),2)*(pow(SPIN,2) - 2*MASS*r +
                pow(r,2))));
    
    P_Ph_d = 0;
    
    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, PTh_d);
    gsl_vector_set(out_state, 4, Ph_d);
    gsl_vector_set(out_state, 5, P_Ph_d); 
       
return 0;
}


int edgb_a2_z1(double t, gsl_vector * in_state, gsl_vector * out_state){

    double r, P_r, Th, P_Th, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    P_Th = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, P_Th_d, Ph_d, P_Ph_d;

    r_d = P_r + ZETA*((-368*pow(MASS,7)*P_r)/(3.*pow(r,7)) + 
      (16*pow(MASS,6)*P_r)/(5.*pow(r,6)) + (2*pow(MASS,5)*P_r)/pow(r,5) + 
      (52*pow(MASS,4)*P_r)/(3.*pow(r,4)) + (pow(MASS,3)*P_r)/pow(r,3) + 
      (pow(MASS,2)*P_r)/pow(r,2)) - (2*MASS*P_r)/r + 
   pow(SPIN,2)*(-((P_r*(-2*MASS*pow(cos(Th),2) - r + pow(cos(Th),2)*r))/
         pow(r,3)) + (ZETA*P_r*(-149940000*pow(MASS,9) + 
           449820000*pow(MASS,9)*pow(cos(Th),2) + 216678000*pow(MASS,8)*r - 
           605934000*pow(MASS,8)*pow(cos(Th),2)*r - 
           47154100*pow(MASS,7)*pow(r,2) + 
           330092700*pow(MASS,7)*pow(cos(Th),2)*pow(r,2) + 
           17575950*pow(MASS,6)*pow(r,3) - 
           57005550*pow(MASS,6)*pow(cos(Th),2)*pow(r,3) - 
           14112390*pow(MASS,5)*pow(r,4) + 
           35060670*pow(MASS,5)*pow(cos(Th),2)*pow(r,4) + 
           7361345*pow(MASS,4)*pow(r,5) - 
           26457285*pow(MASS,4)*pow(cos(Th),2)*pow(r,5) + 
           55626*pow(MASS,3)*pow(r,6) - 
           387378*pow(MASS,3)*pow(cos(Th),2)*pow(r,6) - 
           591696*pow(MASS,2)*pow(r,7) + 
           231588*pow(MASS,2)*pow(cos(Th),2)*pow(r,7) - 242571*MASS*pow(r,8) + 
           562338*MASS*pow(cos(Th),2)*pow(r,8) - 55125*pow(r,9)))/
       (110250.*pow(r,11)));

    P_r_d = (ZETA*(-103040*pow(MASS,10)*pow(P_r,2) + 156864*pow(MASS,9)*pow(P_r,2)*r - 
        4000*pow(ENERGY,2)*pow(MASS,8)*pow(r,2) - 
        79536*pow(MASS,8)*pow(P_r,2)*pow(r,2) + 
        3568*pow(ENERGY,2)*pow(MASS,7)*pow(r,3) + 
        21128*pow(MASS,7)*pow(P_r,2)*pow(r,3) - 
        180*pow(ENERGY,2)*pow(MASS,6)*pow(r,4) - 
        11508*pow(MASS,6)*pow(P_r,2)*pow(r,4) + 
        190*pow(ENERGY,2)*pow(MASS,5)*pow(r,5) + 
        5790*pow(MASS,5)*pow(P_r,2)*pow(r,5) - 
        510*pow(ENERGY,2)*pow(MASS,4)*pow(r,6) - 
        1130*pow(MASS,4)*pow(P_r,2)*pow(r,6) - 
        15*pow(ENERGY,2)*pow(MASS,3)*pow(r,7) + 
        135*pow(MASS,3)*pow(P_r,2)*pow(r,7) - 
        30*pow(MASS,2)*pow(P_r,2)*pow(r,8)))/
    (30.*pow(2*MASS - r,3)*pow(r,8)) + 
   ((2*pow(JZ,2)*pow(csc(Th),2))/pow(r,3) + 
      (2*pow(P_Th,2))/pow(r,3) - (2*MASS*pow(P_r,2))/pow(r,2) - 
      (pow(ENERGY,2)*r)/pow(-2*MASS + r,2) + 
      (pow(ENERGY,2))/(-2*MASS + r))/2. + 
   SPIN*((-2*ENERGY*JZ*(4*pow(MASS,2) - 3*MASS*r))/
       (pow(2*MASS - r,2)*pow(r,3)) + 
      (ENERGY*JZ*ZETA*(6944*pow(MASS,8) - 5616*pow(MASS,7)*r + 
           68*pow(MASS,6)*pow(r,2) - 566*pow(MASS,5)*pow(r,3) + 
           738*pow(MASS,4)*pow(r,4) + 45*pow(MASS,3)*pow(r,5)))/
       (15.*pow(2*MASS - r,3)*pow(r,8))) + 
   pow(SPIN,2)*(((-16*pow(JZ,2)*pow(MASS,3)*pow(csc(Th),2) + 
           16*pow(JZ,2)*pow(MASS,3)*pow(csc(Th),4) + 
           16*pow(MASS,3)*pow(cot(Th),2)*pow(P_Th,2) + 
           18*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),2)*r - 
           24*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),4)*r - 
           24*pow(MASS,4)*pow(cot(Th),2)*pow(P_r,2)*r - 
           24*pow(MASS,2)*pow(cot(Th),2)*pow(P_Th,2)*r - 
           8*pow(ENERGY,2)*pow(MASS,3)*pow(r,2) - 
           5*pow(JZ,2)*MASS*pow(csc(Th),2)*pow(r,2) + 
           12*pow(JZ,2)*MASS*pow(csc(Th),4)*pow(r,2) + 
           44*pow(MASS,3)*pow(cot(Th),2)*pow(P_r,2)*pow(r,2) - 
           8*pow(MASS,3)*pow(csc(Th),2)*pow(P_r,2)*pow(r,2) + 
           12*MASS*pow(cot(Th),2)*pow(P_Th,2)*pow(r,2) + 
           8*pow(ENERGY,2)*pow(MASS,2)*pow(r,3) - 
           2*pow(ENERGY,2)*pow(MASS,2)*pow(cot(Th),2)*pow(r,3) - 
           2*pow(JZ,2)*pow(csc(Th),4)*pow(r,3) - 
           30*pow(MASS,2)*pow(cot(Th),2)*pow(P_r,2)*pow(r,3) + 
           12*pow(MASS,2)*pow(csc(Th),2)*pow(P_r,2)*pow(r,3) - 
           2*pow(cot(Th),2)*pow(P_Th,2)*pow(r,3) + 
           3*pow(ENERGY,2)*MASS*pow(cot(Th),2)*pow(r,4) + 
           9*MASS*pow(cot(Th),2)*pow(P_r,2)*pow(r,4) - 
           6*MASS*pow(csc(Th),2)*pow(P_r,2)*pow(r,4) - 
           pow(cot(Th),2)*pow(P_r,2)*pow(r,5) + 
           pow(csc(Th),2)*pow(P_r,2)*pow(r,5))*pow(sin(Th),2))/
       (pow(r,5)*pow(-2*MASS + r,3)) + 
      (ZETA*(5080320000*pow(JZ,2)*pow(MASS,12)*pow(cot(Th),2)*
            pow(csc(Th),2) - 1693440000*pow(JZ,2)*pow(MASS,12)*
            pow(csc(Th),4) + 5080320000*pow(MASS,12)*pow(cot(Th),2)*
            pow(P_Th,2) - 1693440000*pow(MASS,12)*pow(csc(Th),2)*pow(P_Th,2) + 
           517440000*pow(JZ,2)*pow(MASS,11)*pow(csc(Th),2)*r - 
           13441209600*pow(JZ,2)*pow(MASS,11)*pow(cot(Th),2)*
            pow(csc(Th),2)*r + 
           4480403200*pow(JZ,2)*pow(MASS,11)*pow(csc(Th),4)*r + 
           79168320000*pow(MASS,13)*pow(cot(Th),2)*pow(P_r,2)*r - 
           26389440000*pow(MASS,13)*pow(csc(Th),2)*pow(P_r,2)*r - 
           13441209600*pow(MASS,11)*pow(cot(Th),2)*pow(P_Th,2)*r + 
           4480403200*pow(MASS,11)*pow(csc(Th),2)*pow(P_Th,2)*r - 
           1260672000*pow(JZ,2)*pow(MASS,10)*pow(csc(Th),2)*
            pow(r,2) + 12541603200*pow(JZ,2)*pow(MASS,10)*
            pow(cot(Th),2)*pow(csc(Th),2)*pow(r,2) - 
           4180534400*pow(JZ,2)*pow(MASS,10)*pow(csc(Th),4)*
            pow(r,2) - 255286080000*pow(MASS,12)*pow(cot(Th),2)*pow(P_r,2)*
            pow(r,2) + 87447360000*pow(MASS,12)*pow(csc(Th),2)*pow(P_r,2)*
            pow(r,2) + 12541603200*pow(MASS,10)*pow(cot(Th),2)*pow(P_Th,2)*
            pow(r,2) - 4180534400*pow(MASS,10)*pow(csc(Th),2)*pow(P_Th,2)*
            pow(r,2) - 952560000*pow(ENERGY,2)*pow(MASS,11)*
            pow(cot(Th),2)*pow(r,3) + 
           950443200*pow(JZ,2)*pow(MASS,9)*pow(csc(Th),2)*pow(r,3) + 
           317520000*pow(ENERGY,2)*pow(MASS,11)*pow(csc(Th),2)*
            pow(r,3) - 4982623200*pow(JZ,2)*pow(MASS,9)*
            pow(cot(Th),2)*pow(csc(Th),2)*pow(r,3) + 
           1660874400*pow(JZ,2)*pow(MASS,9)*pow(csc(Th),4)*
            pow(r,3) + 360184708800*pow(MASS,11)*pow(cot(Th),2)*pow(P_r,2)*
            pow(r,3) - 115711310400*pow(MASS,11)*pow(csc(Th),2)*pow(P_r,2)*
            pow(r,3) - 4982623200*pow(MASS,9)*pow(cot(Th),2)*pow(P_Th,2)*
            pow(r,3) + 1660874400*pow(MASS,9)*pow(csc(Th),2)*pow(P_Th,2)*
            pow(r,3) + 1843027200*pow(ENERGY,2)*pow(MASS,10)*
            pow(cot(Th),2)*pow(r,4) - 
           300585600*pow(JZ,2)*pow(MASS,8)*pow(csc(Th),2)*pow(r,4) - 
           645702400*pow(ENERGY,2)*pow(MASS,10)*pow(csc(Th),2)*
            pow(r,4) + 1441712640*pow(JZ,2)*pow(MASS,8)*
            pow(cot(Th),2)*pow(csc(Th),2)*pow(r,4) - 
           480570880*pow(JZ,2)*pow(MASS,8)*pow(csc(Th),4)*pow(r,4) - 
           287371728000*pow(MASS,10)*pow(cot(Th),2)*pow(P_r,2)*pow(r,4) + 
           81027542400*pow(MASS,10)*pow(csc(Th),2)*pow(P_r,2)*pow(r,4) + 
           1441712640*pow(MASS,8)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,4) - 
           480570880*pow(MASS,8)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,4) - 
           204153600*pow(ENERGY,2)*pow(MASS,9)*pow(r,5) - 
           10878000*pow(ENERGY,2)*pow(MASS,9)*pow(cot(Th),2)*pow(r,5) + 
           177752400*pow(JZ,2)*pow(MASS,7)*pow(csc(Th),2)*pow(r,5) + 
           370381200*pow(ENERGY,2)*pow(MASS,9)*pow(csc(Th),2)*
            pow(r,5) - 999742992*pow(JZ,2)*pow(MASS,7)*
            pow(cot(Th),2)*pow(csc(Th),2)*pow(r,5) + 
           333247664*pow(JZ,2)*pow(MASS,7)*pow(csc(Th),4)*pow(r,5) + 
           143242979040*pow(MASS,9)*pow(cot(Th),2)*pow(P_r,2)*pow(r,5) - 
           35248896480*pow(MASS,9)*pow(csc(Th),2)*pow(P_r,2)*pow(r,5) - 
           999742992*pow(MASS,7)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,5) + 
           333247664*pow(MASS,7)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,5) + 
           179692800*pow(ENERGY,2)*pow(MASS,8)*pow(r,6) - 
           1422834000*pow(ENERGY,2)*pow(MASS,8)*pow(cot(Th),2)*
            pow(r,6) - 130771200*pow(JZ,2)*pow(MASS,6)*
            pow(csc(Th),2)*pow(r,6) + 
           9562000*pow(ENERGY,2)*pow(MASS,8)*pow(csc(Th),2)*pow(r,6) + 
           396637704*pow(JZ,2)*pow(MASS,6)*pow(cot(Th),2)*
            pow(csc(Th),2)*pow(r,6) - 
           132212568*pow(JZ,2)*pow(MASS,6)*pow(csc(Th),4)*pow(r,6) - 
           51164569440*pow(MASS,8)*pow(cot(Th),2)*pow(P_r,2)*pow(r,6) + 
           12804322080*pow(MASS,8)*pow(csc(Th),2)*pow(P_r,2)*pow(r,6) + 
           396637704*pow(MASS,6)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,6) - 
           132212568*pow(MASS,6)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,6) - 
           4821600*pow(ENERGY,2)*pow(MASS,7)*pow(r,7) + 
           529716900*pow(ENERGY,2)*pow(MASS,7)*pow(cot(Th),2)*
            pow(r,7) + 26195400*pow(JZ,2)*pow(MASS,5)*pow(csc(Th),2)*
            pow(r,7) - 4876300*pow(ENERGY,2)*pow(MASS,7)*pow(csc(Th),2)*
            pow(r,7) - 22742238*pow(JZ,2)*pow(MASS,5)*pow(cot(Th),2)*
            pow(csc(Th),2)*pow(r,7) + 
           7580746*pow(JZ,2)*pow(MASS,5)*pow(csc(Th),4)*pow(r,7) + 
           17558190540*pow(MASS,7)*pow(cot(Th),2)*pow(P_r,2)*pow(r,7) - 
           5329057380*pow(MASS,7)*pow(csc(Th),2)*pow(P_r,2)*pow(r,7) - 
           22742238*pow(MASS,5)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,7) + 
           7580746*pow(MASS,5)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,7) + 
           14582400*pow(ENERGY,2)*pow(MASS,6)*pow(r,8) + 
           97514160*pow(ENERGY,2)*pow(MASS,6)*pow(cot(Th),2)*pow(r,8) + 
           2116800*pow(JZ,2)*pow(MASS,4)*pow(csc(Th),2)*pow(r,8) - 
           42579120*pow(ENERGY,2)*pow(MASS,6)*pow(csc(Th),2)*pow(r,8) - 
           15053304*pow(JZ,2)*pow(MASS,4)*pow(cot(Th),2)*
            pow(csc(Th),2)*pow(r,8) + 
           5017768*pow(JZ,2)*pow(MASS,4)*pow(csc(Th),4)*pow(r,8) - 
           6152488848*pow(MASS,6)*pow(cot(Th),2)*pow(P_r,2)*pow(r,8) + 
           1944166416*pow(MASS,6)*pow(csc(Th),2)*pow(P_r,2)*pow(r,8) - 
           15053304*pow(MASS,4)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,8) + 
           5017768*pow(MASS,4)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,8) - 
           25578000*pow(ENERGY,2)*pow(MASS,5)*pow(r,9) + 
           31276254*pow(ENERGY,2)*pow(MASS,5)*pow(cot(Th),2)*pow(r,9) - 
           346118*pow(ENERGY,2)*pow(MASS,5)*pow(csc(Th),2)*pow(r,9) + 
           21578193*pow(JZ,2)*pow(MASS,3)*pow(cot(Th),2)*
            pow(csc(Th),2)*pow(r,9) - 
           7192731*pow(JZ,2)*pow(MASS,3)*pow(csc(Th),4)*pow(r,9) + 
           1466237970*pow(MASS,5)*pow(cot(Th),2)*pow(P_r,2)*pow(r,9) - 
           381362490*pow(MASS,5)*pow(csc(Th),2)*pow(P_r,2)*pow(r,9) + 
           21578193*pow(MASS,3)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,9) - 
           7192731*pow(MASS,3)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,9) - 
           1587600*pow(ENERGY,2)*pow(MASS,4)*pow(r,10) - 
           45249858*pow(ENERGY,2)*pow(MASS,4)*pow(cot(Th),2)*
            pow(r,10) + 6865986*pow(ENERGY,2)*pow(MASS,4)*
            pow(csc(Th),2)*pow(r,10) - 
           14653800*pow(JZ,2)*pow(MASS,2)*pow(cot(Th),2)*
            pow(csc(Th),2)*pow(r,10) + 
           4884600*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),4)*pow(r,10) - 
           175000590*pow(MASS,4)*pow(cot(Th),2)*pow(P_r,2)*pow(r,10) + 
           6663030*pow(MASS,4)*pow(csc(Th),2)*pow(P_r,2)*pow(r,10) - 
           14653800*pow(MASS,2)*pow(cot(Th),2)*pow(P_Th,2)*pow(r,10) + 
           4884600*pow(MASS,2)*pow(csc(Th),2)*pow(P_Th,2)*pow(r,10) - 
           2025594*pow(ENERGY,2)*pow(MASS,3)*pow(cot(Th),2)*pow(r,11) + 
           2733198*pow(ENERGY,2)*pow(MASS,3)*pow(csc(Th),2)*pow(r,11) + 
           2811690*pow(JZ,2)*MASS*pow(cot(Th),2)*pow(csc(Th),2)*
            pow(r,11) - 937230*pow(JZ,2)*MASS*pow(csc(Th),4)*
            pow(r,11) + 31140630*pow(MASS,3)*pow(cot(Th),2)*pow(P_r,2)*
            pow(r,11) + 5275290*pow(MASS,3)*pow(csc(Th),2)*pow(P_r,2)*
            pow(r,11) + 2811690*MASS*pow(cot(Th),2)*pow(P_Th,2)*
            pow(r,11) - 937230*MASS*pow(csc(Th),2)*pow(P_Th,2)*pow(r,11) - 
           2249352*pow(ENERGY,2)*pow(MASS,2)*pow(cot(Th),2)*pow(r,12) - 
           132216*pow(ENERGY,2)*pow(MASS,2)*pow(csc(Th),2)*pow(r,12) - 
           12569760*pow(MASS,2)*pow(cot(Th),2)*pow(P_r,2)*pow(r,12) + 
           808920*pow(MASS,2)*pow(csc(Th),2)*pow(P_r,2)*pow(r,12) + 
           1687014*pow(ENERGY,2)*MASS*pow(cot(Th),2)*pow(r,13) - 
           617463*pow(ENERGY,2)*MASS*pow(csc(Th),2)*pow(r,13) + 
           1687014*MASS*pow(cot(Th),2)*pow(P_r,2)*pow(r,13) + 
           154287*MASS*pow(csc(Th),2)*pow(P_r,2)*pow(r,13) - 
           110250*pow(csc(Th),2)*pow(P_r,2)*pow(r,14))*pow(sin(Th),2))/
       (220500.*pow(r,13)*pow(-2*MASS + r,4)));

    Th_d = P_Th/pow(r,2) + pow(SPIN,2)*(-((pow(cos(Th),2)*P_Th)/pow(r,4)) + 
      (MASS*ZETA*(-1 + 3*pow(cos(Th),2))*P_Th*
         (8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
           3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
           887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
           435540*MASS*pow(r,6) + 187446*pow(r,7)))/(110250.*pow(r,12)));

    P_Th_d = (pow(JZ,2)*cot(Th)*pow(csc(Th),2))/pow(r,2) + 
   pow(SPIN,2)*(((4*pow(ENERGY,2)*MASS*cot(Th))/(r*pow(-2*MASS + r,2)) - 
         (pow(JZ,2)*(4*MASS*cot(Th)*pow(csc(Th),2) - 
              2*cot(Th)*pow(csc(Th),2)*r))/((2*MASS - r)*pow(r,4)) - 
         (2*cos(Th)*pow(P_Th,2)*sin(Th))/pow(r,4) - 
         (4*pow(ENERGY,2)*cos(Th)*(2*pow(MASS,2) + MASS*pow(cot(Th),2)*r)*
            sin(Th))/(pow(r,2)*pow(-2*MASS + r,2)) + 
         (pow(P_r,2)*(4*MASS*cos(Th)*sin(Th) - 2*cos(Th)*r*sin(Th)))/
          pow(r,3))/2. + (ZETA*(-(pow(JZ,2)*
               (-141120000*pow(MASS,10)*cot(Th)*pow(csc(Th),2) + 
                 240531200*pow(MASS,9)*cot(Th)*pow(csc(Th),2)*r - 
                 80024000*pow(MASS,8)*cot(Th)*pow(csc(Th),2)*pow(r,2) - 
                 124000*pow(MASS,7)*cot(Th)*pow(csc(Th),2)*pow(r,3) - 
                 30217360*pow(MASS,6)*cot(Th)*pow(csc(Th),2)*pow(r,4) + 
                 8804632*pow(MASS,5)*cot(Th)*pow(csc(Th),2)*pow(r,5) + 
                 2294648*pow(MASS,4)*cot(Th)*pow(csc(Th),2)*pow(r,6) + 
                 766572*pow(MASS,3)*cot(Th)*pow(csc(Th),2)*pow(r,7) + 
                 1256976*pow(MASS,2)*cot(Th)*pow(csc(Th),2)*pow(r,8) - 
                 749784*MASS*cot(Th)*pow(csc(Th),2)*pow(r,9)))/
            (110250.*pow(r,12)*pow(-2*MASS + r,2)) + 
           (MASS*cos(Th)*pow(P_Th,2)*
              (8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 
                3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
                887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 
                435540*MASS*pow(r,6) + 187446*pow(r,7))*sin(Th))/
            (18375.*pow(r,12)) - 
           (pow(ENERGY,2)*cos(Th)*
              (52920000*pow(MASS,10)*pow(cot(Th),2) - 
                17640000*pow(MASS,10)*pow(csc(Th),2) - 
                75499200*pow(MASS,9)*pow(cot(Th),2)*r + 
                27126400*pow(MASS,9)*pow(csc(Th),2)*r + 
                14582400*pow(MASS,8)*pow(r,2) - 
                58543800*pow(MASS,8)*pow(cot(Th),2)*pow(r,2) - 
                5142200*pow(MASS,8)*pow(csc(Th),2)*pow(r,2) - 
                2822400*pow(MASS,7)*pow(r,3) + 
                69783000*pow(MASS,7)*pow(cot(Th),2)*pow(r,3) - 
                5082000*pow(MASS,7)*pow(csc(Th),2)*pow(r,3) - 
                2058000*pow(MASS,6)*pow(r,4) + 
                9833010*pow(MASS,6)*pow(cot(Th),2)*pow(r,4) - 
                4086170*pow(MASS,6)*pow(csc(Th),2)*pow(r,4) - 
                3880800*pow(MASS,5)*pow(r,5) - 
                2356260*pow(MASS,5)*pow(cot(Th),2)*pow(r,5) + 
                1236220*pow(MASS,5)*pow(csc(Th),2)*pow(r,5) - 
                264600*pow(MASS,4)*pow(r,6) - 
                7961679*pow(MASS,4)*pow(cot(Th),2)*pow(r,6) + 
                1499943*pow(MASS,4)*pow(csc(Th),2)*pow(r,6) - 
                630054*pow(MASS,3)*pow(cot(Th),2)*pow(r,7) + 
                533418*pow(MASS,3)*pow(csc(Th),2)*pow(r,7) - 
                562338*pow(MASS,2)*pow(cot(Th),2)*pow(r,8) - 
                33054*pow(MASS,2)*pow(csc(Th),2)*pow(r,8) + 
                562338*MASS*pow(cot(Th),2)*pow(r,9) - 
                205821*MASS*pow(csc(Th),2)*pow(r,9))*sin(Th))/
            (55125.*pow(r,9)*pow(-2*MASS + r,3)) - 
           (pow(ENERGY,2)*(-70560000*pow(MASS,10)*cot(Th)*pow(csc(Th),2) + 
                96745600*pow(MASS,9)*cot(Th)*pow(csc(Th),2)*r + 
                127372000*pow(MASS,8)*cot(Th)*pow(csc(Th),2)*pow(r,2) - 
                129402000*pow(MASS,7)*cot(Th)*pow(csc(Th),2)*pow(r,3) - 
                11493680*pow(MASS,6)*cot(Th)*pow(csc(Th),2)*pow(r,4) + 
                2240080*pow(MASS,5)*cot(Th)*pow(csc(Th),2)*pow(r,5) + 
                12923472*pow(MASS,4)*cot(Th)*pow(csc(Th),2)*pow(r,6) + 
                193272*pow(MASS,3)*cot(Th)*pow(csc(Th),2)*pow(r,7) + 
                1190784*pow(MASS,2)*cot(Th)*pow(csc(Th),2)*pow(r,8) - 
                713034*MASS*cot(Th)*pow(csc(Th),2)*pow(r,9))*pow(sin(Th),2))/
            (110250.*pow(r,9)*pow(-2*MASS + r,3)) - 
           (pow(P_r,2)*(-899640000*pow(MASS,9)*cos(Th)*sin(Th) + 
                1211868000*pow(MASS,8)*cos(Th)*r*sin(Th) - 
                660185400*pow(MASS,7)*cos(Th)*pow(r,2)*sin(Th) + 
                114011100*pow(MASS,6)*cos(Th)*pow(r,3)*sin(Th) - 
                70121340*pow(MASS,5)*cos(Th)*pow(r,4)*sin(Th) + 
                52914570*pow(MASS,4)*cos(Th)*pow(r,5)*sin(Th) + 
                774756*pow(MASS,3)*cos(Th)*pow(r,6)*sin(Th) - 
                463176*pow(MASS,2)*cos(Th)*pow(r,7)*sin(Th) - 
                1124676*MASS*cos(Th)*pow(r,8)*sin(Th)))/(110250.*pow(r,11))))/
       2.);

    Ph_d = (JZ*pow(csc(Th),2))/pow(r,2) + 
   SPIN*((-2*ENERGY*MASS)/((2*MASS - r)*pow(r,2)) - 
      (ENERGY*ZETA*(-496*pow(MASS,7) + 96*pow(MASS,6)*r + 
           70*pow(MASS,5)*pow(r,2) + 132*pow(MASS,4)*pow(r,3) + 
           9*pow(MASS,3)*pow(r,4)))/(15.*pow(r,7)*pow(-2*MASS + r,2))) + 
   pow(SPIN,2)*((JZ*(2*MASS - 2*MASS*pow(csc(Th),2) + pow(csc(Th),2)*r))/
       ((2*MASS - r)*pow(r,4)) + 
      (JZ*ZETA*(105840000*pow(MASS,10)*pow(cot(Th),2) - 
           35280000*pow(MASS,10)*pow(csc(Th),2) + 11760000*pow(MASS,9)*r - 
           180398400*pow(MASS,9)*pow(cot(Th),2)*r + 
           60132800*pow(MASS,9)*pow(csc(Th),2)*r - 
           17404800*pow(MASS,8)*pow(r,2) + 
           60018000*pow(MASS,8)*pow(cot(Th),2)*pow(r,2) - 
           20006000*pow(MASS,8)*pow(csc(Th),2)*pow(r,2) + 
           882000*pow(MASS,7)*pow(r,3) + 
           93000*pow(MASS,7)*pow(cot(Th),2)*pow(r,3) - 
           31000*pow(MASS,7)*pow(csc(Th),2)*pow(r,3) - 
           1764000*pow(MASS,6)*pow(r,4) + 
           22663020*pow(MASS,6)*pow(cot(Th),2)*pow(r,4) - 
           7554340*pow(MASS,6)*pow(csc(Th),2)*pow(r,4) + 
           3733800*pow(MASS,5)*pow(r,5) - 
           6603474*pow(MASS,5)*pow(cot(Th),2)*pow(r,5) + 
           2201158*pow(MASS,5)*pow(csc(Th),2)*pow(r,5) + 
           264600*pow(MASS,4)*pow(r,6) - 
           1720986*pow(MASS,4)*pow(cot(Th),2)*pow(r,6) + 
           573662*pow(MASS,4)*pow(csc(Th),2)*pow(r,6) - 
           574929*pow(MASS,3)*pow(cot(Th),2)*pow(r,7) + 
           191643*pow(MASS,3)*pow(csc(Th),2)*pow(r,7) - 
           942732*pow(MASS,2)*pow(cot(Th),2)*pow(r,8) + 
           314244*pow(MASS,2)*pow(csc(Th),2)*pow(r,8) + 
           562338*MASS*pow(cot(Th),2)*pow(r,9) - 
           187446*MASS*pow(csc(Th),2)*pow(r,9)))/
       (110250.*pow(r,12)*pow(-2*MASS + r,2)));

    P_Ph_d = 0;

    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, P_Th_d);
    gsl_vector_set(out_state, 4, Ph_d);
    gsl_vector_set(out_state, 5, P_Ph_d); 
       
return 0;
}
