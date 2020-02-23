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
/////////  o Coordinates with grENERGYk letters are written as the two-letter abbreviation
/////////    with the first letter capitalized.
/////////    Ex: theta becomes "Th", phi becomes "Ph".
/////////
/////////  Momenta and time derivatives
/////////  o Momenta are written by prepending "P_" to the beginning of the coordinate.
/////////    Ex: The phi momentum becomes "P_Ph"
/////////  o Time derivatives are written by appending "_d" to the end of the momentum
/////////    or coordinate
/////////    Ex: The time derivative of the r coordinate becomes "r_d",
/////////        the time derivative of the phi momentum becomes "P_Ph_d"
///////// 
/////////  Coupling parameters
/////////  o Coupling/expansion parameters are written as the name of their full grENERGYk
/////////    letter, with appropriate capitalization
/////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



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
    
    double r, P_r, Th, P_Th;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    P_Th = gsl_vector_get(in_state, 3);
    
    double r_d, P_r_d, Th_d, P_Th_d;
    
    r_d = 2 * P_r;
    P_r_d = 2 * pow(P_Th, 2) / (pow(r, 3));
    Th_d = 2 * P_Th / pow(r, 2);
    P_Th_d = 0;

    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, P_Th_d);    

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

    double r, P_r, Th, P_Th, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    P_Th = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, P_Th_d, Ph_d, P_Ph_d;
    
    ///////////////////////////////////////////
    //         coordinate derivatives        //
    /////////////////////////////////////////// 
    r_d = (2*P_r*(2*pow(SPIN,2)*MASS*pow(cos(Th),2) - 2*MASS*pow(r,2) + pow(r,3) + pow(SPIN,2)*r*pow(sin(Th),2)))/pow(r,3);
    Th_d = (2*P_Th*(-(pow(SPIN,2)*pow(cos(Th),2)) + pow(r,2)))/pow(r,4);
    Ph_d = (-4*pow(SPIN,2)*JZ*MASS*pow(cot(Th),2) + 2*pow(SPIN,2)*JZ*pow(csc(Th),2)*r + (-4*SPIN*ENERGY*MASS + 4*JZ*MASS*pow(csc(Th),2))*pow(r,2) - 2*JZ*pow(csc(Th),2)*pow(r,3))/((2*MASS - r)*pow(r,4));
    
    //////////////////////////////////////////////
    //          momentum derivatives            //
    //////////////////////////////////////////////
    P_r_d = -((32*pow(SPIN,2)*pow(JZ,2)*pow(MASS,3)*pow(cot(Th),2) - 12*pow(SPIN,2)*pow(MASS,2)*(pow(JZ,2)*(3*pow(cot(Th),2) + pow(csc(Th),2)) + 4*pow(MASS,2)*pow(cos(Th),2)*pow(P_r,2))*r + (2*(pow(JZ,2) + pow(ENERGY,2)*pow(MASS,2) - pow(ENERGY,2)*pow(MASS,2)*cos(2*Th))*pow(csc(Th),2) + (pow(SPIN,2) + 12*pow(MASS,2) - pow(SPIN,2)*cos(2*Th))*pow(P_r,2))*pow(r,5) - 2*MASS*(pow(ENERGY,2) + pow(P_r,2))*pow(r,6) - 2*pow(P_Th,2)*pow(2*MASS - r,3)*(-2*pow(SPIN,2)*pow(cos(Th),2) + pow(r,2)) + 3*MASS*pow(r,4)*(3*pow(SPIN,2)*pow(ENERGY,2) + 4*SPIN*ENERGY*JZ - pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2) - 4*pow(JZ,2)*pow(csc(Th),2) + (-pow(SPIN,2) - 8*pow(MASS,2) + 3*pow(SPIN,2)*cos(2*Th))*pow(P_r,2) - 3*pow(SPIN,2)*pow(ENERGY,2)*pow(sin(Th),2)) + 2*MASS*pow(r,2)*(16*SPIN*ENERGY*JZ*pow(MASS,2) + 5*pow(SPIN,2)*pow(JZ,2)*pow(cot(Th),2) + 7*pow(SPIN,2)*pow(JZ,2)*pow(csc(Th),2) - 8*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),2) + 2*pow(SPIN,2)*pow(MASS,2)*(7 + 11*cos(2*Th))*pow(P_r,2) - 8*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(sin(Th),2)) + pow(r,3)*(-6*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2) - 40*SPIN*ENERGY*JZ*pow(MASS,2) + 2*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(cos(Th),2) - 4*pow(SPIN,2)*pow(JZ,2)*pow(csc(Th),2) + 24*pow(JZ,2)*pow(MASS,2)*pow(csc(Th),2) + 2*pow(MASS,2)*(-3*pow(SPIN,2) + 8*pow(MASS,2) - 15*pow(SPIN,2)*cos(2*Th))*pow(P_r,2) + 22*pow(SPIN,2)*pow(ENERGY,2)*pow(MASS,2)*pow(sin(Th),2)))/(pow(2*MASS - r,3)*pow(r,5)));           
    P_Th_d = (2*(2*pow(SPIN,2)*pow(JZ,2)*MASS*cot(Th)*pow(csc(Th),3) + pow(SPIN,2)*cos(Th)*pow(P_Th,2)*(2*MASS - r) - pow(SPIN,2)*(pow(JZ,2)*cot(Th)*pow(csc(Th),3) + 4*pow(MASS,2)*cos(Th)*pow(P_r,2))*r + 2*MASS*(pow(SPIN,2)*pow(ENERGY,2)*cos(Th) - pow(JZ,2)*cot(Th)*pow(csc(Th),3) + 2*pow(SPIN,2)*cos(Th)*pow(P_r,2))*pow(r,2) + (pow(JZ,2)*cot(Th)*pow(csc(Th),3) - pow(SPIN,2)*cos(Th)*pow(P_r,2))*pow(r,3))*sin(Th))/(pow(r,4)*(-2*MASS + r));
    P_Ph_d = 0;
    
    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, P_Th_d);
    gsl_vector_set(out_state, 4, Ph_d);
    gsl_vector_set(out_state, 5, P_Ph_d); 

return 0;
}


int kerr_eom_exact(double t, gsl_vector * in_state, gsl_vector * out_state){

    double r, P_r, Th, P_Th, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    P_Th = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, P_Th_d, Ph_d, P_Ph_d;
    
    r_d = (2*P_r*(pow(SPIN,2) - 2*MASS*r + pow(r,2)))/(pow(SPIN,2)*
             pow(cos(Th),2) + pow(r,2));
    Th_d = (2*P_Th)/(pow(SPIN,2)*pow(cos(Th),2) + pow(r,2));
    Ph_d = (2*(pow(SPIN,2)*JZ*pow(cot(Th),2) + 2*MASS*(SPIN*ENERGY - JZ*
               pow(csc(Th),2))*r + JZ*pow(csc(Th),2)*pow(r,2)))/
               ((pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*(pow(SPIN,2) - 2*MASS*
               r + pow(r,2)));
    
    
    P_r_d = -((pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*(pow(SPIN,2) - 2*MASS*
               r + pow(r,2))*(2*(-(pow(SPIN,2)*pow(ENERGY,2)*(3 + cos(2*
               Th))) + 2*pow(JZ,2)*pow(csc(Th),2) + 2*pow(P_Th,2))*
               r - 8*pow(ENERGY,2)*pow(r,3) + 4*pow(P_r,2)*(-2*MASS + 2*
               r)*(pow(SPIN,2) - 2*MASS*r + pow(r,2)) + 4*MASS*
               (-pow(P_Th,2) - pow(JZ*csc(Th) - SPIN*ENERGY*sin(Th),2))) -
               (-2*MASS + 2*r)*(pow(SPIN,2)*pow(cos(Th),2) + pow(r,2))*
               (2*pow(SPIN,2)*(-(pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2)) +
               pow(JZ,2)*pow(cot(Th),2) + pow(P_Th,2)) + (-(pow(SPIN,2)*
               pow(ENERGY,2)*(3 + cos(2*Th))) + 2*pow(JZ,2)*pow(csc(Th),2) +
               2*pow(P_Th,2))*pow(r,2) - 2*pow(ENERGY,2)*pow(r,4) +
               2*pow(P_r,2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2) +
               4*MASS*r*(-pow(P_Th,2) - pow(JZ*csc(Th) - SPIN*ENERGY*
               sin(Th),2))) - 2*r*(pow(SPIN,2) - 2*MASS*r + pow(r,2))*
               (2*pow(SPIN,2)*(-(pow(SPIN,2)*pow(ENERGY,2)*pow(cos(Th),2)) +
               pow(JZ,2)*pow(cot(Th),2) + pow(P_Th,2)) + (-(pow(SPIN,2)*
               pow(ENERGY,2)*(3 + cos(2*Th))) + 2*pow(JZ,2)*pow(csc(Th),2) +
               2*pow(P_Th,2))*pow(r,2) - 2*pow(ENERGY,2)*pow(r,4) + 2*
               pow(P_r,2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2) + 4*MASS*
               r*(-pow(P_Th,2) - pow(JZ*csc(Th) - SPIN*ENERGY*
               sin(Th),2))))/(2.*pow(pow(SPIN,2)*pow(cos(Th),2) +
               pow(r,2),2)*pow(pow(SPIN,2) - 2*MASS*r + pow(r,2),2));
    P_Th_d = -((2*MASS*pow(r,3)*(2*pow(JZ,2)*cot(Th)*pow(csc(Th),2) -
                pow(SPIN,2)*pow(ENERGY,2)*sin(2*Th) - 2*pow(SPIN,2)*pow(P_r,2)*
                sin(2*Th)) + pow(r,4)*(-2*pow(JZ,2)*cot(Th)*
                pow(csc(Th),2) + pow(SPIN,2)*pow(P_r,2)*sin(2*Th)) +
                pow(SPIN,4)*(-2*pow(JZ,2)*pow(cos(Th),2)*pow(cot(Th),3) +
                pow(SPIN,2)*pow(P_r,2)*sin(2*Th) + pow(P_Th,2)*
                sin(2*Th)) - 2*pow(SPIN,2)*MASS*r*(2*pow(JZ,2)*cot(Th) -
                2*pow(JZ,2)*pow(cot(Th),3) + 2*pow(SPIN,2)*pow(ENERGY,2)*
                pow(cos(Th),3)*sin(Th) + 2*pow(SPIN,2)*pow(ENERGY,2)*
                cos(Th)*pow(sin(Th),3) - 2*SPIN*ENERGY*JZ*sin(2*Th) +
                2*pow(SPIN,2)*pow(P_r,2)*sin(2*Th) + pow(P_Th,2)*
                sin(2*Th)) + pow(SPIN,2)*pow(r,2)*(-4*pow(JZ,2)*
                pow(cot(Th),3) + 2*(pow(SPIN,2) + 2*pow(MASS,2))*pow(P_r,2)*
                sin(2*Th) + pow(P_Th,2)*sin(2*Th)))/(pow(pow(SPIN,2)*
                pow(cos(Th),2) + pow(r,2),2)*(pow(SPIN,2) - 2*MASS*r +
                pow(r,2))));
    
    P_Ph_d = 0;
    
    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, P_Th_d);
    gsl_vector_set(out_state, 4, Ph_d);
    gsl_vector_set(out_state, 5, P_Ph_d); 
       
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
	 ///////							EdGB METRICS                          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ///////														          ///////
	 ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////

// nENERGYd another auxiliary file where coupling parameters can be specified
int edgb_o2(double t, gsl_vector * in_state, gsl_vector * out_state){

    double r, P_r, Th, P_Th, Ph, P_Ph;
    
    r = gsl_vector_get(in_state, 0);
    P_r = gsl_vector_get(in_state, 1);
    Th = gsl_vector_get(in_state, 2);
    P_Th = gsl_vector_get(in_state, 3);
    Ph = gsl_vector_get(in_state, 4);
    P_Ph = gsl_vector_get(in_state, 5);
    
    double r_d, P_r_d, Th_d, P_Th_d, Ph_d, P_Ph_d;
    
    r_d = (2*P_r)/((pow(ALPHA3,2)*(1840*pow(MASS,5) - 48*pow(MASS,4)*r -
    	  30*pow(MASS,3)*pow(r,2) - 260*pow(MASS,2)*pow(r,3) - 15*MASS*pow(r,4)
    	  - 15*pow(r,5)))/(15.*BETA*KAPPA*pow(MASS,2)*pow(r,5)*pow(-2*MASS + r,2))
    	  + (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))/(-2*MASS*r + pow(r,2) +
    	  pow(SPIN,2)) - (pow(ALPHA3,2)*pow(SPIN,2)*(3675*r*(8000*pow(MASS,9) +
    	  25312*pow(MASS,8)*r - 22664*pow(MASS,7)*pow(r,2) - 724*pow(MASS,6)*
    	  pow(r,3) + 640*pow(MASS,5)*pow(r,4) + 1090*pow(MASS,4)*pow(r,5) -
    	  180*pow(MASS,3)*pow(r,6) + 150*pow(MASS,2)*pow(r,7) - 15*MASS*pow(r,8)
    	  + 15*pow(r,9)) + MASS*(299880000*pow(MASS,9) - 553896000*pow(MASS,8)*r +
    	  404007800*pow(MASS,7)*pow(r,2) - 138548200*pow(MASS,6)*pow(r,3) +
    	  42434430*pow(MASS,5)*pow(r,4) - 26924080*pow(MASS,4)*pow(r,5) + 7433843*
    	  pow(MASS,3)*pow(r,6) + 357018*pow(MASS,2)*pow(r,7) + 224196*MASS*
    	  pow(r,8) - 187446*pow(r,9))*(-1 + 3*pow(cos(Th),2))))/(110250.*BETA*
    	  KAPPA*pow(MASS,4)*pow(2*MASS - r,3)*pow(r,9)));
    
    Th_d = (2*P_Th)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - (pow(ALPHA3,2)*
        (8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 3416700*pow(MASS,5)*
        pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 887110*pow(MASS,3)*pow(r,4) + 
          800733*pow(MASS,2)*pow(r,5) + 435540*MASS*pow(r,6) + 187446*pow(r,7))*
          pow(SPIN,2)*(1 + 3*cos(2*Th)))/(220500.*BETA*KAPPA*pow(MASS,3)*pow(r,8)));
    
    Ph_d = (220500*BETA*KAPPA*pow(MASS,3)*pow(r,8)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))*
     (JZ*(220500*BETA*KAPPA*pow(MASS,4)*pow(r,12) - 110250*BETA*KAPPA*pow(MASS,3)*pow(r,11)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
          7350*pow(ALPHA3,2)*pow(MASS,2)*pow(r,4)*(400*pow(MASS,4) - 96*pow(MASS,3)*r - 66*pow(MASS,2)*pow(r,2) - 130*MASS*pow(r,3) - 5*pow(r,4))*
           (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + pow(ALPHA3,2)*pow(SPIN,2)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))*
           (980000*pow(MASS,7)*r - 13798400*pow(MASS,6)*pow(r,2) + 2660700*pow(MASS,5)*pow(r,3) + 1249500*pow(MASS,4)*pow(r,4) + 
             1487150*pow(MASS,3)*pow(r,5) + 191100*pow(MASS,2)*pow(r,6) + 257250*MASS*pow(r,7) + 18375*pow(r,8) + 
             (8820000*pow(MASS,8) - 8173200*pow(MASS,7)*r - 15803900*pow(MASS,6)*pow(r,2) + 4198950*pow(MASS,5)*pow(r,3) + 
                4061710*pow(MASS,4)*pow(r,4) + 2275145*pow(MASS,3)*pow(r,5) - 164874*pow(MASS,2)*pow(r,6) - 187446*MASS*pow(r,7) - 187446*pow(r,8))*
              (-1 + 3*pow(cos(Th),2)))) - 7350*ENERGY*pow(MASS,2)*pow(r,4)*SPIN*
        (30*BETA*KAPPA*pow(MASS,2)*pow(r,8) + pow(ALPHA3,2)*(400*pow(MASS,4) - 144*pow(MASS,3)*r - 90*pow(MASS,2)*pow(r,2) - 140*MASS*pow(r,3) - 
             9*pow(r,4))*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*pow(sin(Th),2)))/
   (-54022500*pow(MASS,4)*pow(r,5)*pow(SPIN,2)*pow(30*BETA*KAPPA*pow(MASS,2)*pow(r,8) + 
        pow(ALPHA3,2)*(400*pow(MASS,4) - 144*pow(MASS,3)*r - 90*pow(MASS,2)*pow(r,2) - 140*MASS*pow(r,3) - 9*pow(r,4))*
         (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2)*pow(sin(Th),4) + 
     (220500*BETA*KAPPA*pow(MASS,4)*pow(r,12) - 110250*BETA*KAPPA*pow(MASS,3)*pow(r,11)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
        7350*pow(ALPHA3,2)*pow(MASS,2)*pow(r,4)*(400*pow(MASS,4) - 96*pow(MASS,3)*r - 66*pow(MASS,2)*pow(r,2) - 130*MASS*pow(r,3) - 5*pow(r,4))*
         (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + pow(ALPHA3,2)*pow(SPIN,2)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))*
         (980000*pow(MASS,7)*r - 13798400*pow(MASS,6)*pow(r,2) + 2660700*pow(MASS,5)*pow(r,3) + 1249500*pow(MASS,4)*pow(r,4) + 
           1487150*pow(MASS,3)*pow(r,5) + 191100*pow(MASS,2)*pow(r,6) + 257250*MASS*pow(r,7) + 18375*pow(r,8) + 
           (8820000*pow(MASS,8) - 8173200*pow(MASS,7)*r - 15803900*pow(MASS,6)*pow(r,2) + 4198950*pow(MASS,5)*pow(r,3) + 4061710*pow(MASS,4)*pow(r,4) + 
              2275145*pow(MASS,3)*pow(r,5) - 164874*pow(MASS,2)*pow(r,6) - 187446*MASS*pow(r,7) - 187446*pow(r,8))*(-1 + 3*pow(cos(Th),2))))*
      pow(sin(Th),2)*(pow(ALPHA3,2)*(8820000*pow(MASS,7) - 6213200*pow(MASS,6)*r - 3416700*pow(MASS,5)*pow(r,2) - 1855650*pow(MASS,4)*pow(r,3) + 
           887110*pow(MASS,3)*pow(r,4) + 800733*pow(MASS,2)*pow(r,5) + 435540*MASS*pow(r,6) + 187446*pow(r,7))*pow(SPIN,2)*(1 - 3*pow(cos(Th),2))*
         (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 110250*BETA*KAPPA*pow(MASS,3)*pow(r,8)*
         (pow(SPIN,2)*(pow(r,2) + pow(SPIN,2))*pow(cos(Th),2) + r*(pow(r,3) + r*pow(SPIN,2) + 2*MASS*pow(SPIN,2)*pow(sin(Th),2)))));
    
    P_r_d = (pow(P_Th,2)*(2*r + (4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
             (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
             (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2)))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - 
        (4463*pow(ALPHA3,2)*((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
             (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
             (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2)))/(2625.*BETA*KAPPA*pow(MASS,3)*r)))/
    pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - (4463*pow(ALPHA3,2)*
         (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
           (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*
         pow(SPIN,2)*(-1 + 3*pow(cos(Th),2)))/(2625.*BETA*KAPPA*pow(MASS,3)*r),2) + 
   (pow(P_r,2)*((4*pow(ALPHA3,2)*(1 - (368*pow(MASS,5))/(3.*pow(r,5)) + (16*pow(MASS,4))/(5.*pow(r,4)) + (2*pow(MASS,3))/pow(r,3) + 
             (52*pow(MASS,2))/(3.*pow(r,2)) + MASS/r))/(BETA*KAPPA*MASS*pow(1 - (2*MASS)/r,3)*pow(r,4)) + 
        (2*pow(ALPHA3,2)*(1 - (368*pow(MASS,5))/(3.*pow(r,5)) + (16*pow(MASS,4))/(5.*pow(r,4)) + (2*pow(MASS,3))/pow(r,3) + 
             (52*pow(MASS,2))/(3.*pow(r,2)) + MASS/r))/(BETA*KAPPA*pow(MASS,2)*pow(1 - (2*MASS)/r,2)*pow(r,3)) - 
        (pow(ALPHA3,2)*((1840*pow(MASS,5))/(3.*pow(r,6)) - (64*pow(MASS,4))/(5.*pow(r,5)) - (6*pow(MASS,3))/pow(r,4) - 
             (104*pow(MASS,2))/(3.*pow(r,3)) - MASS/pow(r,2)))/(BETA*KAPPA*pow(MASS,2)*pow(1 - (2*MASS)/r,2)*pow(r,2)) + 
        (2*r)/(-2*MASS*r + pow(r,2) + pow(SPIN,2)) - ((-2*MASS + 2*r)*(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))/
         pow(-2*MASS*r + pow(r,2) + pow(SPIN,2),2) - (pow(ALPHA3,2)*pow(SPIN,2)*
           (-(1 + (1600*pow(MASS,9))/(3.*pow(r,9)) + (25312*pow(MASS,8))/(15.*pow(r,8)) - (22664*pow(MASS,7))/(15.*pow(r,7)) - 
                 (724*pow(MASS,6))/(15.*pow(r,6)) + (128*pow(MASS,5))/(3.*pow(r,5)) + (218*pow(MASS,4))/(3.*pow(r,4)) - (12*pow(MASS,3))/pow(r,3) + 
                 (10*pow(MASS,2))/pow(r,2) - MASS/r)/(2.*MASS) - 
             (((-4800*pow(MASS,9))/pow(r,10) - (202496*pow(MASS,8))/(15.*pow(r,9)) + (158648*pow(MASS,7))/(15.*pow(r,8)) + 
                  (1448*pow(MASS,6))/(5.*pow(r,7)) - (640*pow(MASS,5))/(3.*pow(r,6)) - (872*pow(MASS,4))/(3.*pow(r,5)) + (36*pow(MASS,3))/pow(r,4) - 
                  (20*pow(MASS,2))/pow(r,3) + MASS/pow(r,2))*r)/(2.*MASS) + 
             (4463*((64260000*pow(MASS,9))/(4463.*pow(r,10)) - (105504000*pow(MASS,8))/(4463.*pow(r,9)) + (202003900*pow(MASS,7))/(13389.*pow(r,8)) - 
                  (19792600*pow(MASS,6))/(4463.*pow(r,7)) + (35362025*pow(MASS,5))/(31241.*pow(r,6)) - (53848160*pow(MASS,4))/(93723.*pow(r,5)) + 
                  (7433843*pow(MASS,3))/(62482.*pow(r,4)) + (119006*pow(MASS,2))/(31241.*pow(r,3)) + (5338*MASS)/(4463.*pow(r,2)))*(-1 + 3*pow(cos(Th),2)))/
              2625.))/(BETA*KAPPA*pow(MASS,3)*pow(1 - (2*MASS)/r,3)*pow(r,3)) + 
        (6*pow(ALPHA3,2)*pow(SPIN,2)*(-((1 + (1600*pow(MASS,9))/(3.*pow(r,9)) + (25312*pow(MASS,8))/(15.*pow(r,8)) - 
                   (22664*pow(MASS,7))/(15.*pow(r,7)) - (724*pow(MASS,6))/(15.*pow(r,6)) + (128*pow(MASS,5))/(3.*pow(r,5)) + 
                   (218*pow(MASS,4))/(3.*pow(r,4)) - (12*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) - MASS/r)*r)/(2.*MASS) + 
             (4463*(1 - (7140000*pow(MASS,9))/(4463.*pow(r,9)) + (13188000*pow(MASS,8))/(4463.*pow(r,8)) - (28857700*pow(MASS,7))/(13389.*pow(r,7)) + 
                  (9896300*pow(MASS,6))/(13389.*pow(r,6)) - (7072405*pow(MASS,5))/(31241.*pow(r,5)) + (13462040*pow(MASS,4))/(93723.*pow(r,4)) - 
                  (7433843*pow(MASS,3))/(187446.*pow(r,3)) - (59503*pow(MASS,2))/(31241.*pow(r,2)) - (5338*MASS)/(4463.*r))*(-1 + 3*pow(cos(Th),2)))/2625.))/
         (BETA*KAPPA*pow(MASS,2)*pow(1 - (2*MASS)/r,4)*pow(r,5)) + 
        (3*pow(ALPHA3,2)*pow(SPIN,2)*(-((1 + (1600*pow(MASS,9))/(3.*pow(r,9)) + (25312*pow(MASS,8))/(15.*pow(r,8)) - 
                   (22664*pow(MASS,7))/(15.*pow(r,7)) - (724*pow(MASS,6))/(15.*pow(r,6)) + (128*pow(MASS,5))/(3.*pow(r,5)) + 
                   (218*pow(MASS,4))/(3.*pow(r,4)) - (12*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) - MASS/r)*r)/(2.*MASS) + 
             (4463*(1 - (7140000*pow(MASS,9))/(4463.*pow(r,9)) + (13188000*pow(MASS,8))/(4463.*pow(r,8)) - (28857700*pow(MASS,7))/(13389.*pow(r,7)) + 
                  (9896300*pow(MASS,6))/(13389.*pow(r,6)) - (7072405*pow(MASS,5))/(31241.*pow(r,5)) + (13462040*pow(MASS,4))/(93723.*pow(r,4)) - 
                  (7433843*pow(MASS,3))/(187446.*pow(r,3)) - (59503*pow(MASS,2))/(31241.*pow(r,2)) - (5338*MASS)/(4463.*r))*(-1 + 3*pow(cos(Th),2)))/2625.))/
         (BETA*KAPPA*pow(MASS,3)*pow(1 - (2*MASS)/r,3)*pow(r,4))))/
    pow(-((pow(ALPHA3,2)*(1 - (368*pow(MASS,5))/(3.*pow(r,5)) + (16*pow(MASS,4))/(5.*pow(r,4)) + (2*pow(MASS,3))/pow(r,3) + 
             (52*pow(MASS,2))/(3.*pow(r,2)) + MASS/r))/(BETA*KAPPA*pow(MASS,2)*pow(1 - (2*MASS)/r,2)*pow(r,2))) + 
      (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))/(-2*MASS*r + pow(r,2) + pow(SPIN,2)) - 
      (pow(ALPHA3,2)*pow(SPIN,2)*(-((1 + (1600*pow(MASS,9))/(3.*pow(r,9)) + (25312*pow(MASS,8))/(15.*pow(r,8)) - (22664*pow(MASS,7))/(15.*pow(r,7)) - 
                 (724*pow(MASS,6))/(15.*pow(r,6)) + (128*pow(MASS,5))/(3.*pow(r,5)) + (218*pow(MASS,4))/(3.*pow(r,4)) - (12*pow(MASS,3))/pow(r,3) + 
                 (10*pow(MASS,2))/pow(r,2) - MASS/r)*r)/(2.*MASS) + 
           (4463*(1 - (7140000*pow(MASS,9))/(4463.*pow(r,9)) + (13188000*pow(MASS,8))/(4463.*pow(r,8)) - (28857700*pow(MASS,7))/(13389.*pow(r,7)) + 
                (9896300*pow(MASS,6))/(13389.*pow(r,6)) - (7072405*pow(MASS,5))/(31241.*pow(r,5)) + (13462040*pow(MASS,4))/(93723.*pow(r,4)) - 
                (7433843*pow(MASS,3))/(187446.*pow(r,3)) - (59503*pow(MASS,2))/(31241.*pow(r,2)) - (5338*MASS)/(4463.*r))*(-1 + 3*pow(cos(Th),2)))/2625.))/
       (BETA*KAPPA*pow(MASS,3)*pow(1 - (2*MASS)/r,3)*pow(r,3)),2) - 
   JZ*(-((JZ*(-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
              (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
             (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                  (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                  (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                  (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                     (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                     (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
           (-2*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                   pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + 
                (3*pow(ALPHA3,2)*((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*
                   SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) + 
                (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
                (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
              ((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                   pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
             (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
                 (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
                (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                     (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                     (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                     (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                        (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                        (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
              ((4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                     (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                     (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
                 (2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - (4463*pow(ALPHA3,2)*
                   ((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
                     (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
                     (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/(2625.*BETA*KAPPA*pow(MASS,3)*r) + 
                pow(sin(Th),2)*(2*r - (4*MASS*pow(r,2)*pow(SPIN,2)*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + 
                   (2*MASS*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
             ((pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
                 (BETA*KAPPA*MASS*pow(r,4)) - (pow(ALPHA3,2)*((320*pow(MASS,4))/pow(r,5) - (288*pow(MASS,3))/(5.*pow(r,4)) - 
                     (132*pow(MASS,2))/(5.*pow(r,3)) - (26*MASS)/pow(r,2)))/(3.*BETA*KAPPA*MASS*pow(r,3)) - 
                (4*MASS*pow(r,2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + (2*MASS)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
                (4463*pow(ALPHA3,2)*pow(SPIN,2)*((490000*pow(MASS,7))/(13389.*pow(r,8)) - (1971200*pow(MASS,6))/(4463.*pow(r,7)) + 
                     (316750*pow(MASS,5))/(4463.*pow(r,6)) + (119000*pow(MASS,4))/(4463.*pow(r,5)) + (106225*pow(MASS,3))/(4463.*pow(r,4)) + 
                     (9100*pow(MASS,2))/(4463.*pow(r,3)) + (6125*MASS)/(4463.*pow(r,2)) + 
                     ((1680000*pow(MASS,8))/(4463.*pow(r,9)) - (1362200*pow(MASS,7))/(4463.*pow(r,8)) - (2257700*pow(MASS,6))/(4463.*pow(r,7)) + 
                        (499875*pow(MASS,5))/(4463.*pow(r,6)) + (8123420*pow(MASS,4))/(93723.*pow(r,5)) + (2275145*pow(MASS,3))/(62482.*pow(r,4)) - 
                        (54958*pow(MASS,2))/(31241.*pow(r,3)) - MASS/pow(r,2))*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + 
                (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                     (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                     (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                     (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                        (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                        (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(875.*BETA*KAPPA*pow(MASS,3)*pow(r,4)))*
              ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                     (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                     (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
                 (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
                 (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))))/
         pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
                 SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2)) - 
      (ENERGY*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
         (-2*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + 
              (3*pow(ALPHA3,2)*((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*
                 SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) + 
              (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
              (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - (4463*pow(ALPHA3,2)*
                 ((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
                   (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
                   (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/(2625.*BETA*KAPPA*pow(MASS,3)*r) + 
              pow(sin(Th),2)*(2*r - (4*MASS*pow(r,2)*pow(SPIN,2)*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + 
                 (2*MASS*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           ((pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (BETA*KAPPA*MASS*pow(r,4)) - (pow(ALPHA3,2)*((320*pow(MASS,4))/pow(r,5) - (288*pow(MASS,3))/(5.*pow(r,4)) - 
                   (132*pow(MASS,2))/(5.*pow(r,3)) - (26*MASS)/pow(r,2)))/(3.*BETA*KAPPA*MASS*pow(r,3)) - 
              (4*MASS*pow(r,2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + (2*MASS)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*((490000*pow(MASS,7))/(13389.*pow(r,8)) - (1971200*pow(MASS,6))/(4463.*pow(r,7)) + 
                   (316750*pow(MASS,5))/(4463.*pow(r,6)) + (119000*pow(MASS,4))/(4463.*pow(r,5)) + (106225*pow(MASS,3))/(4463.*pow(r,4)) + 
                   (9100*pow(MASS,2))/(4463.*pow(r,3)) + (6125*MASS)/(4463.*pow(r,2)) + 
                   ((1680000*pow(MASS,8))/(4463.*pow(r,9)) - (1362200*pow(MASS,7))/(4463.*pow(r,8)) - (2257700*pow(MASS,6))/(4463.*pow(r,7)) + 
                      (499875*pow(MASS,5))/(4463.*pow(r,6)) + (8123420*pow(MASS,4))/(93723.*pow(r,5)) + (2275145*pow(MASS,3))/(62482.*pow(r,4)) - 
                      (54958*pow(MASS,2))/(31241.*pow(r,3)) - MASS/pow(r,2))*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(875.*BETA*KAPPA*pow(MASS,3)*pow(r,4)))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2) + 
      (JZ*((pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
            (BETA*KAPPA*MASS*pow(r,4)) - (pow(ALPHA3,2)*((320*pow(MASS,4))/pow(r,5) - (288*pow(MASS,3))/(5.*pow(r,4)) - 
                (132*pow(MASS,2))/(5.*pow(r,3)) - (26*MASS)/pow(r,2)))/(3.*BETA*KAPPA*MASS*pow(r,3)) - 
           (4*MASS*pow(r,2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + (2*MASS)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
           (4463*pow(ALPHA3,2)*pow(SPIN,2)*((490000*pow(MASS,7))/(13389.*pow(r,8)) - (1971200*pow(MASS,6))/(4463.*pow(r,7)) + 
                (316750*pow(MASS,5))/(4463.*pow(r,6)) + (119000*pow(MASS,4))/(4463.*pow(r,5)) + (106225*pow(MASS,3))/(4463.*pow(r,4)) + 
                (9100*pow(MASS,2))/(4463.*pow(r,3)) + (6125*MASS)/(4463.*pow(r,2)) + 
                ((1680000*pow(MASS,8))/(4463.*pow(r,9)) - (1362200*pow(MASS,7))/(4463.*pow(r,8)) - (2257700*pow(MASS,6))/(4463.*pow(r,7)) + 
                   (499875*pow(MASS,5))/(4463.*pow(r,6)) + (8123420*pow(MASS,4))/(93723.*pow(r,5)) + (2275145*pow(MASS,3))/(62482.*pow(r,4)) - 
                   (54958*pow(MASS,2))/(31241.*pow(r,3)) - MASS/pow(r,2))*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + 
           (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                   (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(875.*BETA*KAPPA*pow(MASS,3)*pow(r,4))))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))) + 
      (ENERGY*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + (3*pow(ALPHA3,2)*
              ((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*SPIN*pow(sin(Th),2))/
            (5.*BETA*KAPPA*MASS*pow(r,3)) + (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
           (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))))) + 
   ENERGY*((JZ*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
         (-2*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + 
              (3*pow(ALPHA3,2)*((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*
                 SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) + 
              (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
              (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - (4463*pow(ALPHA3,2)*
                 ((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
                   (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
                   (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/(2625.*BETA*KAPPA*pow(MASS,3)*r) + 
              pow(sin(Th),2)*(2*r - (4*MASS*pow(r,2)*pow(SPIN,2)*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + 
                 (2*MASS*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           ((pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (BETA*KAPPA*MASS*pow(r,4)) - (pow(ALPHA3,2)*((320*pow(MASS,4))/pow(r,5) - (288*pow(MASS,3))/(5.*pow(r,4)) - 
                   (132*pow(MASS,2))/(5.*pow(r,3)) - (26*MASS)/pow(r,2)))/(3.*BETA*KAPPA*MASS*pow(r,3)) - 
              (4*MASS*pow(r,2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + (2*MASS)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*((490000*pow(MASS,7))/(13389.*pow(r,8)) - (1971200*pow(MASS,6))/(4463.*pow(r,7)) + 
                   (316750*pow(MASS,5))/(4463.*pow(r,6)) + (119000*pow(MASS,4))/(4463.*pow(r,5)) + (106225*pow(MASS,3))/(4463.*pow(r,4)) + 
                   (9100*pow(MASS,2))/(4463.*pow(r,3)) + (6125*MASS)/(4463.*pow(r,2)) + 
                   ((1680000*pow(MASS,8))/(4463.*pow(r,9)) - (1362200*pow(MASS,7))/(4463.*pow(r,8)) - (2257700*pow(MASS,6))/(4463.*pow(r,7)) + 
                      (499875*pow(MASS,5))/(4463.*pow(r,6)) + (8123420*pow(MASS,4))/(93723.*pow(r,5)) + (2275145*pow(MASS,3))/(62482.*pow(r,4)) - 
                      (54958*pow(MASS,2))/(31241.*pow(r,3)) - MASS/pow(r,2))*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(875.*BETA*KAPPA*pow(MASS,3)*pow(r,4)))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2) + 
      (ENERGY*((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
              (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))*
         (-2*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + 
              (3*pow(ALPHA3,2)*((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*
                 SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) + 
              (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
              (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - (4463*pow(ALPHA3,2)*
                 ((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
                   (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
                   (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/(2625.*BETA*KAPPA*pow(MASS,3)*r) + 
              pow(sin(Th),2)*(2*r - (4*MASS*pow(r,2)*pow(SPIN,2)*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + 
                 (2*MASS*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           ((pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (BETA*KAPPA*MASS*pow(r,4)) - (pow(ALPHA3,2)*((320*pow(MASS,4))/pow(r,5) - (288*pow(MASS,3))/(5.*pow(r,4)) - 
                   (132*pow(MASS,2))/(5.*pow(r,3)) - (26*MASS)/pow(r,2)))/(3.*BETA*KAPPA*MASS*pow(r,3)) - 
              (4*MASS*pow(r,2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + (2*MASS)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*((490000*pow(MASS,7))/(13389.*pow(r,8)) - (1971200*pow(MASS,6))/(4463.*pow(r,7)) + 
                   (316750*pow(MASS,5))/(4463.*pow(r,6)) + (119000*pow(MASS,4))/(4463.*pow(r,5)) + (106225*pow(MASS,3))/(4463.*pow(r,4)) + 
                   (9100*pow(MASS,2))/(4463.*pow(r,3)) + (6125*MASS)/(4463.*pow(r,2)) + 
                   ((1680000*pow(MASS,8))/(4463.*pow(r,9)) - (1362200*pow(MASS,7))/(4463.*pow(r,8)) - (2257700*pow(MASS,6))/(4463.*pow(r,7)) + 
                      (499875*pow(MASS,5))/(4463.*pow(r,6)) + (8123420*pow(MASS,4))/(93723.*pow(r,5)) + (2275145*pow(MASS,3))/(62482.*pow(r,4)) - 
                      (54958*pow(MASS,2))/(31241.*pow(r,3)) - MASS/pow(r,2))*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(875.*BETA*KAPPA*pow(MASS,3)*pow(r,4)))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2) - 
      (JZ*((-9*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,4)) + (3*pow(ALPHA3,2)*
              ((1600*pow(MASS,4))/(9.*pow(r,5)) - (48*pow(MASS,3))/pow(r,4) - (20*pow(MASS,2))/pow(r,3) - (140*MASS)/(9.*pow(r,2)))*SPIN*pow(sin(Th),2))/
            (5.*BETA*KAPPA*MASS*pow(r,3)) + (4*MASS*pow(r,2)*SPIN*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) - 
           (2*MASS*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))) - 
      (ENERGY*((4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (2625.*BETA*KAPPA*pow(MASS,3)*pow(r,2)) - (4463*pow(ALPHA3,2)*
              ((-1470000*pow(MASS,7))/(4463.*pow(r,8)) + (887600*pow(MASS,6))/(4463.*pow(r,7)) + (406750*pow(MASS,5))/(4463.*pow(r,6)) + 
                (1237100*pow(MASS,4))/(31241.*pow(r,5)) - (63365*pow(MASS,3))/(4463.*pow(r,4)) - (266911*pow(MASS,2))/(31241.*pow(r,3)) - 
                (10370*MASS)/(4463.*pow(r,2)))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/(2625.*BETA*KAPPA*pow(MASS,3)*r) + 
           pow(sin(Th),2)*(2*r - (4*MASS*pow(r,2)*pow(SPIN,2)*pow(sin(Th),2))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2) + 
              (2*MASS*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))));
    
    P_Th_d = (pow(P_Th,2)*(-2*pow(SPIN,2)*cos(Th)*sin(Th) + (8926*pow(ALPHA3,2)*
           (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
             (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r)
             )*pow(SPIN,2)*cos(Th)*sin(Th))/(875.*BETA*KAPPA*pow(MASS,3)*r)))/
    pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2) - (4463*pow(ALPHA3,2)*
         (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
           (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*
         pow(SPIN,2)*(-1 + 3*pow(cos(Th),2)))/(2625.*BETA*KAPPA*pow(MASS,3)*r),2) + 
   (pow(P_r,2)*((8926*pow(ALPHA3,2)*(1 - (7140000*pow(MASS,9))/(4463.*pow(r,9)) + (13188000*pow(MASS,8))/(4463.*pow(r,8)) - 
             (28857700*pow(MASS,7))/(13389.*pow(r,7)) + (9896300*pow(MASS,6))/(13389.*pow(r,6)) - (7072405*pow(MASS,5))/(31241.*pow(r,5)) + 
             (13462040*pow(MASS,4))/(93723.*pow(r,4)) - (7433843*pow(MASS,3))/(187446.*pow(r,3)) - (59503*pow(MASS,2))/(31241.*pow(r,2)) - 
             (5338*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*sin(Th))/(875.*BETA*KAPPA*pow(MASS,3)*pow(1 - (2*MASS)/r,3)*pow(r,3)) - 
        (2*pow(SPIN,2)*cos(Th)*sin(Th))/(-2*MASS*r + pow(r,2) + pow(SPIN,2))))/
    pow(-((pow(ALPHA3,2)*(1 - (368*pow(MASS,5))/(3.*pow(r,5)) + (16*pow(MASS,4))/(5.*pow(r,4)) + (2*pow(MASS,3))/pow(r,3) + 
             (52*pow(MASS,2))/(3.*pow(r,2)) + MASS/r))/(BETA*KAPPA*pow(MASS,2)*pow(1 - (2*MASS)/r,2)*pow(r,2))) + 
      (pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))/(-2*MASS*r + pow(r,2) + pow(SPIN,2)) - 
      (pow(ALPHA3,2)*pow(SPIN,2)*(-((1 + (1600*pow(MASS,9))/(3.*pow(r,9)) + (25312*pow(MASS,8))/(15.*pow(r,8)) - (22664*pow(MASS,7))/(15.*pow(r,7)) - 
                 (724*pow(MASS,6))/(15.*pow(r,6)) + (128*pow(MASS,5))/(3.*pow(r,5)) + (218*pow(MASS,4))/(3.*pow(r,4)) - (12*pow(MASS,3))/pow(r,3) + 
                 (10*pow(MASS,2))/pow(r,2) - MASS/r)*r)/(2.*MASS) + 
           (4463*(1 - (7140000*pow(MASS,9))/(4463.*pow(r,9)) + (13188000*pow(MASS,8))/(4463.*pow(r,8)) - (28857700*pow(MASS,7))/(13389.*pow(r,7)) + 
                (9896300*pow(MASS,6))/(13389.*pow(r,6)) - (7072405*pow(MASS,5))/(31241.*pow(r,5)) + (13462040*pow(MASS,4))/(93723.*pow(r,4)) - 
                (7433843*pow(MASS,3))/(187446.*pow(r,3)) - (59503*pow(MASS,2))/(31241.*pow(r,2)) - (5338*MASS)/(4463.*r))*(-1 + 3*pow(cos(Th),2)))/2625.))/
       (BETA*KAPPA*pow(MASS,3)*pow(1 - (2*MASS)/r,3)*pow(r,3)),2) - 
   JZ*((JZ*((8926*pow(ALPHA3,2)*(1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + 
                (1128850*pow(MASS,6))/(13389.*pow(r,6)) - (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - 
                (2275145*pow(MASS,3))/(187446.*pow(r,3)) + (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*pow(SPIN,2)*cos(Th)*sin(Th))/
            (875.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + (4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))) + 
      (ENERGY*((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*cos(Th)*
              sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
           (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))) - 
      (JZ*(-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
            (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
           (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                   (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
         (-2*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 cos(Th)*sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)) + 
           ((8926*pow(ALPHA3,2)*(1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + 
                   (1128850*pow(MASS,6))/(13389.*pow(r,6)) - (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - 
                   (2275145*pow(MASS,3))/(187446.*pow(r,3)) + (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*pow(SPIN,2)*cos(Th)*sin(Th))/
               (875.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + (4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-8926*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*(-1 + 3*pow(cos(Th),2))*sin(Th))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + (8926*pow(ALPHA3,2)*
                 (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
                   (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + 
                   (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*pow(sin(Th),3))/(875.*BETA*KAPPA*pow(MASS,3)*r) + 
              2*cos(Th)*sin(Th)*(pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
              pow(sin(Th),2)*((4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
                 (4*MASS*r*pow(SPIN,4)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2) - 
      (ENERGY*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
         (-2*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 cos(Th)*sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)) + 
           ((8926*pow(ALPHA3,2)*(1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + 
                   (1128850*pow(MASS,6))/(13389.*pow(r,6)) - (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - 
                   (2275145*pow(MASS,3))/(187446.*pow(r,3)) + (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*pow(SPIN,2)*cos(Th)*sin(Th))/
               (875.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + (4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-8926*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*(-1 + 3*pow(cos(Th),2))*sin(Th))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + (8926*pow(ALPHA3,2)*
                 (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
                   (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + 
                   (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*pow(sin(Th),3))/(875.*BETA*KAPPA*pow(MASS,3)*r) + 
              2*cos(Th)*sin(Th)*(pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
              pow(sin(Th),2)*((4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
                 (4*MASS*r*pow(SPIN,4)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2)) + 
   ENERGY*(-((JZ*((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                cos(Th)*sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
             (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))/
         (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))))) - 
      (ENERGY*((-8926*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*(-1 + 3*pow(cos(Th),2))*sin(Th))/
            (2625.*BETA*KAPPA*pow(MASS,3)*r) + (8926*pow(ALPHA3,2)*
              (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
                (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + 
                (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*pow(sin(Th),3))/(875.*BETA*KAPPA*pow(MASS,3)*r) + 
           2*cos(Th)*sin(Th)*(pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
           pow(sin(Th),2)*((4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
              (4*MASS*r*pow(SPIN,4)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2))))/
       (-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
               pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))) + 
      (JZ*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
              pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
         (-2*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 cos(Th)*sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)) + 
           ((8926*pow(ALPHA3,2)*(1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + 
                   (1128850*pow(MASS,6))/(13389.*pow(r,6)) - (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - 
                   (2275145*pow(MASS,3))/(187446.*pow(r,3)) + (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*pow(SPIN,2)*cos(Th)*sin(Th))/
               (875.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + (4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-8926*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*(-1 + 3*pow(cos(Th),2))*sin(Th))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + (8926*pow(ALPHA3,2)*
                 (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
                   (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + 
                   (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*pow(sin(Th),3))/(875.*BETA*KAPPA*pow(MASS,3)*r) + 
              2*cos(Th)*sin(Th)*(pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
              pow(sin(Th),2)*((4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
                 (4*MASS*r*pow(SPIN,4)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2) + 
      (ENERGY*((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
            (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
              (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))))*
         (-2*((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))*
            ((6*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*SPIN*
                 cos(Th)*sin(Th))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (4*MASS*r*SPIN*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4*MASS*r*pow(SPIN,3)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)) + 
           ((8926*pow(ALPHA3,2)*(1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + 
                   (1128850*pow(MASS,6))/(13389.*pow(r,6)) - (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - 
                   (2275145*pow(MASS,3))/(187446.*pow(r,3)) + (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*pow(SPIN,2)*cos(Th)*sin(Th))/
               (875.*BETA*KAPPA*pow(MASS,3)*pow(r,3)) + (4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2))*
            ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*
               (pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))) + 
           (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
               (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
              (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                   (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                   (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                      (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                      (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
            ((-8926*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                   (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                   (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*(-1 + 3*pow(cos(Th),2))*sin(Th))/
               (2625.*BETA*KAPPA*pow(MASS,3)*r) + (8926*pow(ALPHA3,2)*
                 (1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - (81350*pow(MASS,5))/(4463.*pow(r,5)) - 
                   (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + (266911*pow(MASS,2))/(62482.*pow(r,2)) + 
                   (10370*MASS)/(4463.*r))*pow(SPIN,2)*cos(Th)*pow(sin(Th),3))/(875.*BETA*KAPPA*pow(MASS,3)*r) + 
              2*cos(Th)*sin(Th)*(pow(r,2) + pow(SPIN,2) + (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2))) + 
              pow(sin(Th),2)*((4*MASS*r*pow(SPIN,2)*cos(Th)*sin(Th))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) + 
                 (4*MASS*r*pow(SPIN,4)*cos(Th)*pow(sin(Th),3))/pow(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2),2)))))/
       pow(-pow((3*pow(ALPHA3,2)*(1 - (400*pow(MASS,4))/(9.*pow(r,4)) + (16*pow(MASS,3))/pow(r,3) + (10*pow(MASS,2))/pow(r,2) + (140*MASS)/(9.*r))*
               SPIN*pow(sin(Th),2))/(5.*BETA*KAPPA*MASS*pow(r,3)) - (2*MASS*r*SPIN*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)),2) + 
         (-1 - (pow(ALPHA3,2)*(1 - (80*pow(MASS,4))/pow(r,4) + (96*pow(MASS,3))/(5.*pow(r,3)) + (66*pow(MASS,2))/(5.*pow(r,2)) + (26*MASS)/r))/
             (3.*BETA*KAPPA*MASS*pow(r,3)) + (2*MASS*r)/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)) - 
            (4463*pow(ALPHA3,2)*pow(SPIN,2)*(-0.09802823213085368 - (70000*pow(MASS,7))/(13389.*pow(r,7)) + (985600*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (63350*pow(MASS,5))/(4463.*pow(r,5)) - (29750*pow(MASS,4))/(4463.*pow(r,4)) - (106225*pow(MASS,3))/(13389.*pow(r,3)) - 
                 (4550*pow(MASS,2))/(4463.*pow(r,2)) - (6125*MASS)/(4463.*r) + 
                 (1 - (210000*pow(MASS,8))/(4463.*pow(r,8)) + (194600*pow(MASS,7))/(4463.*pow(r,7)) + (1128850*pow(MASS,6))/(13389.*pow(r,6)) - 
                    (99975*pow(MASS,5))/(4463.*pow(r,5)) - (2030855*pow(MASS,4))/(93723.*pow(r,4)) - (2275145*pow(MASS,3))/(187446.*pow(r,3)) + 
                    (27479*pow(MASS,2))/(31241.*pow(r,2)) + MASS/r)*(-1 + 3*pow(cos(Th),2))))/(2625.*BETA*KAPPA*pow(MASS,3)*pow(r,3)))*
          ((-4463*pow(ALPHA3,2)*(1 + (210000*pow(MASS,7))/(4463.*pow(r,7)) - (443800*pow(MASS,6))/(13389.*pow(r,6)) - 
                 (81350*pow(MASS,5))/(4463.*pow(r,5)) - (309275*pow(MASS,4))/(31241.*pow(r,4)) + (63365*pow(MASS,3))/(13389.*pow(r,3)) + 
                 (266911*pow(MASS,2))/(62482.*pow(r,2)) + (10370*MASS)/(4463.*r))*pow(SPIN,2)*(-1 + 3*pow(cos(Th),2))*pow(sin(Th),2))/
             (2625.*BETA*KAPPA*pow(MASS,3)*r) + pow(sin(Th),2)*(pow(r,2) + pow(SPIN,2) + 
               (2*MASS*r*pow(SPIN,2)*pow(sin(Th),2))/(pow(r,2) + pow(SPIN,2)*pow(cos(Th),2)))),2));
        
    P_Ph_d = 0.;
    
    gsl_vector_set(out_state, 0, r_d);
    gsl_vector_set(out_state, 1, P_r_d);
    gsl_vector_set(out_state, 2, Th_d);
    gsl_vector_set(out_state, 3, P_Th_d);
    gsl_vector_set(out_state, 4, Ph_d);
    gsl_vector_set(out_state, 5, P_Ph_d); 
       
	return 0;
    }



