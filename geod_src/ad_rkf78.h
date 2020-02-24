//////////////////////////
// GLOBALS
/////////////////////////
/* set these in ctypes: https://stackoverflow.com/questions/30338896/os-x-python-ctypes-global-variables
double MASS;
double SPIN;
double Jz;
double Ee;
*/

////////////////////////////////
//  FUNCTIONS
///////////////////////////////

int *vect_add(gsl_vector * sum, int vect_len, int count, ...);
int populate_history(int dof, int step, double clock, gsl_vector * state);
void write_history(int step);
int arr2vec(int len, double in_arr[len], gsl_vector * out_vec);
void print_vec(gsl_vector * vec);
gsl_vector * cur_state;
double (*int_func)(int (*f) (double, gsl_vector *, gsl_vector *),
                          int dof,
                          gsl_vector * in_vec,
                          gsl_vector * out_vec,
                          double x0,
                          double h);
int (*f) (double, gsl_vector *, gsl_vector *);
#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 10
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.
#define HIRES_SCALE_FACTOR 100
#define SLOWDOWN_DIST 10

#define hist_len 1000
#define GRAVITY 9.8
#define TPI 6.283185307179586

void get_eom(char eom_name[]);
int make_section;
int lines_written;
double poincare_condition[2];
double poincare_tolerance;
double history[hist_len][7];
_Bool poincare_yes;
_Bool add_to_history;

char *filename;
