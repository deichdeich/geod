int stepper(int (*f) (double, gsl_vector *, gsl_vector *),
	     int dof, gsl_vector * in_state_vec, gsl_vector * out_state_vec,
	     double x, double h, double xmax, double *h_next,
	     double tolerance);
int stepper2(int (*f) (double, gsl_vector *, gsl_vector *),
	     int dof,
	     gsl_vector * in_state_vec,
	     gsl_vector * out_state_vec,
	     double x,
	     double h,
	     double xmax,
	     double tolerance);
int check_err(double tol, int verbose);
void poincare_check(gsl_vector * state);
double poincare_check2(gsl_vector * state, double h);
double poincare_check_bisect(int (*f) (double, gsl_vector *, gsl_vector *), gsl_vector *state, double h, int dof, double x);
double bisect(int (*f) (double, gsl_vector *, gsl_vector *), double prev_dist, double h, gsl_vector *bisect_out, int dof, double x);
double dist1;
