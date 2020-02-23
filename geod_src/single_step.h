double single_stepRK4(int (*f) (double, gsl_vector *, gsl_vector *),
                          int dof,
                          gsl_vector * in_vec,
                          gsl_vector * out_vec,
                          double x0,
                          double h);
double single_stepRKF78(int (*f) (double, gsl_vector *, gsl_vector *),
		    int dof, gsl_vector * in_vec, gsl_vector * out_vec, double x0,
		    double h);