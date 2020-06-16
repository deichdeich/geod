int double_pendulum(double t, gsl_vector * in_state, gsl_vector * out_state);
int test_sine(double t, gsl_vector * in_state, gsl_vector * out_state);
int test_high_pow(double t, gsl_vector * in_state, gsl_vector * out_state);
int polar2d_lagr(double t, gsl_vector * in_state, gsl_vector * out_state);
int polar2d_ham(double t, gsl_vector * in_state, gsl_vector * out_state);
int kerr_eom_SR_o2(double t, gsl_vector * in_state, gsl_vector * out_state);
int kerr_eom_exact(double t, gsl_vector * in_state, gsl_vector * out_state);
int edgb_a2_z1(double t, gsl_vector * in_state, gsl_vector * out_state);

static double csc(double arg);
static double cot(double arg);
