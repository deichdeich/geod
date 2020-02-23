 void integration_vector_init(int dof);
 void integration_vector_free();
 void integration_vector_zero(int dof);
 struct {
    gsl_vector *k1_in;
    gsl_vector *k2_in;
    gsl_vector *k3_in;
    gsl_vector *k4_in;
    gsl_vector *k5_in;
    gsl_vector *k6_in;
    gsl_vector *k7_in;
    gsl_vector *k8_in;
    gsl_vector *k9_in;
    gsl_vector *k10_in;
    gsl_vector *k11_in;
    gsl_vector *k12_in;
    gsl_vector *k13_in;
    gsl_vector *k1_out;
    gsl_vector *k2_out;
    gsl_vector *k3_out;
    gsl_vector *k4_out;
    gsl_vector *k5_out;
    gsl_vector *k6_out;
    gsl_vector *k7_out;
    gsl_vector *k8_out;
    gsl_vector *k9_out;
    gsl_vector *k10_out;
    gsl_vector *k11_out;
    gsl_vector *k12_out;
    gsl_vector *k13_out;

    gsl_vector *k3k1_out;
    gsl_vector *k3k2_out;

    gsl_vector *k4k1_out;
    gsl_vector *k4k3_out;

    gsl_vector *k5k1_out;
    gsl_vector *k5k3_out;
    gsl_vector *k5k4_out;

    gsl_vector *k6k1_out;
    gsl_vector *k6k4_out;
    gsl_vector *k6k5_out;

    gsl_vector *k7k1_out;
    gsl_vector *k7k4_out;
    gsl_vector *k7k5_out;
    gsl_vector *k7k6_out;

    gsl_vector *k8k1_out;
    gsl_vector *k8k5_out;
    gsl_vector *k8k6_out;
    gsl_vector *k8k7_out;

    gsl_vector *k9k1_out;
    gsl_vector *k9k4_out;
    gsl_vector *k9k5_out;
    gsl_vector *k9k6_out;
    gsl_vector *k9k7_out;
    gsl_vector *k9k8_out;

    gsl_vector *k10k1_out;
    gsl_vector *k10k4_out;
    gsl_vector *k10k5_out;
    gsl_vector *k10k6_out;
    gsl_vector *k10k7_out;
    gsl_vector *k10k8_out;
    gsl_vector *k10k9_out;

    gsl_vector *k11k1_out;
    gsl_vector *k11k4_out;
    gsl_vector *k11k5_out;
    gsl_vector *k11k6_out;
    gsl_vector *k11k7_out;
    gsl_vector *k11k8_out;
    gsl_vector *k11k9_out;
    gsl_vector *k11k10_out;

    gsl_vector *k12k1_out;
    gsl_vector *k12k6_out;
    gsl_vector *k12k7_out;
    gsl_vector *k12k8_out;
    gsl_vector *k12k9_out;
    gsl_vector *k12k10_out;

    gsl_vector *k13k1_out;
    gsl_vector *k13k4_out;
    gsl_vector *k13k5_out;
    gsl_vector *k13k6_out;
    gsl_vector *k13k7_out;
    gsl_vector *k13k8_out;
    gsl_vector *k13k9_out;
    gsl_vector *k13k10_out;
    gsl_vector *k13k12_out;

    gsl_vector *c_1_11_vec;
    gsl_vector *c_6_vec;
    gsl_vector *c_7_8_vec;
    gsl_vector *c_9_10_vec;
    gsl_vector *c_tot_vec;

    gsl_vector *err_vec;
    gsl_vector *ek12;
    gsl_vector *ek13;

    gsl_vector *temp_in_vec;
    gsl_vector *temp_out_vec;
    gsl_vector *zero_vec;

} integration_vectors;