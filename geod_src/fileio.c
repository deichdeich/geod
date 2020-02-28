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


void write_history(int step){
    int line;
    FILE *outfile = fopen(filename, "ab");
    for(line = 0; line <= step; line++){
        fwrite(history[line], sizeof(history[line][0]), 7, outfile);
        //printf("\t printing line %d\n", line);

    }
    fclose(outfile);
}

void clear_file(){
    FILE *outfile = fopen(filename, "w");
    fprintf(outfile, "");
    fclose(outfile);
}

int metadata(double t0, double tmax, double tol, double h, char int_func[]){
    FILE *outfile;
    outfile = fopen(filename, "w");
    fprintf(outfile, "# SIMULATION PARAMETERS\n");
    fprintf(outfile, "# t0: %0.1f\n", t0);
    fprintf(outfile, "# tmax: %0.1f\n", tmax);
    fprintf(outfile, "# tolerance: %0.10f\n", tol);
    fprintf(outfile, "# h: %0.10f\n", h);
    fprintf(outfile, "# integrator: %s\n", int_func);
    fprintf(outfile, "#\n");
    fprintf(outfile, "# PHYSICAL PARAMETERS\n");
    fprintf(outfile, "# M: %0.4f\n", MASS);
    fprintf(outfile, "# a: %0.4f\n", SPIN);
    fprintf(outfile, "# Jz: %0.4f\n", JZ);
    fprintf(outfile, "# energy: %0.4f\n", ENERGY);
    fprintf(outfile, "#\n");
    fclose(outfile);
return 0;
}
