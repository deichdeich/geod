#include <stdio.h>
#include <math.h>

int main(){
    double da = 0.1;
    long double lda = 0.1;
    double fact = 1. / 9.;
    printf("%.2le\n", fact);
    printf("%.2le\n", pow(da, fact));
    printf("%.2le\n", da * pow(10, -25));
}
