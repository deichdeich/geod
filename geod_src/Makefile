geod : ad_rkf78.o eoms.o integration_vectors.o steppers.o fileio.o single_step.o get_eom.o
	gcc -w -o geod ad_rkf78.o eoms.o integration_vectors.o fileio.o single_step.o get_eom.o \
	steppers.o -lgslcblas -lgsl -lm

#rkf78test : rkf78test.o integration_vectors.o single_step.o
#	gcc -w -o rkf78test rkf78test.o integration_vectors.o single_step.o \
#	-lgslcblas -lgsl -lm

ad_rkf78.o : ad_rkf78.c ad_rkf78.h eoms.h integration_vectors.h single_step.h steppers.h \
	fileio.h definitions.h
	gcc -w -c -Ofast -I/usr/local/include ad_rkf78.c -lm

eoms.o : eoms.c ad_rkf78.h eoms.h integration_vectors.h single_step.h steppers.h definitions.h \
	fileio.h
	gcc -w -c -Ofast -I/usr/local/include eoms.c -lm

integration_vectors.o : integration_vectors.c ad_rkf78.h eoms.h integration_vectors.h \
	single_step.h steppers.h fileio.h definitions.h
	gcc -w -c -Ofast -I/usr/local/include integration_vectors.c -lm

steppers.o : steppers2.c ad_rkf78.h eoms.h integration_vectors.h single_step.h steppers.h \
	fileio.h definitions.h
	gcc -w -c -Ofast -I/usr/local/include steppers2.c -lm -o steppers.o
    
single_step.o : single_step.c ad_rkf78.h eoms.h integration_vectors.h single_step.h \
	steppers.h fileio.h definitions.h
	gcc -w -c -Ofast -I/usr/local/include single_step.c -lm

fileio.o : fileio.c ad_rkf78.h eoms.h integration_vectors.h single_step.h steppers.h \
	fileio.h definitions.h
	gcc -w -c -Ofast -I/usr/local/include fileio.c -lm

rkf78test.o : rkf78test.c single_step.h integration_vectors.h
	gcc -w -c -Ofast -I/usr/local/include rkf78test.c -lm

get_eom.o : ad_rkf78.h
	gcc -w -c -Ofast -I/usr/local/include get_eom.c -lm
    
clean :
	rm geod ad_rkf78.o eoms.o integration_vectors.o steppers.o fileio.o single_step.o
