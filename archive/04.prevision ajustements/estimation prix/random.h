#include <stdlib.h>
#include <math.h>

double drand() {
	return (double)rand()/RAND_MAX;
}

double randn(double mu,double sigma) {
	double r=0;
	int n=12;
	for(int i=0;i<n;i++) {
		r+=drand()-0.5;
	}
	r*=sigma;
	r+=mu;
	return r;
}
