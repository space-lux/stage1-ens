#include <stdlib.h>
#include <math.h>

float drand() {
	return (float)rand()/RAND_MAX;
}

float randn(float mu,float sigma) {
	float r=0;
	int n=12;
	for(int i=0;i<n;i++) {
		r+=drand()-0.5;
	}
	r*=sigma;
	r+=mu;
	return r;
}
