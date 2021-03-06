#include <stdlib.h>
#include <mpi.h>

typedef struct{
	unsigned int len;
	double * data;
} vec;


void vec_zero(vec* v) {
	for(int i=0;i<v->len;i++) {
		v->data[i]=0;
	}
}

vec* vec_new(unsigned int len) {
	vec* v=malloc(sizeof(vec));
	v->data=malloc(len*sizeof(double));
	v->len=len;
	return v;
}

vec* vec_new_fromtab(double* tab,unsigned int len) {
	vec* v=malloc(sizeof(vec));
	v->data=malloc(len*sizeof(double));
	v->len=len;
	for(int i=0;i<v->len;i++) {
		v->data[i]=tab[i];
	}
	return v;
}

double vec_sum(vec* v) {
	double r=0;
	
	for(int i=0;i<v->len;i++) {
		r+=v->data[i];
	}
	return r;
}

double vec_avg(vec* v) {
	return vec_sum(v)/((double)v->len);
}

vec* vec_copy(vec* v) {
	vec* nv=vec_new(v->len);
	
	for(int i=0;i<v->len;i++) {
		nv->data[i]=v->data[i];
	}

	return nv;
}

void vec_free(vec* v) {
	free(v->data);
	free(v);
}

void vec_copyover(vec* dest, vec* src) {
	dest->len=src->len;
	dest->data=realloc(dest->data,src->len*sizeof(double));
	
	for(int i=0;i<src->len;i++) {
		dest->data[i]=src->data[i];
	}
}

void vec_mult(vec* v,double a) {
	for(int i=0;i<v->len;i++) {
		v->data[i]*=a;
	}
}

void vec_add(vec* a, vec* b) {
	for(int i=0;i<a->len;i++) {
		a->data[i]+=b->data[i];
	}
}

void vec_sub(vec* a, vec* b) {
	for(int i=0;i<a->len;i++) {
		a->data[i]-=b->data[i];
	}
}

void vec_clamp(vec* v,double vmin,double vmax) {
	double s=vec_sum(v);
	
	if(s>vmax) {
		vec_mult(v,vmax/s);
	}
	if(s<vmin) {
		vec_mult(v,vmin/s);
	}
}

void vec_scatter(vec* source,vec* dest,int root,MPI_Comm com) {
	MPI_Scatter(source->data,1,MPI_DOUBLE,&(dest->data[root]),1,MPI_DOUBLE,root,com);
}
