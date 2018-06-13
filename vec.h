#include <stdlib.h>

typedef struct{
	unsigned int len;
	float * data;
} vec;


void vec_zero(vec* v) {
	for(int i=0;i<v->len;i++) {
		v->data[i]=0;
	}
}

vec* vec_new(unsigned int len) {
	vec* v=malloc(sizeof(vec));
	v->data=malloc(len*sizeof(float));
	v->len=len;
	return v;
}

vec* vec_new_fromtab(float* tab,unsigned int len) {
	vec* v=malloc(sizeof(vec));
	v->data=malloc(len*sizeof(float));
	v->len=len;
	for(int i=0;i<v->len;i++) {
		v->data[i]=tab[i];
	}
	return v;
}

float vec_sum(vec* v) {
	float r=0;
	
	for(int i=0;i<v->len;i++) {
		r+=v->data[i];
	}
	return r;
}

float vec_avg(vec* v) {
	return vec_sum(v)/((float)v->len);
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
	dest->data=realloc(dest->data,src->len*sizeof(float));
	
	for(int i=0;i<src->len;i++) {
		dest->data[i]=src->data[i];
	}
}

void vec_mult(vec* v,float a) {
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

void vec_clamp(vec* v,float vmin,float vmax) {
	float s=vec_sum(v);
	
	if(s>vmax) {
		vec_mult(v,vmax/s);
	}
	if(s<vmin) {
		vec_mult(v,vmin/s);
	}
}
