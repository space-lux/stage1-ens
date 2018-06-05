#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vec.h"
#include "random.h"

// uniquement pour avoir un "passé" du prix, pour faire

#define ITERS 150//itérations "globales" (échanges d'info)
#define LITERS 20//itérations pour la minimisation locale
#define RHO 1.0

int NNODES;
double u;
double a;
double b;
double p0;
double pmax;
double pmin;
double pr;
int pfxe;
double pi=0.0;
vec* pis=NULL;
vec* qis=NULL;
vec* yis=NULL;
vec* uis=NULL;
vec* us=NULL;
double avg;
double lamb;


#include "defs_300.h"
double prs[]={0,1,0.9,1.3,2,2,0,0,0,0,0,0,0,0,-0.2};
double reserve=0.0;//0.1
double r;
int n_situations=10;

//char* names[]={"charbon","eolienne","eolienne","eolienne","eolienne","eolienne","barrage","panneau","panneau","datacenter","logement","usine","tram 1","tram 2","hopital"};


double clamp(double val,double mi,double ma){
	return fmin(fmax(val,mi),ma);
}

void agent_min(vec* pis,double a,double x0,double alpha,vec* y,double xmax,double xmin) {
	double u=0;
	double z=0;
	double rho=1;
	double mu;
	double xavg;
	unsigned int N=yis->len;
	
	
	for(int iter=0;iter<LITERS;iter++) {
		xavg=vec_avg(pis);
		for(unsigned int i=0;i<N;i++) {
			mu=pis->data[i]-xavg+z-u;
			pis->data[i]=(2*alpha*yis->data[i]+ mu*rho)/(2*alpha+rho);
		}
		xavg=vec_avg(pis);
		mu=u+xavg;
		
		z=clamp((2*a*x0+ mu*rho)/(2*a*((double)N)+rho),xmin,xmax);
		u=u+xavg-z;
	}
	
	vec_clamp(pis,xmin,xmax);
}



int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    //int node;
    int i;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	NNODES=world_size;
	
	
	/*FILE *prix;
	if(world_rank==0) {
		prix=fopen("price.csv","w");
	}*/
	
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	srand(time(NULL)+world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    //printf("Bongeourre depuis %s, numero %d sur %d processus\n",processor_name, world_rank, world_size);

    // Finalize the MPI environment.
	
	if(world_rank==0) {
		printf("ADMM rules\n");
	}
	
	pis=vec_new(NNODES);
	qis=vec_new(NNODES);
	yis=vec_new(NNODES);
	uis=vec_new(NNODES);
	us=vec_new(n_situations);
	vec_zero(pis);
	vec_zero(uis);
	vec_zero(yis);
	vec_zero(qis);
	vec_zero(us);
		
	
	pfxe=pfxes[world_rank];
	
	for(int sit=0;sit<n_situations;sit++) {
		vec_zero(pis);
		vec_zero(uis);
		vec_zero(yis);
		vec_zero(qis);
		
		a=as[world_rank];
		b=bs[world_rank];
		p0=p0s[world_rank];
		pmax=pmaxs[world_rank];
		pmin=pmins[world_rank];
		
		pi=0;
		avg=0;
		u=1.0;
		
	
		
		//second marché
		
		pr=pmin+drand()*(pmax-pmin);// un peu d'aléatoire pour avoir queleque chose de "représentatif"
		//pr=clamp(randn(p0,1.5),pmin,pmax);
		if(pfxe) {
			pmin=pr;//si on subit, alors on subit
			pmax=pr;
		}
		
		for(i=0;i<ITERS;i++) {
			vec_copyover(yis,qis);
			vec_add(yis,uis);//y=q+u
			
			agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//awiwiwi
			
			for(int node=0;node<NNODES;node++) {
				vec_scatter(pis,qis,node,MPI_COMM_WORLD);//y'a un bug
			}//transmission de pt (dans q)
			
			vec_sub(qis,pis);
			vec_mult(qis,-0.5);//q=(p-pt)/2
			
			vec_add(uis,qis);
			vec_sub(uis,pis);//u:=u+q-p
		}
		
		pi=vec_sum(pis);
		printf("Valeur de l'élément %d : %f (%f<p<%f)\n",world_rank,pi,pmin,pmax);
		
		MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(world_rank==0) {
			printf("Valeur totale : %f\n",avg);
			printf("Variable duale : %f\n",uis->data[0]);
			
			us->data[sit]=uis->data[0];
			
			//fprintf(prix,"%f\n",uis->data[0]);
			/*
			//récap
			printf("\nRécap\n\n");
			printf("réserve : %f\n",reserve);
			printf("prix premier marché : %f\n",u);
			printf("prix second marché : %f\n",uis->data[0]);
			printf("(diff relative : %f)\n",(uis->data[0]-u)/uis->data[0]);*/
		}
	}
	
	if(world_rank==0) {
		double u_avg=vec_avg(us);
		double u2_avg=0;
		
		for(int sit=0;sit<n_situations;sit++) {
			u2_avg+=us->data[sit]*us->data[sit];
		}
		u2_avg/=(double)n_situations;
		
		double variance=u2_avg-u_avg*u_avg;
	
	
		printf("Moyenne : %f\nVariance : %f\n",u_avg,variance);
		//fclose(prix);
	}
	
	vec_free(uis);
	vec_free(qis);
	vec_free(yis);
	vec_free(pis);
	vec_free(us);
	
    MPI_Finalize();
	return 0;
}
