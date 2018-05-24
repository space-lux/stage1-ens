#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"

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
double avg;
double lamb;


double as[]={	0.6,1,1,1,1,1,0.01,0.9,0.9,5,1,5,6,10,200};
double bs[]={1000,10,10,10,10,10,100,11,11,20,10,8,9,9,30};
double p0s[]={	8.5,1.9,1.9,1.9,1.9,1.9,3,		1.9,1.9,-10,-7.5,	-9,	-0.12,	-0.12,	-0.2};
double pmaxs[]={10,	2,	2,	2,	2,	2,	16,		2,	2,	0,		0,		0,	0,		0,		-0.2};
double pmins[]={0,	0,	0,	0,	0,	0,	-10,	0,	0,	-15,  -10,	-10,	-0.5,		-0.5,	-0.2};
double prs[]={0,1,0.9,1.3,2,2,0,0,0,0,0,0,0,0,-0.2};
int pfxes[]={0,1,1,1,1,1,0,0,0,0,0,0,0,0,1};
double reserve=0.0;//0.1
double r;

char* names[]={"charbon","eolienne","eolienne","eolienne","eolienne","eolienne","barrage","panneau","panneau","datacenter","logement","usine","tram 1","tram 2","hopital"};


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
    int node;
    int i;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	NNODES=world_size;
	
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
	pr=prs[world_rank];
	pfxe=pfxes[world_rank];
	avg=0;
	u=1.0;
	
	
	//if(!pfxe) {
		r=(pmax-pmin)*reserve;
		pmax-=r;
		pmin+=r;
	//}
	
	//premier marché
	for(i=0;i<ITERS;i++) {
		vec_copyover(yis,qis);
		vec_add(yis,uis);//y=q+u
		
		agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//awiwiwi
		
		for(node=0;node<NNODES;node++) {
			vec_scatter(pis,qis,node,MPI_COMM_WORLD);
		}//transmission de pt (dans q)
		
		vec_sub(qis,pis);
		vec_mult(qis,-0.5);//q=(p-pt)/2
		
		vec_add(uis,qis);
		vec_sub(uis,pis);//u:=u+q-p
	}
	
	pi=vec_sum(pis);
	printf("Valeur de l'élément %d (%s) : %f (%f<p<%f)\n",world_rank,(names[world_rank]),pi,pmin,pmax);
	
	MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if(world_rank==0) {
		printf("Valeur totale : %f\n",avg);
		u=uis->data[0];//vu que c'est la même partout...
		printf("Variable duale : %f\n",u);//vu que c'est la même partout...
		printf("\nSecond marché\n\n");
	}
	
	
	vec_zero(pis);
	vec_zero(uis);
	vec_zero(yis);
	vec_zero(qis);
	
	
	//second marché
	if(pfxe) {
		pmin=pr;//si on subit, alors on subit
		pmax=pr;
	} else {
		pmax=pmaxs[world_rank];
		pmin=pmins[world_rank];
	}
	
	for(i=0;i<ITERS;i++) {
		vec_copyover(yis,qis);
		vec_add(yis,uis);//y=q+u
		
		agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//awiwiwi
		
		for(node=0;node<NNODES;node++) {
			vec_scatter(pis,qis,node,MPI_COMM_WORLD);
		}//transmission de pt (dans q)
		
		vec_sub(qis,pis);
		vec_mult(qis,-0.5);//q=(p-pt)/2
		
		vec_add(uis,qis);
		vec_sub(uis,pis);//u:=u+q-p
	}
	
	pi=vec_sum(pis);
	printf("Valeur de l'élément %d (%s) : %f (%f<p<%f)\n",world_rank,(names[world_rank]),pi,pmin,pmax);
	
	MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if(world_rank==0) {
		printf("Valeur totale : %f\n",avg);
		printf("Variable duale : %f\n",uis->data[0]);
		
		//récap
		printf("\nRécap\n\n");
		printf("réserve : %f\n",reserve);
		printf("prix premier marché : %f\n",u);
		printf("prix second marché : %f\n",uis->data[0]);
		printf("(diff relative : %f)\n",(uis->data[0]-u)/uis->data[0]);
	}
	
	
    MPI_Finalize();
	return 0;
}
