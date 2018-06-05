#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vec.h"
#include "random.h"

#define ITERS 150//itérations "globales" (échanges d'info)
#define LITERS 20//itérations pour la minimisation locale
#define RHO 1.0//rho=1 pour que la variable duale réduite (u) soit égale à la variable duale (le prix)

int NNODES;
double u;
double a;
double b;
double p0;
double pmax;
double pmin;
double pr;
double fi;
int pfxe;
double pi=0.0;
vec* pis=NULL;
vec* pis_anticip=NULL;
vec* qis=NULL;
vec* yis=NULL;
vec* uis=NULL;
vec* lambda_es=NULL;
vec* fis_global=NULL;
double avg;
double lamb;
double lambda_e;


double as[]={	0.6,1,1,1,1,1,0.01,0.9,0.9,5,1,5,6,10,200};
double bs[]={1000,10,10,10,10,10,100,11,11,20,10,8,9,9,30};
double p0s[]={	8.5,1.9,1.9,1.9,1.9,1.9,3,		1.9,1.9,-10,-7.5,	-9,	-0.12,	-0.12,	-0.2};
double pmaxs[]={10,	2,	2,	2,	2,	2,	16,		2,	2,	0,		0,		0,	0,		0,		-0.2};
double pmins[]={0,	0,	0,	0,	0,	0,	-10,	0,	0,	-15,  -10,	-10,	-0.5,		-0.5,	-0.2};
double prs[]={0,1,0.9,1.3,2,2,0,0,0,0,0,0,0,0,-0.2};
int pfxes[]={0,1,1,1,1,1,0,0,0,0,0,0,0,0,1};
double ecart_type=0.1;// écart-type (inverse de la fiabilité) subi par les agents pas fiables
double ecart_type_step=0.1;
double ecart_type_max=1.61;
unsigned int n_lambda_e=100;
int n_situations=50;//nombre de situations par niveau de fiabilité

//pour s'y retrouver à la lecture des résultats :
char* names[]={"charbon","eolienne","eolienne","eolienne","eolienne","eolienne","barrage","panneau","panneau","datacenter","logement","usine","tram 1","tram 2","hopital"};


double f(double p) {
	return a*(p-p0)*(p-p0)+b;
}

double clamp(double val,double mi,double ma){
	return fmin(fmax(val,mi),ma);
}

void agent_min(vec* pis,double a,double x0,double alpha,vec* y,double xmax,double xmin) {//algo adapté du matlab de RLGL - smiley clin d'oeil
	double u=0;
	double z=0;
	double rho=1;
	double mu;
	double xavg;
	unsigned int N=yis->len;
	
	//résolution du problème local de partage par ADMM
	
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
	
	vec_clamp(pis,xmin,xmax); // on borne p_i = somme(p_ij)
}

void agent_min_anticip(vec* pis,double a,double x0,double alpha,vec* y,double xmax,double xmin) {
	
	// calculd'une valeur moyenne de p_i pour plusieurs valeurs de 
	
	vec_zero(pis_anticip);
	
	for(unsigned int i_le=0;i_le<n_lambda_e;i_le++) {//pour chaque valeur de lambda
		lambda_e=lambda_es->data[i_le];
		agent_min(pis,a,x0-lambda_e/(2*a),alpha,y,xmax,xmin); // on peut montrer que pour notre problème, la fonction de coût anticipée correspond au polynôme centré non pas en x0 mais en x0-(lambda/2a)
		double pi=vec_sum(pis);
		double pe=x0-pi+lambda_e/(2*a);//p^epsilon
		double pmin=xmin;
		double pmax=xmax;
		if((pi+pe)>xmax) {//on adapte les bornes de p_i pour p^*_i
			pmax=pi-(pi+pe-xmax)/2;
		}
		if((pi+pe)<xmin) {
			pmin=pi-(pi+pe-xmin)/2;
		}
		vec_clamp(pis,pmin,pmax);//projection orthogonale sur xmin<xi*+xie<xmax
		vec_add(pis_anticip,pis);
	}
	vec_copyover(pis,pis_anticip);
	vec_mult(pis,((double)1)/((double)n_lambda_e));// on divise la somme par le nombre d'éléments -> moyenne
	vec_clamp(pis,xmin,xmax);// des fois qu'on sorte des bornes...
}



int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

	
    // Get the number of processes
    int world_size;
    int node;
    int i;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	NNODES=world_size;//étape franchement inutile, vestige d'un ancien programme, mais l'enlever obligerait à changer pas mal de notations
	
	
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	srand(time(NULL)+world_rank);
    
	FILE *prix;
	FILE *cout;
	if(world_rank==0) {
		prix=fopen("price.csv","w");
		cout=fopen("cout.csv","w");
	}

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    //printf("Bongeourre depuis %s, numero %d sur %d processus\n",processor_name, world_rank, world_size);

    // Finalize the MPI environment.
	
	if(world_rank==0) {
		printf("ADMM rules\n");// l'ADMM pèse
	}
	
	pmax=pmaxs[world_rank];
	pmin=pmins[world_rank];
	fis_global=vec_new(NNODES);
	
	pis=vec_new(NNODES);
	lambda_es=vec_new(n_lambda_e);
	pis_anticip=vec_new(NNODES);
	for(unsigned int i_le=0;i_le<n_lambda_e;i_le++) {// génération des valeurs de lambda : randn(moyenne,ecart_type)
		lambda_es->data[i_le]=clamp(randn(0.125,0.02),pmin,pmax);
	}
	qis=vec_new(NNODES);
	yis=vec_new(NNODES);
	uis=vec_new(NNODES);
	
	a=as[world_rank];
	b=bs[world_rank];
	p0=p0s[world_rank];
	pfxe=pfxes[world_rank];
	
	
	for(ecart_type=ecart_type;ecart_type<=ecart_type_max;ecart_type+=ecart_type_step) {//pour différentes valeurs de fiabilité
		for(int situation=0;situation<n_situations;situation++) {//on se place dans plusieurs situations, pour avoir une certaine représentativité
			vec_zero(pis);
			vec_zero(uis);
			vec_zero(yis);
			vec_zero(qis);
				
			pi=0;
			avg=0;
			u=1.0;
			//pr=pmin+drand()*(pmax-pmin);
			pr=clamp(randn(p0,ecart_type),pmin,pmax);//on génère aléatoirement la valeur de puissance subie, en fonction de la fiabilité
			

			//premier marché : calcul de p^*_i et p^epsilon_i anticipé
			for(i=0;i<ITERS;i++) {
				vec_copyover(yis,qis);
				vec_add(yis,uis);//y=q+u
				
				agent_min_anticip(pis,a,p0,RHO/2,yis,pmax,pmin);//calcul de p
				
				for(node=0;node<NNODES;node++) {
					vec_scatter(pis,qis,node,MPI_COMM_WORLD);
				}//transmission de pt (dans q)
				
				vec_sub(qis,pis);
				vec_mult(qis,-0.5);//q=(p-pt)/2
				
				vec_add(uis,qis);
				vec_sub(uis,pis);//u:=u+q-p
			}
			
			
			/*if(world_rank==0) {
				printf("prix : %f\npuissance : %f\n",lambda_es->data[i_le],vec_sum(pis));
			}*/
			
			pi=vec_sum(pis);
			printf("Valeur de l'élément %d (%s) : %f (%f<p<%f)\n",world_rank,(names[world_rank]),pi,pmin,pmax);
			
			MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			
			
			//calcul de f_i(p^*_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
			fi=f(pi);
			MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			
			
			u=uis->data[0];
			
			if(world_rank==0) {
				printf("Valeur totale : %f\n",avg);
				printf("Variable duale : %f\n",u);
				//enregistrement : tous les f_i(p^*_i)
				fprintf(cout,"%f",ecart_type);
				for(int j=0;j<NNODES;j++) {
					fprintf(cout,";%f",fis_global->data[j]);
				}
			}
			
				
			vec_zero(pis);
			vec_zero(uis);
			vec_zero(yis);
			vec_zero(qis);
			
			if(world_rank==0) {
				printf("\nSecond marché\n\n");
			}
			
			//second marché : calcul de p^epsilon_i
			pmax=pmaxs[world_rank];
			pmin=pmins[world_rank];
			if(pfxe) {
				//pr=pmin+drand()*(pmax-pmin);
				pmin=pr;//si on subit, alors on subit
				pmax=pr;
			}
			
			for(i=0;i<ITERS;i++) {
				vec_copyover(yis,qis);
				vec_add(yis,uis);//y=q+u
				
				agent_min(pis,a,p0-pi,RHO/2,yis,pmax-pi,pmin-pi);//calcul de p
				//agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//calcul de p
				
				for(node=0;node<NNODES;node++) {
					vec_scatter(pis,qis,node,MPI_COMM_WORLD);
				}//transmission de pt (dans q)
				
				vec_sub(qis,pis);
				vec_mult(qis,-0.5);//q=(p-pt)/2
				
				vec_add(uis,qis);
				vec_sub(uis,pis);//u:=u+q-p
			}
			
			pi=vec_sum(pis)+pi;
			//pi=vec_sum(pis);
			printf("Valeur de l'élément %d (%s) : %f (%f<p<%f)\n",world_rank,(names[world_rank]),pi,pmin,pmax);
			
			MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			
			//calcul de f_i(p^*_i+p^epsilon_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
			fi=f(pi);
			MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			
			if(world_rank==0) {
				printf("Valeur totale : %f\n",avg);
				printf("Variable duale : %f\n",uis->data[0]);
				
				//enregistrement fichier : reserve, u premier marché, u second marché
				fprintf(prix,"%f;%f;%f\n",ecart_type,u,uis->data[0]);
				
				//enregistrement : tous les f_i(p^*_i+p^epsilon_i)
				for(int j=0;j<NNODES;j++) {
					fprintf(cout,";%f",fis_global->data[j]);
				}
				fprintf(cout,"\n");
				
				//récap
				printf("\nRécap\n\n");
				printf("écart-type (1/fiabilité) : %f\n",ecart_type);
				printf("prix premier marché : %f\n",u);
				printf("prix second marché : %f\n",uis->data[0]);
			}
		}
	}
	
	if(world_rank==0) {
		fclose(prix);
		fclose(cout);
	}
	
    MPI_Finalize();
	return 0;
}
