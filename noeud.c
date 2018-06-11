#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "vec.h"
#include "random.h"

#define ITERS 200//itérations "globales" (échanges d'info)
#define LITERS 30//itérations pour la minimisation locale
#define RHO 1.0//rho=1 pour que la variable duale réduite (u) soit égale à la variable duale (le prix)

int NNODES=175;
double u;
double a;
double b;
double p0;
double pmax;
double pmin;
double pr;
double fi;
int pfxe;
vec* fis_global=NULL;
double pi=0.0;
double avg;
double lamb;
double lambda_e;
double reserve=0.0;// réserve en proportion de la plage de puissance disponible
double reserve_step=0.002;
double reserve_max=0.3;
double r;
vec** pis_g=NULL;

pthread_barrier_t barriere;
pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;

//#include "defs_15.h"//le problème originel
#include "defs_175.h"//le problème de Thomas
double prs[]={0,1,0.9,1.3,2,2,0,0,0,0,0,0,0,0,-0.2};
double ecart_type=0.1;//valeur d'écart-type : inverse de la fiabilité des agents qui subissent
double ecart_type_step=0.1;
double ecart_type_max=1.71;
unsigned int n_lambda_e=200;//nombre d'échantillons de prix générés
int n_situations=200;//nombre de situations par niveau de fiabilité
    
FILE *f_2e_marche;
FILE *f_simple;
FILE *f_reserve;
FILE *f_prevision;

//pour s'y retrouver à la lecture des résultats :
//char* names[]={"charbon","eolienne","eolienne","eolienne","eolienne","eolienne","barrage","panneau","panneau","datacenter","logement","usine","tram 1","tram 2","hopital"};


double f(double p,int n) {
	return as[n]*(p-p0s[n])*(p-p0s[n])+bs[n];
}

double clamp(double val,double mi,double ma){
	return fmin(fmax(val,mi),ma);
}

void agent_min(vec* pis,double a,double x0,double alpha,vec* yis,double xmax,double xmin) {//algo adapté du matlab de RLGL - smiley clin d'oeil
	double u=0;
	double z=0;
	double rho=1;
	double mu;
	double xavg;
	unsigned int N=pis->len;
	
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

void agent_min_anticip(vec* pis,double a,double x0,double alpha,vec* y,double xmax,double xmin,vec* pis_anticip,vec *lambda_es) {
	
	// calcul d'une valeur moyenne de p_i pour plusieurs valeurs de lambda
	
	vec_zero(pis_anticip);
	
	double pmin;
	double pmax;
	
	for(unsigned int i_le=0;i_le<n_lambda_e;i_le++) {//pour chaque valeur de lambda
		lambda_e=lambda_es->data[i_le];
		agent_min(pis,a,x0-lambda_e/(2*a),alpha,y,xmax,xmin); // on peut montrer que pour notre problème, la fonction de coût anticipée correspond au polynôme centré non pas en x0 mais en x0-(lambda/2a)
		double pi=vec_sum(pis);
		double pe=x0-pi+lambda_e/(2*a);//p^epsilon
		pmin=xmin;
		pmax=xmax;
		if((pi+pe)>xmax) {//La projection orthogonale ayant l'air d'être source d'opérations interdites, on se contente de borner bêtement pi+pe
			//pmax=pi-(pi+pe-xmax)/2;
			vec_mult(pis,xmax/(pi+pe));
		}
		if((pi+pe)<xmin) {
			//pmin=pi-(pi+pe-xmin)/2;
			vec_mult(pis,xmin/(pi+pe));
		}
		//vec_clamp(pis,pmin,pmax);//projection orthogonale sur xmin<xi*+xie<xmax
		vec_add(pis_anticip,pis);
	}
	vec_copyover(pis,pis_anticip);
	vec_mult(pis,((double)1)/((double)n_lambda_e));// on divise la somme par le nombre d'éléments -> moyenne
	vec_clamp(pis,xmin,xmax);// des fois qu'on sorte des bornes...
}

void travail_noeud(int world_rank) {
	double pmax;
	double pmin;
	double p0;
	double pr;
	double pi;
	double fi;
	double u;
	int pfxe;
	vec* pis=NULL;
	vec* pis_anticip=NULL;
	vec* qis=NULL;
	vec* yis=NULL;
	vec* uis=NULL;
	vec* lambda_es=NULL;
	

    int node;
    int i;

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
	pis_g[world_rank]=pis;
	
	lambda_es=vec_new(n_lambda_e);
	pis_anticip=vec_new(NNODES);
	for(unsigned int i_le=0;i_le<n_lambda_e;i_le++) {// génération des valeurs de lambda : randn(moyenne,ecart_type)
		lambda_es->data[i_le]=randn(-9.9,5.7);
	}
	qis=vec_new(NNODES);
	yis=vec_new(NNODES);
	uis=vec_new(NNODES);
	
	a=as[world_rank];
	b=bs[world_rank];
	p0=p0s[world_rank];
	pfxe=pfxes[world_rank];
	
	
	
	//premier marché naïf
	
	vec_zero(pis);
	vec_zero(uis);
	vec_zero(yis);
	vec_zero(qis);
		
	pi=0;
	
	if(world_rank==0) {
		printf("\nPremier marché naïf\n\n");
	}

	//premier marché : calcul de p^*_i
	for(i=0;i<ITERS;i++) {
		vec_copyover(yis,qis);
		vec_add(yis,uis);//y=q+u
		
		agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//calcul de p
	
		pthread_barrier_wait(&barriere);
		
		for(node=0;node<NNODES;node++) {
			qis->data[node]=pis_g[node]->data[world_rank];
		}
		
		pthread_barrier_wait(&barriere);

		
		vec_sub(qis,pis);
		vec_mult(qis,-0.5);//q=(p-pt)/2
		
		vec_add(uis,qis);
		vec_sub(uis,pis);//u:=u+q-p
	}
	
	
	/*if(world_rank==0) {
		printf("prix : %f\npuissance : %f\n",lambda_es->data[i_le],vec_sum(pis));
	}*/
	
	pi=vec_sum(pis);
	printf("Valeur de l'élément %d : %f (%f<p<%f)\n",world_rank,pi,pmin,pmax);
	
	if(world_rank==0){
		avg=0;
	}
	pthread_barrier_wait(&barriere);	
	pthread_mutex_lock(&mutex);
	avg+=pi;
	pthread_mutex_unlock(&mutex);
	//MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	
	//calcul de f_i(p^*_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
	fi=f(pi,world_rank);
	fis_global->data[world_rank]=fi;
	//MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	
	u=uis->data[0];
	
	if(world_rank==0) {
		printf("Valeur totale : %f\n",avg);
		printf("Variable duale : %f\n",u);
		//enregistrement : tous les f_i(p^*_i) et u*
		fprintf(f_simple,"%f",uis->data[0]);
		for(int j=0;j<NNODES;j++) {
			fprintf(f_simple,";%f",fis_global->data[j]);
		}
	}
	
	
	
	//premier marché avec réserve
	
	for(reserve=0;reserve<reserve_max;reserve+=reserve_step) {
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
		pfxe=pfxes[world_rank];
		
		
		//if(!pfxe) { // est-ce qu'on veut imposer une réserve aux agents qui subissent de toute façon ?
			r=(pmax-pmin)*reserve;
			pmax-=r;
			pmin+=r;
		//}
		
		//premier marché
		for(i=0;i<ITERS;i++) {
			vec_copyover(yis,qis);
			vec_add(yis,uis);//y=q+u
			
			agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//awiwiwi
			
			pthread_barrier_wait(&barriere);
			
			for(node=0;node<NNODES;node++) {
				qis->data[node]=pis_g[node]->data[world_rank];
			}
			
			pthread_barrier_wait(&barriere);
			
			vec_sub(qis,pis);
			vec_mult(qis,-0.5);//q=(p-pt)/2
			
			vec_add(uis,qis);
			vec_sub(uis,pis);//u:=u+q-p
		}
		
		pi=vec_sum(pis);
		printf("Valeur de l'élément %d : %f (%f<p<%f)\n",world_rank,pi,pmin,pmax);
		
		if(world_rank==0){
			avg=0;
		}
		pthread_barrier_wait(&barriere);	
		pthread_mutex_lock(&mutex);
		avg+=pi;
		pthread_mutex_unlock(&mutex);
		//MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		
		
		//calcul de f_i(p^*_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
		fi=f(pi,world_rank);
		fis_global->data[world_rank]=fi;
		//MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		
		u=uis->data[0];
		
		if(world_rank==0) {
			printf("Valeur totale : %f\n",avg);
			printf("Variable duale : %f\n",u);
			
			fprintf(f_reserve,"%f;%f",reserve,u);// enregistrement de reserve, u, et des f_i(p^*_i)
			for(int j=0;j<NNODES;j++) {
				fprintf(f_reserve,";%f",fis_global->data[j]);
			}
			fprintf(f_reserve,"\n");
		}
	}
	
	
	
	//premier marché avec prévision
	
	vec_zero(pis);
	vec_zero(uis);
	vec_zero(yis);
	vec_zero(qis);
		
	pi=0;
	
	if(world_rank==0) {
		printf("\nPremier marché prévisionnel\n\n");
	}

	//premier marché : calcul de p^*_i et p^epsilon_i anticipé
	for(i=0;i<ITERS;i++) {
		vec_copyover(yis,qis);
		vec_add(yis,uis);//y=q+u
		
		agent_min_anticip(pis,a,p0,RHO/2,yis,pmax,pmin,pis_anticip,lambda_es);//calcul de p
		
		pthread_barrier_wait(&barriere);
		
		for(node=0;node<NNODES;node++) {
			qis->data[node]=pis_g[node]->data[world_rank];
		}
		
		pthread_barrier_wait(&barriere);

		
		vec_sub(qis,pis);
		vec_mult(qis,-0.5);//q=(p-pt)/2
		
		vec_add(uis,qis);
		vec_sub(uis,pis);//u:=u+q-p
	}
	
	
	/*if(world_rank==0) {
		printf("prix : %f\npuissance : %f\n",lambda_es->data[i_le],vec_sum(pis));
	}*/
	
	pi=vec_sum(pis);
	printf("Valeur de l'élément %d : %f (%f<p<%f)\n",world_rank,pi,pmin,pmax);
	
	
	if(world_rank==0){
		avg=0;
	}
	pthread_barrier_wait(&barriere);	
	pthread_mutex_lock(&mutex);
	avg+=pi;
	pthread_mutex_unlock(&mutex);
	//MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	
	//calcul de f_i(p^*_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
	fi=f(pi,world_rank);
	fis_global->data[world_rank]=fi;
	//MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	
	u=uis->data[0];
	
	if(world_rank==0) {
		printf("Valeur totale : %f\n",avg);
		printf("Variable duale : %f\n",u);
		//enregistrement : tous les f_i(p^*_i) et u*
		fprintf(f_prevision,"%f",uis->data[0]);
		for(int j=0;j<NNODES;j++) {
			fprintf(f_prevision,";%f",fis_global->data[j]);
		}
	}
	
	
	//second marché (toujours le même)
	
	if(world_rank==0) {
		printf("\nSecond marché\n\n");
	}
	
	for(ecart_type=ecart_type;ecart_type<=ecart_type_max;ecart_type+=ecart_type_step) {//pour différentes valeurs de fiabilité
		for(int situation=0;situation<n_situations;situation++) {//on se place dans plusieurs situations, pour avoir une certaine représentativité
				
			vec_zero(pis);
			vec_zero(uis);
			vec_zero(yis);
			vec_zero(qis);
			
			
			//second marché : calcul de p^epsilon_i
			pmax=pmaxs[world_rank];
			pmin=pmins[world_rank];
			//pr=pmin+drand()*(pmax-pmin);
			pr=clamp(randn(p0,ecart_type),pmin,pmax);//on génère aléatoirement la valeur de puissance subie, en fonction de la fiabilité
			if(pfxe) {
				//pr=pmin+drand()*(pmax-pmin);
				pmin=pr;//si on subit, alors on subit
				pmax=pr;
			}
			
			for(i=0;i<ITERS;i++) {
				vec_copyover(yis,qis);
				vec_add(yis,uis);//y=q+u
				
				agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//calcul de p
				//agent_min(pis,a,p0,RHO/2,yis,pmax,pmin);//calcul de p
				
				
				pthread_barrier_wait(&barriere);
				
				for(node=0;node<NNODES;node++) {
					qis->data[node]=pis_g[node]->data[world_rank];
				}
				
				pthread_barrier_wait(&barriere);
				
				/*for(node=0;node<NNODES;node++) {
					vec_scatter(pis,qis,node,MPI_COMM_WORLD);
				}//transmission de pt (dans q)
				*/
				vec_sub(qis,pis);
				vec_mult(qis,-0.5);//q=(p-pt)/2
				
				vec_add(uis,qis);
				vec_sub(uis,pis);//u:=u+q-p
			}
			
			//pi=vec_sum(pis)+pi;
			pi=vec_sum(pis);
			printf("Valeur de l'élément %d : %f (%f<p<%f)\n",world_rank,pi,pmin,pmax);
			
			if(world_rank==0){
				avg=0;
			}
			pthread_barrier_wait(&barriere);	
			pthread_mutex_lock(&mutex);
			avg+=pi;
			pthread_mutex_unlock(&mutex);
			//MPI_Reduce(&pi,&avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			
			
			//calcul de f_i(p^*_i) et transmission au process 0, qui ne fera qu'enregister dans un fichier : pas indispensable au fonctionnement
			fi=f(pi,world_rank);
			fis_global->data[world_rank]=fi;
			//MPI_Gather(&fi,1,MPI_DOUBLE,fis_global->data,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			
			if(world_rank==0) {
				printf("Valeur totale : %f\n",avg);
				printf("Variable duale : %f\n",uis->data[0]);
				
				//enregistrement fichier : reserve, u second marché
				fprintf(f_2e_marche,"%f;%f",ecart_type,uis->data[0]);
				
				//enregistrement : tous les f_i(p^*_i+p^epsilon_i)
				for(int j=0;j<NNODES;j++) {
					fprintf(f_2e_marche,";%f",fis_global->data[j]);
				}
				fprintf(f_2e_marche,"\n");
			}
		}
	}
}

void *thread_noeud(void *arg) {
	int n=(int)arg;
	
	travail_noeud(n);
	
	pthread_exit(NULL);
}

int main(int argc, char** argv) {

	pthread_t *tid;
    int node;
        
	f_2e_marche=fopen("2e_marche.csv","w");
	f_simple=fopen("simple.csv","w");
	f_reserve=fopen("reserve.csv","w");
	f_prevision=fopen("prevision.csv","w");
	
	pthread_barrier_init(&barriere,NULL,NNODES);
	
	srand(time(NULL));
	
	tid=calloc(sizeof(pthread_t),NNODES);
	
	pis_g=malloc(NNODES*sizeof(pis_g));
	
	for(node=0;node<NNODES;node++) {
		pthread_create(&tid[node],NULL,thread_noeud,(void *)node);
	}
	
	for(node=0;node<NNODES;node++) {
		pthread_join(tid[node],NULL);
	}
	
	fclose(f_2e_marche);
	fclose(f_reserve);
	fclose(f_prevision);
	fclose(f_simple);
	
	return 0;
}
