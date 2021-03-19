#include "global_variables.h"
#include "time_evolution.h"
#include "pfr_cufft1.h"
#include "misc.h"
#include "inhom_v3.h"

double  cr,nr,mur,hr;
double g2;
cufftDoubleComplex ck,nk,muk;
float randnum[L][M][N];
void GetEintC(int timestep)
{
  int g_i,g_j,g_k,gijk;
  int x_i,x_j,x_k,xijk;
  int l;
  int i;
  double n;
  //Evolving composition c(k)
  compute_mur_c(mu_cr);
  cufftrc3k(&mu_cr[0][0][0],&mu_ck[0][0][0]);
#pragma omp parallel for private(g_i,g_j,g_k,gijk,l,temp1)
#pragma acc data present(c_k,g_mod2,mu_ck)
#pragma acc parallel loop collapse(3) private(muk,g2,ck,g_i,g_j,g_k)
//#pragma acc kernels
  fgijk
	ck=c_k[g_i][g_j][g_k];
  	g2=g_mod2[g_i][g_j][g_k];
        muk=mu_ck[g_i][g_j][g_k];
  	c_k[g_i][g_j][g_k].x=(ck.x-Mc*g2*muk.x*dt)/            \
						(1+kc*Mc*g2*g2*dt);
  	c_k[g_i][g_j][g_k].y=(ck.y-Mc*g2*muk.y*dt)/            \
						(1+kc*Mc*g2*g2*dt);
  efgijk
  //--------------------------------------------------

  //Evolving eta eta(k)
  compute_mur_eta(mu_r);
  for(i=0;i<V;i++){
  cufftrc3k(&mu_r[i][0][0][0],&mu_k[i][0][0][0]); 
#pragma omp parallel for private(g_i,g_j,g_k,l)
#pragma acc data present(eta_k,g_mod2,mu_k)
#pragma acc parallel loop collapse(3) private(muk,g2,nk,g_i,g_j,g_k)
fgijk
	nk=eta_k[i][g_i][g_j][g_k];
  	g2=g_mod2[g_i][g_j][g_k];
    muk=mu_k[i][g_i][g_j][g_k];
	eta_k[i][g_i][g_j][g_k].x=(nk.x-Meta*(muk.x)*dt)/(1+Meta*keta*g2*dt);
	eta_k[i][g_i][g_j][g_k].y=(nk.y-Meta*(muk.y)*dt)/(1+Meta*keta*g2*dt);
  efgijk
  cufftcr3k(&eta_k[i][0][0][0],&eta_r[i][0][0][0]);
}
  //Converting back to real space
  cufftcr3k(&c_k[0][0][0],&c_r[0][0][0]);
//  cufftcr3k(&eta_k[0][0][0],&eta_r[0][0][0]);
  
#ifdef NOISE_ON
  if(timestep<=NRAND)
  {
//if(L!=2 && M!=2)//for 3D case
//{
  //Adding noise 
#pragma acc data create(randnum) present(eta_r)
  for(i=0;i<V;i++)
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	eta_r[i][x_i][x_j][x_k]+=Namp*(randnum[x_i][x_j][x_k]-0.5);
  efor_xijk
}
//}
/*
else if(L==2 && M!=2)//for 2D case
{
  //Adding noise
#pragma acc data create(randnum) present(eta_r)
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	if(x_i==0)
	{eta_r[x_i][x_j][x_k]+=Namp*(randnum[x_i][x_j][x_k]-0.5);
	eta_r[1][x_j][x_k]=eta_r[x_i][x_j][x_k];
	}
  efor_xijk
		}//end of data region
}
else if(L==2 && M==2)//for 1D case
{
  //Adding noise
  for_xijk
	if(x_i==0 && x_j==0)
	{eta_r[x_i][x_j][x_k]+=Namp*(unirand()-0.5);
	eta_r[0][1][x_k]=eta_r[x_i][x_j][x_k];
	eta_r[1][0][x_k]=eta_r[x_i][x_j][x_k];
	eta_r[1][1][x_k]=eta_r[x_i][x_j][x_k];
	}
  efor_xijk

}
  }
  else
  {
#pragma acc data create(randnum) present(eta_r)
		{
curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  for_xijk
	eta_r[x_i][x_j][x_k]+=Namp2*(randnum[x_i][x_j][x_k]-0.5);
  efor_xijk
}
 */ 
  }
  
#endif  

//cut off eta field
#pragma acc data present(eta_r)
  for(i=0;i<V;i++)
  {
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  fxijk
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]>1.0)?1.0:eta_r[i][x_i][x_j][x_k];
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]<0.0)?0.0:eta_r[i][x_i][x_j][x_k];
  efxijk
//modified to generate void structure
/*
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
  fxijk
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]>0.5)?1.0:eta_r[i][x_i][x_j][x_k];
  	eta_r[i][x_i][x_j][x_k]=(eta_r[i][x_i][x_j][x_k]<0.5)?0.0:eta_r[i][x_i][x_j][x_k];
  efxijk
*/
cufftrc3k(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);	 
}

#ifdef MISFIT_COUPLED_WITH_ETA
 //Updating h(eta)
// printf("misfit coupled with eta\n");
//#pragma omp parallel for private(x_i,x_j,x_k,xijk,n)
#pragma acc data present(heta_r,eta_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k,n)
for_xijk
heta_r[x_i][x_j][x_k]=0;
  for(i=0;i<V;i++)
  {
	n=eta_r[i][x_i][x_j][x_k];
	heta_r[x_i][x_j][x_k]+=n*n*(3-2*n);
//	heta_r[x_i][x_j][x_k]+=eta_r[i][x_i][x_j][x_k];
}
efor_xijk
// printf("misfit coupled with eta end\n");
  cufftrc3k(&heta_r[0][0][0],&heta_k[0][0][0]);
 
#endif

#ifdef MISFIT_COUPLED_WITH_CONC
//Updating h(c)
#pragma acc data present(c_r,hc_r)
#pragma acc parallel loop collapse(3) private(cr,x_i,x_j,x_k)
for_xijk
cr=c_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=1-cr*cr*(3-2*cr);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
cufftrc3k(&hc_r[0][0][0],&hc_k[0][0][0]);	
#endif
/*
#ifdef MISFIT_COUPLED_WITH_ETA
//Updating h(c)
#pragma acc data present(eta_r,hc_r)
#pragma acc parallel loop collapse(3) private(nr,x_i,x_j,x_k)
for_xijk
nr=eta_r[x_i][x_j][x_k];
hc_r[x_i][x_j][x_k]=nr*nr*(3-2*nr);
//hc_r[x_i][x_j][x_k]=1/3*c_r[x_i][x_j][x_k]-1.767;
efor_xijk
cufftrc3k(&hc_r[0][0][0],&hc_k[0][0][0]);	
#endif
*/
}//end of GetEintC

static inline void compute_mur_c(double mu[L][M][N])
{
int x_i,x_j,x_k,xijk;	
int g_i,g_j,g_k,gijk;
int j;
//inhom_solver();// to compute the local stress

#pragma omp parallel for private(x_i,x_j,x_k,xijk,j)
#pragma acc data present(eta_r,c_r,mu,Sig,E_oo,Sij0)
#pragma acc parallel loop collapse(3) private(nr,hr,cr,x_i,x_j,x_k,j)
for_xijk
	hr=0;
	for(j=0;j<V;j++)
	{
	nr=eta_r[j][x_i][x_j][x_k];
	hr+=nr*nr*(3-2*nr);
	}
	cr=c_r[x_i][x_j][x_k];
	//mu from chemical free energy
	  /*
		mu[x_i][x_j][x_k]=hr*A1* \
	( cr*(2+cr*(-6+cr*4)) + B1*6*cr*(cr-1) ) \
		+ (1-hr)*A2*6*cr*(1-cr);
*/
		mu[x_i][x_j][x_k]=A1*(hr* \
        (0.0312*cr*cr+4.9*cr-130.75  ) \
		+ (1-hr)*(2.9*cr-103.67));

#ifdef MISFIT_COUPLED_WITH_CONC
		voigt e_i;//elastic strain
		voigt sigma;// local stress
			real hrp=6*cr*(cr-1);
			//real hrp=1/3;
		V6_loop{sigma[mi]=Sig[xijk][mi];}
		V6_loop{
				e_i[mi]=0.0;
				V6p_loop{e_i[mi]+=Sij0[mi][mip]*sigma[mip];}
			}
		real mu_elas=0.0;
		V6_loop{
				mu_elas+=(-sigma[mi]*E_oo[mi]);
				V6p_loop{mu_elas+=(-0.5*Cij2[mi][mip]*e_i[mi]*e_i[mip]);}
		}
		mu_elas*=hrp;
		mu[x_i][x_j][x_k]+=mu_elas;
#endif

efor_xijk
}//end of compute_mu_c()

static inline void compute_mur_eta(double mu[V][L][M][N])
{
int x_i,x_j,x_k,xijk;	
int g_i,g_j,g_k,gijk;
int i;
  // inhom_solver();// to compute the local stress

//Computing mu_eta
for(i=0;i<V;i++)
{
#pragma acc data present(c_r,eta_r,mu,Sig,E_oo,Sij0)
#pragma acc parallel loop collapse(3) private(nr,cr,x_i,x_j,x_k,hr)
for_xijk
	cr=c_r[x_i][x_j][x_k];
	nr=eta_r[i][x_i][x_j][x_k];
/*
#ifdef ALLOW_APB
		mu[x_i][x_j][x_k]= \
			+6*nr*(1-fabs(nr))*(                       	   \
				A1*( cr*cr*(1+cr*(-2+cr)) + B1 )       \
				-(A2+A1*B1)*( cr*cr*(3-2*cr) )   \
				)                                  \
			+W*( nr*(2-6*fabs(nr)+4*nr*nr) ) ;
#else
*/
/*
		mu[x_i][x_j][x_k]= \
			+6*nr*(1-nr)*(                       	   \
				( cr*cr*(1+cr*(-2+cr)) + B1 )       \
				-(A2+A1*B1)*( cr*cr*(3-2*cr) )   \
				)                                  \
			+W*( nr*(2+nr*(-6+4*nr)) ) ;
*/
		mu[i][x_i][x_j][x_k]= \
			+6*nr*(1-nr)*(                       	   \
			A1*( 0.0104*cr*cr*cr+cr*cr-27.08*cr+136.67      \
				 )   \
				)                                  \
			+W*( nr*(2+nr*(-6+4*nr)) ) ;
#ifdef MISFIT_COUPLED_WITH_ETA
		voigt e_i;//elastic strain
		voigt sigma;// local stress
			real hrp=6*nr*(1-nr);
			//real hrp=1/3;
		V6_loop{sigma[mi]=Sig[xijk][mi];}
		V6_loop{
				e_i[mi]=0.0;
				V6p_loop{e_i[mi]+=Sij0[mi][mip]*sigma[mip];}
			}
		real mu_elas=0.0;
		V6_loop{
				mu_elas+=(-sigma[mi]*E_oo[mi]);
				V6p_loop{mu_elas+=(-0.5*Cij2[mi][mip]*e_i[mi]*e_i[mip]);}
		}
		mu_elas*=hrp;
		mu[i][x_i][x_j][x_k]+=mu_elas;
#endif
efor_xijk
}
}//end of compute_mu_eta()

