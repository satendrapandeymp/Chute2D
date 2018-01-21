#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define max(a,b) a > b ? a : b
#define min(a,b) a < b ? a : b


static double x[12000][4],v[12000][4],w[12000][4],fm[12000][4],torque[12000][4],acc[4];
static double m[12000],R[12000],mup[12000],m_inv[12000],I_m[12000],I_m_inv[12000];
static double s[12000][300][4],s_prev[12000][300][4];
static double Vel_species[21][250][4],RotVel_species[21][250][4],Temp_species[21][250][4],n_species[21][250];
static int cell[300][50][10],cell2d[300][100],mapc[650000],head[50000],llist[12000];
static int Nneb[12000],Nneb_prev[12000],neb[12000][300],neb_prev[12000][300],simflag;

static int Nbinx, Nbiny, Nbinz, Ncell, Nparts, Nboun, Npartstot, Ntmp, N_species, N_max, Nspecies[21];
static int iterin, iprint,iter_species;
static long int Ncontacts, count, binmax[21];
static double Lx,Ly,Lz,Lbinx,Lbiny,Lbinz,tcol,theta,rng,rngsq,mu,kl,kn,kt,el,deltin,lne2;
static double step, t, dt, pi, dely, tfinal, R_largest, rng_species[21], rhoratio[21], Rratio[21],muratio[21];
FILE *outf,*outf1,*outf2,*outf3;


// Including Headers
#include "./Headers/Rest.h"
#include "./Headers/ForceCal.h"
#include "./Headers/Move.h"
#include "./Headers/Makemaps.h"
#include "./Headers/Heaps.h"

void main(void)
{
	int i,count1;
	double delt1;
	struct timeval t1,t2;
	unsigned long start_sec, end_sec;
	long diff;

	delt1 = 0.1;									// Time step to take a image for properties calculation
	count = 1;									// Related to images from force calculation part (to average data)
	count1 = 1;									// Related to images from properties calculation part (to average data)
	Ncontacts=0;

	outf1 = fopen("./Output/posvel.dat","a");						// Particles position in every delt1 time step, getting 10 images in a time-step
	outf2 = fopen("./Output/stressinput.dat","a");						// Force b/w colliding particles in every delt2 time step, getting 20 images in a time-step

	Initialize();

	(void)gettimeofday(&t1, NULL);
	start_sec = t1.tv_sec;
	while (t<tfinal)
	{
		Heap();
		if(t>count1*delt1)
		{
			Profile();								// Averaging flow properties in layers

			for(i=Nboun+1;i<=Npartstot;i++)
			{
				fprintf(outf1,"%f %f %f %f %f %f %f %f %f %f %f %f %f\n",x[i][1],x[i][2],x[i][3],v[i][1],v[i][2],v[i][3],w[i][1],w[i][2],w[i][3],R[i],m[i],mup[i]);
			}
			count1++;
		}

		if(floor(t/1)>iprint)
		{
			outf=fopen("./Output/pos_t.dat","w");
			for(i=1;i<=Npartstot;i++) fprintf(outf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",x[i][1],x[i][2],x[i][3],v[i][1],v[i][2],v[i][3],w[i][1],w[i][2],w[i][3],R[i],m[i],mup[i]);
			fclose(outf);
			iprint++;
		}
	}
	fclose(outf1);
	fclose(outf2);

	(void)gettimeofday(&t2, NULL);
	end_sec = t2.tv_sec;
	diff =end_sec-start_sec;

	("Reached at the end at time =%f and took %d secs\n",t,diff);
	outf=fopen("./Output/time.dat","w");
	fprintf(outf,"Time Taken =%ld secs\n",diff);
	fclose(outf);

	outf=fopen("./Data/pos.dat","w");
	for(i=1;i<=Npartstot;i++) fprintf(outf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",x[i][1],x[i][2],x[i][3],v[i][1],v[i][2],v[i][3],w[i][1],w[i][2],w[i][3],R[i],m[i],mup[i]);
	fclose(outf);

	WriteAll();

	outf=fopen("./Output/iter_contacts.par","w");						// Outputing images in posvel and stressinput file with total number of contacts in stressinput file
	fprintf(outf,"%d %ld %ld\n",count1-1,count-1,Ncontacts);
	fclose(outf);
}
