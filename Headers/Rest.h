#ifndef REST_H_INCLUDED
#define REST_H_INCLUDED


void SortNeb()
{
	int i,j,k,n,length;
	int A[300];
	double swapped;

	for (i=1;i<=Npartstot;++i)
	{
		for (n=1;n<=Nneb[i];++n) A[n]=neb[i][n];
		swapped=1;
		length=Nneb[i]-1;
		if(length>0)
		{
			while (swapped==1)
			{
				swapped=0;
				for(j=1;j<=length;++j)
				{
					if(A[j]>A[j+1])
					{
						k=A[j];
						A[j]=A[j+1];
						A[j+1]=k;
						swapped=1;
					}
				}
			}
		}

		for (n=1;n<=Nneb[i];++n) neb[i][n]=A[n];					// Arranging neighbours in ascending order of their ID
/*		if(Nneb[i]>100)
		{
			printf("Nneb=%d for i=%d\n",Nneb[i],i);
		}
*/	}
}


void LinkList(void)
{
	int i, icell;

	if(simflag==3)
	{
		for (i=1;i<=Ncell;++i) head[i]=0;
		for (i=1;i<=Npartstot;++i)
		{
			icell=1+floor(x[i][3]/Lbinz)+floor(x[i][1]/Lbinx)*Nbinz+floor(x[i][2]/Lbiny)*Nbinx*Nbinz;	// Getting cell particle associated to
			llist[i]=head[icell];										// llist of a particle - contains number of last particles in that cell
			head[icell]=i;											// Last particle in the cell
		}
	}

	if(simflag==2)
	{
		for (i=1;i<=Ncell;++i) head[i]=0;
		for (i=1;i<=Npartstot;++i)
		{
			icell=1+floor(x[i][1]/Lbinx)+floor(x[i][2]/Lbiny)*Nbinx;
			llist[i]=head[icell];
			head[icell]=i;
		}
	}
}


void AdjustSP()
{
	int i,n,p;
	int j_curr,j_prev;

	for (i=1;i<=Npartstot;++i)
	{
		for(n=1;n<=Nneb[i];++n)
		{
			j_prev=1;
			j_curr=1;
			while(j_curr<=Nneb[i] && j_prev<=Nneb_prev[i])
			{
				if (neb_prev[i][j_prev]<neb[i][j_curr]) 	// CASE 1: Current value is greater than previous value
				{
					j_prev++;
				}
				else if (neb_prev[i][j_prev]>neb[i][j_curr])	// CASE 2: Previous value is greater than current value
				{
					for (p=1; p<=3; p++) s[i][j_curr][p]=0;
					j_curr++;
				}
				else 						// CASE 3: Both values are equal
				{
					for (p=1; p<=3; p++) s[i][j_curr][p]=s_prev[i][j_prev][p];
					j_prev++;
					j_curr++;
				}
			}

			if (j_curr<=Nneb[i]) 					// Some elements in current are remaining
			{
				for(; j_curr<=Nneb[i]; j_curr++)
				for (p=1; p<=3; p++) s[i][j_curr][p]=0;
			}
		}
	}
}

void NebList(void)
{
	int i,j,icell,jcell0,jcell,n,ncheck;
	double dsq,dx1,dx2,dx3;
	double range;

	for (i=1;i<=Npartstot;++i)
	{
		for (j=1;j<=Nneb[i];++j) neb[i][j]=0;
		Nneb[i]=0;
	}
	for (icell=1;icell<=Ncell;++icell)
	{
		i=head[icell];									// Starting from last particle in the 'icell'
		while (i>0)
		{
	        	j=llist[i];								// Ptcle in its list
			while(j>0)
			{
				range = rng*(R[i]+R[j]);					// Range varying with particle size
				rngsq = range*range;
				dx1=x[i][1]-x[j][1];
				dx2=x[i][2]-x[j][2];
				dx3=x[i][3]-x[j][3];

				if(dx1>Lx-range) dx1=dx1-Lx;					// Getting distance b/w two particles, i.e. i & j
				else if (dx1<-Lx+range) dx1=dx1+Lx;
				if(dx3>Lz-range) dx3=dx3-Lz;
				else if (dx3<-Lz+range) dx3=dx3+Lz;
				dsq=dx1*dx1+dx2*dx2+dx3*dx3;

				if(dsq<rngsq)							// if within the range, adding in the neighbourlist of particle with smaller id number
				{
					if (i>j)
					{
						Nneb[i]=Nneb[i]+1;				// Number of neighbours of i
						neb[i][Nneb[i]]=j;				// Id of neighbour
	  				}
	  				else
					{
						Nneb[j]=Nneb[j]+1;
						neb[j][Nneb[j]]=i;
					}
				}
				j=llist[j];							// Next particle to check with i, that is in list of j
			}

			ncheck=0.5*(pow(3,simflag)-1);

			jcell0=ncheck*(icell-1);						// checking with particles in neighbouring cell
			for (n=1;n<=ncheck;++n)
			{
				jcell=mapc[jcell0+n];
				if (jcell>0) j=head[jcell];
				else j=0;

				while (j>0)
				{
					range = rng*(R[i]+R[j]);
					rngsq = range*range;
					dx1=x[i][1]-x[j][1];
					dx2=x[i][2]-x[j][2];
					dx3=x[i][3]-x[j][3];
					if(dx1>Lx-range) dx1=dx1-Lx;
					else if (dx1<-Lx+range) dx1=dx1+Lx;
					if(dx3>Lz-range) dx3=dx3-Lz;
					else if (dx3<-Lz+range) dx3=dx3+Lz;
					dsq=dx1*dx1+dx2*dx2+dx3*dx3;
					if(dsq<rngsq)
					{
						if (i>j)
						{
							Nneb[i]=Nneb[i]+1;
							neb[i][Nneb[i]]=j;
						}
						else
						{
							Nneb[j]=Nneb[j]+1;
							neb[j][Nneb[j]]=i;
						}
					}
					j=llist[j];
				}
			}
			i=llist[i];								// Next ptcle to check for, one in list of i
		}
	}
	SortNeb();
}


void Boundary(void)
{
	int i;
	for (i=Nboun+1;i<=Npartstot;++i)
	{
		if(x[i][1]>Lx) x[i][1]=x[i][1]-Lx;					// As PBC applied, if ptcle crossing a boundary, reentering it from other side
		else if(x[i][1]<0) x[i][1]=x[i][1]+Lx;

		if(x[i][2]<0)
		{
			printf("Ptcle %d penetrated the base\n",i);
			exit(1);
		}

		if(x[i][2]>(Ly-R[i]) && v[i][2]>0)					// Rebounding ptcle back in box without any loss of energy
		{
			v[i][2]=-v[i][2];
			printf("Ptcle %i reached to the top at time %f\n",i,t);
		}

		if(x[i][3]>Lz) x[i][3]=x[i][3]-Lz;
		else if(x[i][3]<0) x[i][3]=x[i][3]+Lz;
	}
}


void Profile(void)
{
	int i,j,k,p,binnum;
	double yavg[N_species+1], vavg[N_species+1], ydev[N_species+1];

	k = Nboun;
	for(i=1;i<=N_species;++i)
	{
		for(j=k+1;j<=k+Nspecies[i];++j)							// COM of a species
		{
			yavg[i] += x[j][2];
			vavg[i] += v[j][2];
		}
		yavg[i] = yavg[i]/Nspecies[i];
		vavg[i] = vavg[i]/Nspecies[i];

		for(j=k+1;j<=k+Nspecies[i];++j)
		{
			ydev[i] += (x[j][2]-yavg[i])*(x[j][2]-yavg[i]);
		}
		ydev[i] = ydev[i]/Nspecies[i];

		k = k + Nspecies[i];
	}

	outf=fopen("./Output/Ptcle_com.txt","a");
	fprintf(outf,"%f",t);
	for(i=1;i<=N_species;++i)
	{
		fprintf(outf," %f %f %f",yavg[i],vavg[i],ydev[i]);
	}
	fprintf(outf,"\n");
	fclose(outf);

	k = Nboun;
	for(i=1;i<=N_species;++i)								// Averaging flow properties in different bins according to height
	{
		for(j=k+1;j<=k+Nspecies[i];++j)
		{
			binnum=floor(x[j][2]/dely)+1;
			for(p=1;p<=3;++p)
			{
				Vel_species[i][binnum][p]=Vel_species[i][binnum][p]+v[j][p];
				RotVel_species[i][binnum][p]=RotVel_species[i][binnum][p]+w[j][p];
				Temp_species[i][binnum][p]=Temp_species[i][binnum][p]+v[j][p]*v[j][p];
			}
			n_species[i][binnum]=n_species[i][binnum]+1;

			if(binnum>binmax[i])
				binmax[i] = binnum;
		}
		k = k + Nspecies[i];
	}
	iter_species = iter_species + 1;
}


void WriteAll(void)
{
	int i,j,p;
	double Vint;
	char fname[50];

	Vint=Lx*Lz*dely;

	for(j=1;j<=N_species;++j)
	{
		for(i=1;i<=binmax[j];i++)
		{
			if((n_species[j][i]>0))
			{
				for(p=1;p<=3;++p)
				{
					Vel_species[j][i][p]=Vel_species[j][i][p]/n_species[j][i];
					RotVel_species[j][i][p]=RotVel_species[j][i][p]/n_species[j][i];
					Temp_species[j][i][p]=Temp_species[j][i][p]/n_species[j][i]-Vel_species[j][i][p]*Vel_species[j][i][p];
				}
				n_species[j][i]=n_species[j][i]/Vint/(iter_species-1);
			}
		}
	}

	for(j=1;j<=N_species;++j)
	{
		sprintf(fname,"./Output/profile_%d.txt",j);
		outf2 = fopen(fname,"w");
		for(i=1;i<=binmax[j];i++)
		{
			if(n_species[j][i]>0)
				fprintf(outf2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",(i-0.5)*dely,Vel_species[j][i][1],Vel_species[j][i][2],Vel_species[j][i][3],RotVel_species[j][i][1],RotVel_species[j][i][2],RotVel_species[j][i][3],Temp_species[j][i][1],Temp_species[j][i][2],Temp_species[j][i][3],n_species[j][i]);
		}
		fclose(outf2);
	}
}


void Initialize(void)
{
	FILE *inp;
	int i,p,k;
	double Lmin,gamma,rng_large,Rbasis,Mbasis,mubasis;

	iprint=iter_species=0;

	inp = fopen("./Data/ch3d.par","r");
	fscanf(inp,"%i%*s\n", &simflag);
	fscanf(inp,"%i%lf%lf%lf%*s%*s%*s%*s\n", &N_species, &Rbasis, &Mbasis, &mubasis);
	fscanf(inp,"%lf%lf%lf%lf%lf%lf%lf%i%i%lf%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s\n", &el, &mu, &kl, &theta, &Lx, &Ly, &Lz, &Nboun, &Nparts, &tfinal);
	fscanf(inp,"%i%lf%lf%*s%*s%*s\n", &Ntmp, &rng, &dely);

	fscanf(inp,"%*s%*s%*s%*s\n");
	for(i=1;i<=N_species;i++)
	{
		fscanf(inp,"%lf%lf%lf%i", &rhoratio[i], &Rratio[i], &muratio[i], &Nspecies[i]);
	}
	fclose(inp);

	Npartstot=Nparts+Nboun;									// Total Particles in the system
	inp = fopen("./Data/pos.dat","r");
	for (i=1;i<=Npartstot;++i)
		fscanf(inp, "%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg%lg", &x[i][1], &x[i][2], &x[i][3],&v[i][1], &v[i][2], &v[i][3],&w[i][1], &w[i][2], &w[i][3], &R[i],&m[i],&mup[i]);
	fclose(inp);

	for(i=1;i<=Npartstot;++i)								// To get moment of inertia of all particles
	{
		m_inv[i]=1/m[i];
		if(simflag==2) I_m[i]=(0.5*m[i]*R[i]*R[i]);
		else if(simflag==3) I_m[i]=(0.4*m[i]*R[i]*R[i]);
		I_m_inv[i]=1/I_m[i];
	}

	N_max = 0;										// To allocate memory to variables in forcecalc function
	R_largest = R[1];									// For getting the maximum range required in the code
	printf("rng=%f ",rng);
	for(i=1;i<=N_species;i++)
	{
		if(Nspecies[i]>N_max)
			N_max = Nspecies[i];

		if((Rratio[i]*Rbasis)>R_largest)
			R_largest = (Rratio[i]*Rbasis);
	}

	rng_large = rng;
	for(i=1;i<=N_species;i++)
	{
		rng_species[i] = rng*((Rratio[i]*Rbasis)+R_largest);
		printf("rng-%d=%f ",i,rng_species[i]);

		if(rng_species[i]>rng_large)
			rng_large = rng_species[i];
	}
	printf("\n");

	pi=atan(1.0)*4;										// Value of pi
	lne2=log(el)*log(el);
	gamma=sqrt((4*kl*lne2)/(0.5*(pi*pi+lne2)));						// gamma corresponding to monodispersed case with m = 1
	tcol=pi/sqrt(kl/0.5-gamma*gamma/4);							// Collision time
	step=tcol/50;										// Step for calculation of forces and updating position and velocity

	Lmin=rng*rng_large;									// Getting sizes of bins

	Nbinx=floor(Lx/Lmin);
	Lbinx=Lx/Nbinx;

	Nbiny=floor(Ly/Lmin);
	Lbiny=Ly/Nbiny;

	if(simflag==2)	Nbinz=1;
	else if(simflag==3) Nbinz=floor(Lz/Lmin);
	Lbinz=Lz/Nbinz;

	kn=kl;											// Normal and tangential spring constant
	kt=2*kn/7;
	Ncell=Nbinx*Nbiny*Nbinz;

	acc[1] = sin(theta*pi/180.0);								// Acceleration due to gravity
	acc[2] = -cos(theta*pi/180.0);
	acc[3] = 0.0;

	for (i=1;i<=Npartstot;++i)
	{
		for (p=1;p<=3;++p)
		{
			for (k=1;k<=Nneb[i];++k) s[i][k][p]=0;					// Initializing slip of particle
		}
	}

	iterin=1;
	deltin=0.1;
	MakeMap();
	LinkList();
	NebList();
}


#endif // TEST_H_INCLUDED



