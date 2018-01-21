#ifndef HEAPS_H_INCLUDED
#define HEAPS_H_INCLUDED


void Heap(void)
{
	int i,p,n,itmp;
	double ke,re;

	ke=re=0.0;

	if(simflag==2)
	{
		for (itmp=1;itmp<=Ntmp;++itmp) 							// Ntmp times, updating positions of particles before updating Neblist
		{
			dt = step;
			Move2D();
			t = t + dt;
			for (i=Nboun+1;i<=Npartstot;++i)
			{
			    	for (p=1;p<=2;++p)
				{
					ke=ke+0.5*m[i]*v[i][p]*v[i][p];				// kinetic energy of system
				}
				re=re+0.5*I_m[i]*w[i][3]*w[i][3];				// rotational kinetic energy of system
			}
		}
	}
	else if(simflag==3)
	{
		for(itmp=1;itmp<=Ntmp;++itmp)
		{
			dt = step;
			Move3D();
			t = t + dt;
			for (i=Nboun+1;i<=Npartstot;++i)
			{
			    	for(p=1;p<=3;++p)
				{
					ke=ke+0.5*m[i]*v[i][p]*v[i][p];
					re=re+0.5*0.4*m[i]*R[i]*R[i]*w[i][p]*w[i][p];	//re=1/2*Iw^2
				}
			}
		}
	}

  	ke=ke/Nparts/Ntmp;									// Average kinetic energy per particle
	re=re/Nparts/Ntmp;
	if(t>iterin*deltin)
	{
		iterin=iterin+1;
		printf("time %lg dt=%f & ke=%f & re=%f\n",t,dt,ke,re);
		outf=fopen("./Output/KEREnv.txt","a");
		fprintf(outf,"%f %f\n",ke,re);
		fclose(outf);
	}

	Boundary();
	LinkList();

	for(i=1;i<=Npartstot;++i)
	{
		Nneb_prev[i]=Nneb[i];
		for(n=1;n<=Nneb[i];++n)
		{
			neb_prev[i][n]=neb[i][n];
			for(p=1;p<=3;++p) s_prev[i][n][p]=s[i][n][p];
		}
	}

	NebList();
	AdjustSP();
}



#endif // TEST_H_INCLUDED
