#ifndef FORCECAL_H_INCLUDED
#define FORCECAL_H_INCLUDED


void ForceCalc(double xp[][4],double vp[][4],double fmp[][4],double wp[][4], double sp[][300][4], double torquep[][4])
{
	int i,j,k,l,p,n,slip;
	double vt[4],vn[4],tij[4],dxij[4],xcontact[4],dvij[4],mij[4],vmid[4];
	double d,al,dal,st,st_prev,sdotr,vtang,NForce,TForce,meff,gamn,gamt,delt2;
	double dinv, stinv, ktinv, st_previnv,muij;
	double tijdotvt, tijdotvt_cal, vtanginv;
	double fnavg[N_species+1], ftavg[N_species+1], fn[N_species+1][N_max+1], ft[N_species+1][N_max+1];

	delt2 = 0.05;

	k = Nboun;
	for(i=1;i<=N_species;i++)
	{
		for(j=k+1;j<=k+Nspecies[i];++j)
		{
			fn[i][j-k] = ft[i][j-k] = 0;
			for (p=1;p<=3;++p)
			{
				fmp[j][p]=m[j]*acc[p];
				torquep[j][p]=0;
			}
		}
		k = k + Nspecies[i];
	}

	ktinv=1/kt;

	l = Nboun;
	for(k=1;k<=N_species;++k)
	{
		for (i=l+1;i<=l+Nspecies[k];++i)
		{
			for(n=1;n<=Nneb[i];++n)							// For neighbours of a ptcle checking which one are in contact
			{

				j=neb[i][n];
				dxij[1]=xp[i][1]-xp[j][1];
				if(dxij[1]>Lx-rng*(R[i]+R[j])) dxij[1]=dxij[1]-Lx;
				else if (dxij[1]<-Lx+rng*(R[i]+R[j])) dxij[1]=dxij[1]+Lx;
				dxij[2]=xp[i][2]-xp[j][2];
				dxij[3]=xp[i][3]-xp[j][3];
				if(dxij[3]>Lz-rng*(R[i]+R[j])) dxij[3]=dxij[3]-Lz;
				else if (dxij[3]<-Lz+rng*(R[i]+R[j])) dxij[3]=dxij[3]+Lz;
				d=sqrt(dxij[1]*dxij[1]+dxij[2]*dxij[2]+dxij[3]*dxij[3]);
				dinv=1/d;

				for (p=1;p<=3;++p) dvij[p]=vp[i][p]-vp[j][p];

				if (d<(R[i]+R[j]))
				{
					al=((R[i]+R[j])-d);					// overlap
					dal=vtang=st=0;
					for (p=1;p<=3;++p)
					{
						mij[p]=dxij[p]*dinv;				// unit normal from i to j
						dal=dal+mij[p]*(vp[i][p]-vp[j][p]);		// Amplitude of relative velocity
					}
					for (p=1;p<=3;++p)
					{
						vn[p]=dal*mij[p];				// Normal component of relative velocity
						vt[p]=(vp[i][p]-vp[j][p])-vn[p];		// Tangential component of relative velocity
					}

					vt[1]=vt[1]-((R[i]*wp[i][2]+R[j]*wp[j][2])*dxij[3]-(R[i]*wp[i][3]+R[j]*wp[j][3])*dxij[2])*dinv;		// including component of
					vt[2]=vt[2]-((R[i]*wp[i][3]+R[j]*wp[j][3])*dxij[1]-(R[i]*wp[i][1]+R[j]*wp[j][1])*dxij[3])*dinv;		// rotational velocities in
		  			vt[3]=vt[3]-((R[i]*wp[i][1]+R[j]*wp[j][1])*dxij[2]-(R[i]*wp[i][2]+R[j]*wp[j][2])*dxij[1])*dinv;		// tangential velocity

					if(j<Nboun) meff=m[i];					// Effective mass
					else meff=m[i]*m[j]/(m[i]+m[j]);
					if (fabs(mup[i]-mup[j])>0.01) muij=mup[i];
					else muij=max(mup[i],mup[j]);

					gamn=sqrt((4*kl*lne2)/(meff*(pi*pi+lne2)));		// normal gamma corresponding to colliding ptcle
					gamt=0;
					sdotr=0;

					for(p=1;p<=3;++p)
					{
						sp[i][n][p] += vt[p]*dt*0.5;			// slip in ptcle due to collision
						sdotr += sp[i][n][p]*dxij[p];			// dot product of slip with vector joining two ptcles
					}
					for(p=1;p<=3;++p)
					{
						sp[i][n][p] -= sdotr*dxij[p]*dinv*dinv;		// final slip, after removing deformation contribtion
						st += sp[i][n][p]*sp[i][n][p];
					}
					st=sqrt(st);						// Magnitude of slip

					stinv=1/st;
					for(p=1;p<=3;++p)
					{
						if(st>1e-8) tij[p] = sp[i][n][p]*stinv;		// unit vector of slip
						else tij[p] = 0;
					}

					for(p=1;p<=3;++p) vtang += vt[p]*vt[p];

					vtang=sqrt(vtang);					// magnitude of relative tangential velocity
					vtanginv=1/vtang;
					NForce=kn*al-gamn*meff*dal;				// Normal force
					TForce=-kt*st-gamt*meff*vtang;				// Tangential force
					slip=0;

					if(fabs(TForce)>muij*fabs(NForce))
					{
						st_prev=st;
						st_previnv=1/st_prev;
						st=(-TForce/fabs(TForce)*muij*fabs(NForce)-gamt*meff*vtang)*ktinv;
						TForce=TForce/fabs(TForce/(muij*NForce));		// Tforce if its magnitude is more than Nforce

						for (p=1;p<=3;++p)
						{
							sp[i][n][p] *= st*st_previnv;		// New slip
							tij[p]=vt[p]*vtanginv;			// New slip unit vector
						}
						slip=1;
					}

					for (p=1;p<=3;++p)
					{
				  		fmp[i][p]=fmp[i][p]+mij[p]*NForce+tij[p]*TForce;	// directional forces on each particle
						fmp[j][p]=fmp[j][p]-mij[p]*NForce-tij[p]*TForce;
					}
					tijdotvt=0;

					for (p=1;p<=3;++p)
					{
						tijdotvt += tij[p]*vt[p]*slip;				// dot product of slip unit vector with relative tangential velocity
					}

					if (tijdotvt<0)
					{
						printf("WAIT: tijdotvt negative!!!! i=%d j=%d t=%f\n",i,j,t);
					}

					torquep[i][1]=torquep[i][1]-R[i]*dinv*TForce*(dxij[2]*tij[3]-dxij[3]*tij[2]);	// torque due to collision
					torquep[i][2]=torquep[i][2]-R[i]*dinv*TForce*(dxij[3]*tij[1]-dxij[1]*tij[3]);
					torquep[i][3]=torquep[i][3]-R[i]*dinv*TForce*(dxij[1]*tij[2]-dxij[2]*tij[1]);

					torquep[j][1]=torquep[j][1]-R[j]*dinv*TForce*(dxij[2]*tij[3]-dxij[3]*tij[2]);
					torquep[j][2]=torquep[j][2]-R[j]*dinv*TForce*(dxij[3]*tij[1]-dxij[1]*tij[3]);
					torquep[j][3]=torquep[j][3]-R[j]*dinv*TForce*(dxij[1]*tij[2]-dxij[2]*tij[1]);

					if(t>count*delt2)
					{
//------------------------------------------- y-directional force ------------------------------------------------
						fn[k][i-l] += NForce*mij[2];
						ft[k][i-l] += TForce*tij[2];  			// y-component of force contributes to drag
// ---------------------------------------------------------------------------------------------------------------
						tijdotvt_cal = 0;
						for (p=1;p<=3;++p)
						{
							xcontact[p] = x[j][p]+R[j]*mij[p];		// Point of contact
							vmid[p] = (v[i][p]+v[j][p])/2;			// Velocity of point of contact
							tijdotvt_cal += tij[p]*vt[p];			// dot product of slip unit vector with relative tangential velocity
						}

						fprintf(outf2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n",xcontact[1],xcontact[2],xcontact[3],vmid[1],vmid[2],vmid[3],dxij[1],dxij[2],dxij[3],dvij[1],dvij[2],dvij[3],tij[1],tij[2],tij[3],NForce,TForce,tijdotvt_cal,slip,meff);
						Ncontacts += 1;
					}
		      		}
				else for(p=1;p<=3;++p) sp[i][n][p]=0;
			}
		}
		l = l + Nspecies[k];
	}

	if(t>count*delt2)
	{
		k = Nboun;
		outf=fopen("./Output/yforce.txt","a");				// Average y-component of form and friction drag
		fprintf(outf,"%f",t);
		for(j=1;j<=N_species;++j)
		{
			fnavg[j]=ftavg[j]=fnavg[j]=ftavg[j]=0;

			for(i=k+1;i<=k+Nspecies[j];++i)
			{
				fnavg[j] += fn[j][i-k];
				ftavg[j] += ft[j][i-k];
			}
			k = k + Nspecies[j];

			fnavg[j] = fnavg[j]/Nspecies[j];
			ftavg[j] = ftavg[j]/Nspecies[j];

			fprintf(outf," %lf %lf", fnavg[j], ftavg[j]);
		}
		fprintf(outf,"\n");
		fclose(outf);
		count++;
	}
}


#endif // FORCECAL_H_INCLUDED
