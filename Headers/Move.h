#ifndef MOVE_H_INCLUDED
#define MOVE_H_INCLUDED

void Move3D(void)									// Applying velocity-verlet algorithm
{
	short i, p;
	double dt2;

	dt2 = dt*0.5;
	ForceCalc(x,v,fm,w,s,torque);
	for (i=Nboun+1;i<=Npartstot;++i)
	{
		for(p=1;p<=3;p++)							// Updating position and velocities
		{
			x[i][p] = x[i][p] + (v[i][p] + 0.5*dt*fm[i][p]*m_inv[i])*dt;
			v[i][p] = v[i][p] + fm[i][p]*dt2*m_inv[i]; 			//Note dt2 here
			w[i][p] = w[i][p] + torque[i][p]*dt2*I_m_inv[i];
		}
	}
	ForceCalc(x,v,fm,w,s,torque);
	for (i=Nboun+1;i<=Npartstot;++i)
	{
		for(p=1;p<=3;p++)
		{
			v[i][p] = v[i][p] + fm[i][p]*dt2*m_inv[i];
			w[i][p] = w[i][p] + torque[i][p]*dt2*I_m_inv[i];
		}
	}
}

void Move2D(void)									// Applying velocity-verlet algorithm
{
	short i, p;
	double dt2;

	dt2 = 0.5*dt;
	ForceCalc(x,v,fm,w,s,torque);
	for(i=Nboun+1;i<=Npartstot;++i)
	{
		for(p=1;p<=2;p++)							// Updating position and velocities
		{
			x[i][p] = x[i][p] + (v[i][p] + 0.5*dt*fm[i][p]*m_inv[i])*dt;
			v[i][p] = v[i][p] + fm[i][p]*dt2*m_inv[i]; 			//Note dt2 here
		}
		w[i][3] = w[i][3] + torque[i][3]*dt2*I_m_inv[i];
	}
	ForceCalc(x,v,fm,w,s,torque);
	for (i=Nboun+1;i<=Npartstot;++i)
	{
		for(p=1;p<=2;p++)
		{
			v[i][p] = v[i][p] + fm[i][p]*dt2*m_inv[i];
		}
		w[i][3] = w[i][3] + torque[i][3]*dt2*I_m_inv[i];
	}
}


#endif // FORCECAL_H_INCLUDED
