#ifndef MAKEMAPS_H_INCLUDED
#define MAKEMAPS_H_INCLUDED


void MakeMap()
{

	if(simflag==3)
	{
		int i, j, k, idx;
		for (j=1;j<=Nbiny;++j)
		{
			for (i=1;i<=Nbinx;++i)
			{
				for (k=1;k<=Nbinz;++k)
				{
					cell[i][j][k]=k+(i-1)*Nbinz+(j-1)*Nbinz*Nbinx;		// Numbering of cells, i.e. bins
				}
			}
		}

		for (j=1;j<=Nbiny;++j)
		{
			for (i=1;i<=Nbinx;++i)
			{
				for (k=1;k<=Nbinz;++k)
				{
					idx=13*(cell[i][j][k]-1);

					if(k==1)						//apply PBC in Z direction
					{
						if(j==1)
						{
							mapc[idx+1]=-1;
							mapc[idx+2]=-1;
							mapc[idx+3]=-1;
						}
						else
						{
							if(i!=1) mapc[idx+1]=cell[i-1][j-1][Nbinz];
							mapc[idx+2]=cell[i][j-1][Nbinz];
							if(i!=Nbinx) mapc[idx+3]=cell[i+1][j-1][Nbinz];
							else mapc[idx+3]=cell[1][j-1][Nbinz];
						}

						if(i==1)
						{
							if(j!=1)mapc[idx+1]=cell[Nbinx][j-1][Nbinz];
							mapc[idx+4]=cell[Nbinx][j][Nbinz];
							mapc[idx+5]=cell[Nbinx][j+1][Nbinz];
						}
						else
						{
							mapc[idx+4]=cell[i-1][j][Nbinz];
							mapc[idx+5]=cell[i-1][j+1][Nbinz];
						}
						if (i==Nbinx) mapc[idx+6]=cell[1][j][Nbinz];
						else mapc[idx+6]=cell[i+1][j][Nbinz];
						mapc[idx+7]=cell[i][j][Nbinz];
						mapc[idx+8]=cell[i][j+1][Nbinz];
						if (i==Nbinx) mapc[idx+9]=cell[1][j+1][Nbinz];
						else mapc[idx+9]=cell[i+1][j+1][Nbinz];
					}

					else						//if k!=1
					{
						if(i==1)
						{
							if(j!=1) mapc[idx+1]=cell[Nbinx][j-1][k-1];
							mapc[idx+4]=cell[Nbinx][j][k-1];
							mapc[idx+5]=cell[Nbinx][j+1][k-1];
						}
						else
						{
							if(j!=1)mapc[idx+1]=cell[i-1][j-1][k-1];
							mapc[idx+4]=cell[i-1][j][k-1];
							mapc[idx+5]=cell[i-1][j+1][k-1];
						}

						if(j==1)
						{
							mapc[idx+1]=-1;
							mapc[idx+2]=-1;
							mapc[idx+3]=-1;
						}

						else
						{
							if(i!=1) mapc[idx+1]=cell[i-1][j-1][k-1];
							mapc[idx+2]=cell[i][j-1][k-1];
							if (i==Nbinx) mapc[idx+3]=cell[1][j-1][k-1];
							else mapc[idx+3]=cell[i+1][j-1][k-1];
						}

						if (i==Nbinx) mapc[idx+6]=cell[1][j][k-1];
						else mapc[idx+6]=cell[i+1][j][k-1];
						mapc[idx+7]=cell[i][j][k-1];
						mapc[idx+8]=cell[i][j+1][k-1];
						if (i==Nbinx) mapc[idx+9]=cell[1][j+1][k-1];
						else mapc[idx+9]=cell[i+1][j+1][k-1];
					}

					if(j==1)
					{
						mapc[idx+10]=-1;
						mapc[idx+11]=-1;
						mapc[idx+12]=-1;
					}
					else
					{
						if(i==1) mapc[idx+10]=cell[Nbinx][j-1][k];
						else mapc[idx+10]=cell[i-1][j-1][k];
						mapc[idx+11]=cell[i][j-1][k];
						if(i==Nbinx) mapc[idx+12]=cell[1][j-1][k];
						else mapc[idx+12]=cell[i+1][j-1][k];
					}
					if(i==1) mapc[idx+13]=cell[Nbinx][j][k];
					else mapc[idx+13]=cell[i-1][j][k];
				}//i
			}//j
		}//k
	}//simflag


	if(simflag==2)
	{
		int i, j, idx;
		for (j=1;j<=Nbiny;++j)
		{
			for (i=1;i<=Nbinx;++i)
			{
				cell2d[i][j]=i+(j-1)*Nbinx;					// numbering of cells in 2D
			}
		}
		for (j=1;j<=Nbiny;++j)
		{
			for (i=1;i<=Nbinx;++i)
			{
				idx=4*(cell2d[i][j]-1);

				if(j==1)
				{
					mapc[idx+1]=-1;
					mapc[idx+2]=-1;
					mapc[idx+3]=-1;
				}
				else
				{
					if(i!=1) mapc[idx+1]=cell2d[i-1][j-1];
					else mapc[idx+1]=cell2d[Nbinx][j-1];
					mapc[idx+2]=cell2d[i][j-1];
					if(i!=Nbinx) mapc[idx+3]=cell2d[i+1][j-1];
					else mapc[idx+3]=cell2d[1][j-1];
				}

				if(i==Nbinx) mapc[idx+4]=cell2d[1][j];
				else mapc[idx+4]=cell2d[i+1][j];
			}
		}
	}//simflag

}//program



#endif // TEST_H_INCLUDED



