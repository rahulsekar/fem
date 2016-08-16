#include <vector>
#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;

long double Length,ElementLength;
long double E,I; 
long double **Forces;
long double **K,*F,*U,*Reactions;

int NumberElements,NumberNodes,TotalDOF,NumberBoundary,NumberForces;
bool *Boundary;

vector<long double> gauss_elem(vector< vector<long double> > final);

char title[1000],dummy[1000];

int main(int argc,char **argv)
{
  FILE *in,*out;
  in = fopen(argv[1],"r");
  out = fopen(argv[2],"w");
  
  /* Take Inputs */
  fgets(title,1000,in);
  fgets(dummy,1000,in);
  fscanf(in,"%Lf",&Length);
  fgets(dummy,1000,in);
  fgets(dummy,1000,in);
  fscanf(in,"%Lf %Lf",&E,&I);
  fgets(dummy,1000,in);
  fgets(dummy,1000,in);
  fscanf(in,"%d",&NumberElements);

  NumberNodes = NumberElements+1;
  TotalDOF = 2*NumberNodes;
  
  Boundary = new bool[TotalDOF];
  for(int i=0;i<TotalDOF;i++)
    Boundary[i]=0;

  fgets(dummy,1000,in);
  fgets(dummy,1000,in);
  fscanf(in,"%d",&NumberBoundary);
  fgets(dummy,1000,in);
  fgets(dummy,1000,in);
  for(int i=0;i<NumberBoundary;i++)
    {
      int nodeno,mmnt,frce;
      fscanf(in,"%d %d %d",&nodeno,&frce,&mmnt);
      Boundary[2*nodeno] = frce;
      Boundary[2*nodeno+1] = mmnt;
    }

  ElementLength = Length/NumberElements;
  fgets(dummy,1000,in);
  fgets(dummy,1000,in);
  fscanf(in,"%d",&NumberForces);
  fgets(dummy,1000,in);
  fgets(dummy,1000,in);

  Forces = new long double* [NumberForces];
  for(int i=0;i<NumberForces;i++)
    {
      Forces[i] = new long double [5];  // (type,x1,x2,a,b) or (type,x,p)
      int type;
      fscanf(in,"%d",&type);
      if(!type)
	Forces[i][0]=0,
	  fscanf(in,"%Lf %Lf",&Forces[i][1],&Forces[i][2]);
      else
	Forces[i][0]=1,
	  fscanf(in,"%Lf %Lf %Lf %Lf",&Forces[i][1],&Forces[i][2],&Forces[i][3],&Forces[i][4]);
    }

  /* Stiffness Matrix */
  
  K = new long double* [TotalDOF];
  for(int i=0;i<TotalDOF;i++)
    {
      K[i] = new long double[TotalDOF];
      for(int j=0;j<TotalDOF;j++)
	K[i][j]=0;
    }
  
  for(int i=0;i<NumberElements;i++)
    {
      long double temp = 2*E*I/pow(ElementLength,3);

      K[2*i][2*i]     += 6*temp;
      K[2*i][2*i+1]   += -3*ElementLength*temp;
      K[2*i][2*i+2]   += -6*temp;
      K[2*i][2*i+3]   += -3*ElementLength*temp;

      K[2*i+1][2*i]   += -3*ElementLength*temp;
      K[2*i+1][2*i+1] += 2*pow(ElementLength,2)*temp;
      K[2*i+1][2*i+2] += 3*ElementLength*temp;
      K[2*i+1][2*i+3] += pow(ElementLength,2)*temp;

      K[2*i+2][2*i]   += -6*temp;
      K[2*i+2][2*i+1] += 3*ElementLength*temp;
      K[2*i+2][2*i+2] += 6*temp;
      K[2*i+2][2*i+3] += 3*ElementLength*temp;

      K[2*i+3][2*i]   += -3*ElementLength*temp;
      K[2*i+3][2*i+1] += pow(ElementLength,2)*temp;
      K[2*i+3][2*i+2] += 3*ElementLength*temp;
      K[2*i+3][2*i+3] += 2*pow(ElementLength,2)*temp;
    }

  /* Force Matrix */

  F = new long double [TotalDOF];
  bool *done = new bool [NumberForces];
  
  for(int i=0;i<NumberForces;i++)
    done[i]=0;
  for(int i=0;i<TotalDOF;i++)
    F[i]=0;
  
  for(int j=0;j<NumberForces;j++)
    for(int i=0;i<NumberNodes;i++)
      if(Forces[j][0]==0)
      if(abs(Forces[j][1] - i*ElementLength)<1e-9)
	F[2*i] += Forces[j][2],
	  done[j]=1;
	
  for(int i=0;i<NumberElements;i++)
    for(int j=0;j<NumberForces;j++)
      if(!done[j])
	if(Forces[j][0] == 0)  // Point Load.. Shift the point load to left node & add moments
	  {
	    if(Forces[j][1] > i*ElementLength and Forces[j][1] < (i+1)*ElementLength)
	      {
		F[2*i] += Forces[j][2];
		F[2*i+1] -= Forces[j][2]*(Forces[j][1] - i*ElementLength);
		F[2*i+3] += Forces[j][2]*(Forces[j][1] - i*ElementLength);
	      }
	  }
	else                   // Line Load
	  {
	    if(Forces[j][1]<=i && Forces[j][2]>=i+1)
	      {
		long double left = Forces[j][3]+((long double)(i)-Forces[j][1])*(Forces[j][4]-Forces[j][3])/(Forces[j][2]-Forces[j][1]);
		long double right = Forces[j][3]+((long double)(i+1)-Forces[j][1])*(Forces[j][4]-Forces[j][3])/(Forces[j][2]-Forces[j][1]);
		long double avg = (left+right)/2;
		F[2*i] += avg*ElementLength/2;
		F[2*i+2] += avg*ElementLength/2;
		F[2*i+1] -= avg*pow(ElementLength,2)/12;
		F[2*i+3] += avg*pow(ElementLength,2)/12;
	      }
	  }



  /* Get ready to solve */
  
  vector<vector<long double> > matrix;
  for(int i=0;i<TotalDOF;i++)
    if(!Boundary[i])
      {
	vector<long double> row;
	for(int j=0;j<TotalDOF;j++)
	  if(!Boundary[j])
	    row.push_back(K[i][j]);
	row.push_back(F[i]);
	matrix.push_back(row);
      }
  
  /* Solve */
  
  vector<long double> unknown_disp = gauss_elem(matrix);
  
  U = new long double [TotalDOF];
  Reactions = new long double [TotalDOF];

  int count=0;
  for(int i=0;i<TotalDOF;i++)
    if(!Boundary[i])
      U[i] = unknown_disp[count++];
    else
      U[i] = 0;
  for(int i=0;i<TotalDOF;i++)
    if(Boundary[i]==1)
      {
	Reactions[i]=0;
	for(int j=0;j<TotalDOF;j++)
	  if(Boundary[j] == 0)
	    Reactions[i] += K[i][j]*U[j];
      }
    else
      Reactions[i]=0;
  

  /* Print the output */

  fprintf(out,"%s\n",title);
  fprintf(out,"Node  Displacement   Angle   Reaction    Moment\n");
  for(int i=0;i<NumberNodes;i++)
    {
      fprintf(out,"%d  %Lf  %Lf ",i,U[2*i],U[2*i+1]);
      if(!Boundary[2*i])
	fprintf(out," - ");
      else
	fprintf(out," %Lf ",Reactions[2*i]);
      
      if(!Boundary[2*i+1])
	fprintf(out," - ");
      else
	fprintf(out," %Lf ",Reactions[2*i+1]);
      fprintf(out,"\n");
    }
  
  
  fclose(in);
  fclose(out);
  return 0;
}

vector<long double> gauss_elem(vector< vector<long double> > final)
{
  int nv = final[0].size()-1;

  for(int i=0;i<nv;i++)
    {
      int ind=i;
      for(int k=i+1;k<nv;k++)
	if(abs(final[k][i]) > abs(final[ind][i]))
	  ind = k;
      vector<long double> t =final[i];
      final[i] = final[ind];
      final[ind] = t;
      for(int j=i+1;j<nv;j++)
	{
	  
	  if(abs(final[j][i]) < 0.0000001)
	    {
	      final[j][i] = 0;
	      continue;
	    }
	  long double mul = final[i][i]/final[j][i];
	  for(int k=i;k<=nv;k++)
	    {
	      final[j][k] = final[j][k]*mul - final[i][k];
	      if(abs(final[j][k]) <1e-9)
		final[j][k] =0;
	    }
	}
    }
  vector<long double> ans;
  ans.resize(nv);
  for(int i=0;i<nv;i++)
    ans[i]=0;

  for(int i=nv-1;i>=0;i--)
    {
      long double c = final[i][nv];
      for(int j=i+1;j<nv;j++)
	c-= final[i][j]*ans[j];
      ans[i] = c/final[i][i];
      if(abs(ans[i]) < 1e-9)
	ans[i] = 0;
    }
  return ans;
}
