#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <stack>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <fstream>
using namespace std;

typedef long long ll;
typedef long double ld;
#define VI vector<int>
#define VS  vector<string>
#define PII pair<int,int>
#define FOR(i,n) for(int i=0;i<n;i++)
#define FO(i,a,b) for(int i=a;i<b;i++)
#define FOX(i,a) for(int i=0;i<a.size();i++)
#define X first
#define Y second
#define PB push_back
#define MP make_pair

ld length, length_each_elem;
int num_elem;

int *boundary;
/*
  0 - free node
  1 - Hinged
  2 - Fixed
*/

// E & I for each of num_elems
ld *E,*I; 


int num_forces;
ld **forces;
/*
  0 - Point load
      x - location
      p - value

  1 - Distributed load
      n1 - Starting node
      n2 - Ending node
      a  - Force at n1
      b  - Force at n2
*/  

ld **K,*F,*U,*Reactions;
int *u;

ld interpolate(ld *a,ld n)
{
  return a[3] + (n-a[1])*((a[4]-a[3])/(a[2]-a[1]));
}


vector<ld> gauss_elem(vector< vector<ld> > final);

int main(int argc,char **argv)
{
  ifstream in(argv[1]);
  ofstream out(argv[2]);
  

  /////////////////////// Input Routines  \\\\\\\\\\\\\\\\\\\\\\\

  in>>length>>num_elem;

  //num_elem+1 nodes
  boundary = new int[num_elem+1];
  FOR(i,num_elem+1)
    in>>boundary[i];

  // E & I for num_elems 
  E = new ld [num_elem];
  I = new ld [num_elem];
  FOR(i,num_elem)
    in>>E[i]>>I[i];
  // Length of each elem = (Total length)/(num_elem)
  length_each_elem = length/num_elem;


  in>>num_forces;
  forces = new ld* [num_forces];
  FOR(i,num_forces)
    {
      forces[i] = new ld [5];  // (type,x1,x2,a,b) or (type,x,p)
      int type;
      in>>type;
      if(!type)
	forces[i][0]=0,
	  in>>forces[i][1]>>forces[i][2];
      else
	forces[i][0]=1,
	  in>>forces[i][1]>>forces[i][2]>>forces[i][3]>>forces[i][4],
	  forces[i][1]--,forces[i][2]--;
    }

  ///////////////////// End of Inputs  \\\\\\\\\\\\\\\\\\\\\\\\

  ///////////////////// Stiffness Matrix \\\\\\\\\\\\\\\\\\\\\\
  
  // Note: For X elements, size of stiffness matrix is (2X+2) by (2X+2)

  K = new ld* [2*(num_elem+1)];
  FOR(i,2*(num_elem+1))
    {
      K[i] = new ld[2*(num_elem+1)];
      FOR(j,2*(num_elem+1))
	K[i][j]=0;
    }
  
  FOR(i,num_elem)
    {
      ld temp = 2*E[i]*I[i]/pow(length_each_elem,3);

      K[2*i][2*i]     += 6*temp;
      K[2*i][2*i+1]   += -3*length_each_elem*temp;
      K[2*i][2*i+2]   += -6*temp;
      K[2*i][2*i+3]   += -3*length_each_elem*temp;

      K[2*i+1][2*i]   += -3*length_each_elem*temp;
      K[2*i+1][2*i+1] += 2*pow(length_each_elem,2)*temp;
      K[2*i+1][2*i+2] += 3*length_each_elem*temp;
      K[2*i+1][2*i+3] += pow(length_each_elem,2)*temp;

      K[2*i+2][2*i]   += -6*temp;
      K[2*i+2][2*i+1] += 3*length_each_elem*temp;
      K[2*i+2][2*i+2] += 6*temp;
      K[2*i+2][2*i+3] += 3*length_each_elem*temp;

      K[2*i+3][2*i]   += -3*length_each_elem*temp;
      K[2*i+3][2*i+1] += pow(length_each_elem,2)*temp;
      K[2*i+3][2*i+2] += 3*length_each_elem*temp;
      K[2*i+3][2*i+3] += 2*pow(length_each_elem,2)*temp;
    }
  //////////////// Phew! Stiffness matrix done.. \\\\\\\\\\\\\\\\\

  ///////////////////// Force Matrix \\\\\\\\\\\\\\\\\\\

  F = new ld [ 2*(num_elem+1)];
  bool *done = new bool [num_forces];
  
  FOR(i,2*(num_elem+1))
    F[i] = 0;
  FOR(i,num_forces)
    done[i]=0;
  
  FOR(j,num_forces)
    if(forces[j][0]==0)
      FOR(i,num_elem+1)
	if(abs(forces[j][1] - i*length_each_elem)<1e-9)
	  F[2*i] += forces[j][2],
	    done[j]=1;
	
  
  FOR(i,num_elem)
    FOR(j,num_forces)
    if(!done[j])
      if(forces[j][0] == 0)  // Point Load.. Shift the point load to left node & add moments
	{
	  if(forces[j][1] > i*length_each_elem and forces[j][1] < (i+1)*length_each_elem)
	    {
	      F[2*i] += forces[j][2];
	      F[2*i+1] -= forces[j][2]*(forces[j][1] - i*length_each_elem);
	      F[2*i+3] += forces[j][2]*(forces[j][1] - i*length_each_elem);
	    }
	}
      else                   // Line Load
	{
	  if(forces[j][1]<=i && forces[j][2]>=i+1)
	    {
	      ld average = (interpolate(forces[j],i)+interpolate(forces[j],i+1))/2;
	      F[2*i] += average*length_each_elem/2;
	      F[2*i+2] += average*length_each_elem/2;
	      F[2*i+1] -= average*pow(length_each_elem,2)/12;
	      F[2*i+3] += average*pow(length_each_elem,2)/12;
	    }
	}
  ///////////////////// Done Force Matrix \\\\\\\\\\\\\\\\\\	\

  ///////////////////// Displacement Matrix \\\\\\\\\\\\\\\\\\
  
  u = new int [2*(num_elem+1)];
  FOR(i,num_elem+1)
    if(boundary[i]==0)
      u[2*i] = u[2*i+1]=1;
    else if(boundary[i] == 1)
      u[2*i] = 0, u[2*i+1] = 1;
    else
      u[2*i] = u[2*i+1] = 0;
      
  ///////////////////// Done Displacement Matrix \\\\\\\\\\\\\\\


  //////////////////// Get ready to solve \\\\\\\\\\\\\\\\\\\\\

  vector<vector<ld> > matrix;
  FOR(i,2*(num_elem+1))
    if(u[i])
      {
	vector<ld> row;
	FOR(j,2*(num_elem+1))
	  if(u[j])
	    row.PB(K[i][j]);
	row.PB(F[i]);
	matrix.PB(row);
      }
  //////////////////////// Solve it! \\\\\\\\\\\\\\\\\\\\\\	\
  
  vector<ld> unknown_disp = gauss_elem(matrix);
  
  U = new ld [2*(num_elem+1)];
  Reactions = new ld [2*(num_elem+1)];

  int count=0;
  FOR(i,2*(num_elem+1))
    if(u[i])
      U[i] = unknown_disp[count++];
    else
      U[i] = 0;
  FOR(i,2*(num_elem+1))
    if(u[i]==0)
      {
	Reactions[i]=0;
	FOR(j,2*(num_elem+1))
	  if(u[j] == 1)
	    Reactions[i] += K[i][j]*U[j];
      }
    else
      Reactions[i]=0;
  

  //////////////////////// Done! Print the output \\\\\\\\\\\\\\\\\	\

  
  out<<setw(5)<<"Node"<<setw(10)<<"U"<<setw(10)<<"Theta"<<setw(10)<<"R"<<setw(10)<<"M"<<endl<<endl;
  
  FOR(i,num_elem+1)
    {
      out<<setw(5)<<i+1<<setw(10)<<U[2*i]<<setw(10)<<U[2*i+1];
      if(u[2*i])
	out<<setw(10)<<"";
      else
	out<<setw(10)<<Reactions[2*i];

      if(u[2*i+1])
	out<<setw(10)<<"";
      else
	out<<setw(10)<<Reactions[2*i+1];
      out<<endl;
    }
								    
								    
  in.close();
  out.close();
   
  //////////////////////// Bye Bye \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  return 0;
}

vector<ld> gauss_elem(vector< vector<ld> > final)
{
  int nv = final[0].size()-1;
  /*  FOX(i,final)
    {
      FOX(j,final[0])
	cout<<final[i][j]<<" ";
      cout<<endl;
    }
  */
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
  vector<ld> ans;
  ans.resize(nv);
  FOR(i,nv)
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
