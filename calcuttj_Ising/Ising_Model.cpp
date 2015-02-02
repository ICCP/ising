#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstdlib> 
#include <time.h> 
#include <stdio.h>
#include <vector>

using namespace std;

int main()
{ 
  int N=100;
  double kT;
  double max_kT=6;
  ofstream output_file;
  output_file.open("M_AvevskT.txt");
  int max_samples=100000;
  vector<vector<int> > lattice(N,vector<int>(N,1) );
  srand(time(NULL)); 
  bool temp_check = false; 
  for(kT = 1;kT<max_kT;kT+=.1)
    {
      double beta = 1/(kT);
      printf("Temp: %f\n",kT);
      
      

      double mag_avg=0;
      double mag_array[50];
      for(int runs =0; runs<50;runs++)
	{

	  for (int count_1=0; count_1<N;count_1++)
	    {
	      for (int count_2 = 0; count_2<N;count_2++)
		{lattice[count_1][count_2] = 1;}	   
	    }
	  
	  for(int sample=0;sample<max_samples;sample++)
	  {
	    int indx_1 = rand() % N;
	    int indx_2 = rand() % N;
	    double delta_E = 0;
	    
	      
	      if (indx_1 == 0)
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[N-1][indx_2]);}
	      else
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1-1][indx_2]);}
	      
	      if (indx_1 == N-1)
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[0][indx_2]);}
	      else
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1+1][indx_2]);}
	      
	      if (indx_2 == 0)
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1][N-1]);}
	      else
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1][indx_2-1]);}
	      
	      if (indx_2 == N-1)
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1][0]);}
	      else
		{delta_E+=(2*lattice[indx_1][indx_2]*lattice[indx_1][indx_2+1]);}
	      

	    
	    /*if (indx_1 != 0)
	      {delta_E+=(2*lattice_array[indx_1][indx_2]*lattice_array[indx_1-1][indx_2]);}
	    if (indx_1 !=(N-1))
	      {delta_E+=(2*lattice_array[indx_1][indx_2]*lattice_array[indx_1+1][indx_2]);}
	    if (indx_2 !=0)
	      {delta_E+=(2*lattice_array[indx_1][indx_2]*lattice_array[indx_1][indx_2-1]);}
	    if (indx_2 !=(N-1))
	      {delta_E+=(2*lattice_array[indx_1][indx_2]*lattice_array[indx_1][indx_2+1]);}
	    */
	    if (delta_E<0)
	      {lattice[indx_1][indx_2] = -1*lattice[indx_1][indx_2];}
	    else
	      {
		double R = ((double) rand()/(RAND_MAX));
		double boltz_factor = exp(-1*beta*delta_E);

		if (boltz_factor > R)
		  {lattice[indx_1][indx_2] = -1*lattice[indx_1][indx_2];}
	      }	  
	  }
	  double spin_sum = 0;
	  for(int count_1 = 0; count_1<N; count_1++)
	    {
	      for(int count_2 = 0; count_2<N;count_2++)
		{spin_sum += lattice[count_1][count_2];}
	    }
	  mag_avg +=spin_sum/(50*N*N);
	  mag_array[runs] = spin_sum/(N*N); 	  
	}
      double variance = 0;
      for(int counter = 0; counter < 50; counter++)
	{
	  variance += pow((mag_array[counter]-mag_avg),2)*.02;
	}
      double std_dev = pow(variance,.5);
      if(temp_check ==false)
	{
	  if (mag_avg-std_dev*.5 < 0)
	    {
	      cout << "FOUND CURIE POINT: " << kT <<"\n";
	      temp_check = true;
	    }
	}
      output_file << kT << " "<< mag_avg << " " << mag_avg+std_dev*.5<<" "<<mag_avg-std_dev*.5 << "\n";
    }
  output_file.close();
}


