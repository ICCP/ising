#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <fstream>

int main()
{
  const int N = 100;
  double T_initial = 0.1;
  double T_final = 6.0;
  double T_step_size  = 0.1;
  const int T_steps = (T_final - T_initial) / (T_step_size);
  int Init;
  const int iterations = 15000;
  const int repetitions = 50;
  double mag = 0;
  double mag_ave = 0;
  double mag_dev = 0;
  double mag_devplus = 0;
  double mag_devminus = 0;
  double delta_E;
  double boltzmann = 0;
  int lattice[N][N];
  double results[T_steps][4];

  srand(time(NULL));

  //  std::cout << "Size of lattice";
  //  std::cin >> N;

  std::cout << "For initialization of the Spins enter: Down = -1, Random = 0, Up = 1" << std::endl;
  std::cin >>  Init;
  for (int w = 0 ; w <= T_steps ; w++)
    {
      double T = T_initial+T_step_size*w;
      results[w][0] = T;
      std::vector<double> TemporaryMag;      

      for (int reps = 1 ; reps <= repetitions ; reps++)
	{
	  for (int indx_1 = 0; indx_1 < N; indx_1++)            // fill the NxN matrix with 1's
	    {
	      for (int indx_2 = 0; indx_2 < N; indx_2++)
		{
		  if (Init == 1)
		    {lattice[indx_1][indx_2] = 1;}
		  if (Init == -1)
		    {lattice[indx_1][indx_2] = -1;}
		  if (Init == 0)
		    {
		      int rando_calrissian = int(rand()) % 2;
		      if (rando_calrissian == 0)
			{lattice[indx_1][indx_2] = -1;}
		      else
			{lattice[indx_1][indx_2] = 1;}
		    }
		}
	    }

	  for (int X = 0 ; X < iterations; X++)              // Repeat for X iterations where X is defined above
	    { 
	      int indx_i = rand() % N;                      // Choose random i index
	      int indx_j = rand() % N;                      // Choose random j index
	      delta_E = 0;
	  
	      if (indx_i !=0)
		{delta_E+=(2*lattice[indx_i][indx_j]*lattice[indx_i-1][indx_j]);}
	      if (indx_i !=N-1)
		{delta_E+=(2*lattice[indx_i][indx_j]*lattice[indx_i+1][indx_j]);}
	      if (indx_j !=0)
		{delta_E+=(2*lattice[indx_i][indx_j]*lattice[indx_i][indx_j-1]);}
	      if (indx_j !=N-1)
		{delta_E+=(2*lattice[indx_i][indx_j]*lattice[indx_i][indx_j+1]);}

	      if (delta_E < 0)                              // Check for Spin Inversion
		{
		  lattice[indx_i][indx_j] = -1*lattice[indx_i][indx_j];
		}
	      else
		{
		  double R = ((double) rand() / (RAND_MAX));
		  double beta = 1/(T);
		  boltzmann = exp(-1*beta*delta_E);
		  if (boltzmann > R)
		    {
		      lattice[indx_i][indx_j] = -1*lattice[indx_i][indx_j];
		    }
		}
	      double total_sum = 0;                        // Sum up all the Spins
	      for (int indx_1 = 0; indx_1 < N; indx_1++)
		{
		  for (int indx_2 = 0; indx_2 < N; indx_2++)
		    {total_sum += lattice[indx_1][indx_2];}
		}
	      mag = total_sum/(N*N);
	    }

	  TemporaryMag.push_back(mag);

	} //end of repetitions
 
      mag_ave = 0;
      mag_dev = 0;

      for (int i=0;i<TemporaryMag.size();i++)
	{
	  mag_ave += TemporaryMag[i];
	}
      mag_ave /= TemporaryMag.size();

      for (int i=0 ; i < TemporaryMag.size() ; i++)
	{
	  mag_dev += pow( pow(TemporaryMag[i]-mag_ave,2) / TemporaryMag.size()  ,0.5);
	}
      mag_dev = pow(mag_dev / TemporaryMag.size() , 0.5);

      mag_devplus = mag_ave + mag_dev;
      mag_devminus = mag_ave - mag_dev;

      results[w][1] = mag_ave;
      results[w][2] = mag_devminus;
      results[w][3] = mag_devplus;
    } //end of temps
  std::ofstream myfile;
  myfile.open ("ising.txt");
  for (int i = 0; i <= T_steps; i++)     
    {
      for (int j = 0; j < 4; j++)
	{
	  myfile << results[i][j] << ' ';
	}
      myfile << std::endl;
    }
  return 0;
}
