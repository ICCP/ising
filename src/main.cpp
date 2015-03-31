/*
 * File:   main.cpp
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 5, 2015, 6:03 PM
 */

#include <cstdlib>
#include "SimulationParametersStruct.cpp"
#include "CSimulation.h"
#include <fstream>
#include <streambuf>
#include <iostream>
#include <vector>
#include <string.h>
#include <limits>

using namespace std;
/*
 *
 */

unsigned int calc_timesteps(SimulationParameters *SimPar)
{
	unsigned int t = 1;
	t *= SimPar->HEIGHT;
	t *= SimPar->WIDTH;
	t *= SimPar->temperature;
	t /= 2;
	return t;
}

void print_header(SimulationParameters *SimPar, ostream &head_stream)
{
	head_stream << "Simulation Parameters\n";
	head_stream << SimPar->WIDTH  << " WIDTH \n";
	head_stream << SimPar->HEIGHT  << " HEIGHT \n";
	head_stream << SimPar->MAXTIME  << " MAXTIME\n";
	head_stream << SimPar->TIMESTEP  << " TIMESTEP\n";
	head_stream << SimPar->MAXTIMESTEPS  << " MAXTIMESTEPS\n";
	head_stream << SimPar->COUPLING_CONSTANT  << " COUPLING_CONSTANT\n";
	head_stream << SimPar->temperature  << " temperature\n";
	head_stream << SimPar->BOLTZMAN_CONSTANT  << " BOLTZMAN_CONSTANT\n";
	head_stream << SimPar->PRINT_FREQ  << " PRINT_FREQ\n";
	head_stream << SimPar->MIN_TEMP  << " MIN_TEMP\n";
	head_stream << SimPar->MAX_TEMP  << " MAX_TEMP\n";
	head_stream << SimPar->TEMP_STEP  << " TEMP_STEP\n";
	head_stream << SimPar->CRIT_MODE  << " CRIT_MODE\n";
	head_stream << SimPar->CRIT_REPEATS  << " CRIT_REPEATS\n";
	head_stream << SimPar->var_timestep  << " var_timestep\n";
	return;
}

int main(int argc, char** argv) {
	 srand (time(NULL));
	//output stream for position
	streambuf * pos_buf;
	ofstream pos_of;
	streambuf * mag_buf;
	ofstream mag_of;
	streambuf * en_buf;
	ofstream en_of;
	streambuf * head_buf;
	ofstream head_of;
	if (false){
		//set to cout
		pos_buf = cout.rdbuf();
		mag_buf = cout.rdbuf();
		en_buf = cout.rdbuf();
		head_buf = cout.rdbuf();
	}
	else{
		//set to filename
		pos_of.open("data/raw_positions.txt");
		pos_buf = pos_of.rdbuf();
		mag_of.open("data/mag.txt");
		mag_buf = mag_of.rdbuf();
		en_of.open("data/energy.txt");
		en_buf = en_of.rdbuf();
		head_of.open("data/parameters.txt");
		head_buf = head_of.rdbuf();
	}


	ostream pos_stream(pos_buf);
	ostream mag_stream(mag_buf);
	ostream en_stream(en_buf);
	ostream head_stream(head_buf);
	//TODO allow for outputting to file
	//TODO allow for command line parameters


	SimulationParameters* SimPar = new SimulationParameters();

	print_header(SimPar, head_stream);
	head_of.close();

	if(SimPar->CRIT_MODE){
		float TempCritsAvg;
		for(unsigned short r =0;r<SimPar->CRIT_REPEATS;r++){
			cout << "\nRepeat: " << r << "/" << SimPar->CRIT_REPEATS << "\n";
			vector<float> mags;
			for(float T = SimPar->MIN_TEMP;T<SimPar->MAX_TEMP;T+=SimPar->TEMP_STEP){
				cout << "Temperature: " << T << endl;
				SimPar->temperature = T;

				//TODO SimPar->MAXTIMESTEPS=calc_timesteps(SimPar);
				cout << "MAXTIMESTEPS = " << SimPar->MAXTIMESTEPS << "\n";
				CSimulation Simulation(SimPar);
				CLattice lattice(SimPar->HEIGHT,SimPar->WIDTH, SimPar->COUPLING_CONSTANT);
				CObject *lattice2 = &lattice;
				Simulation.addObject(lattice2);

				Simulation.run(pos_stream,mag_stream, SimPar);
				pos_of.flush();
				mag_of.flush();
				en_of.flush();
				mags.push_back(Simulation.calcMag());
				//cout << "Final Energy: " << Simulation.calcEnergy() << "\n";

			}

			//calculate crit temp
			unsigned short min_sus_idx = 0;//index in terms of mag
			float min_sus = numeric_limits<float>::infinity();
			for(unsigned short i = 1; i<mags.size()-1;i++){
				float new_sus = (mags[i+1]-mags[i-1])/(2*SimPar->TEMP_STEP);
				if(new_sus < min_sus){
					min_sus = new_sus;
					min_sus_idx = i;
				}
			}
			TempCritsAvg+=(SimPar->MIN_TEMP+SimPar->TEMP_STEP*min_sus_idx);
			cout << (SimPar->MIN_TEMP+SimPar->TEMP_STEP*min_sus_idx) << "\n";
		}
		TempCritsAvg/=SimPar->CRIT_REPEATS;
		cout << "Final Average Critical Temp: " << TempCritsAvg << endl;

	}
	else{
		for(float T = SimPar->MIN_TEMP;T<SimPar->MAX_TEMP;T+=SimPar->TEMP_STEP)
		{
			SimPar->temperature = T;
			cout << "\nTemperature: " << T << endl;
			en_stream << T << "\n";

			//TODO SimPar->MAXTIMESTEPS=calc_timesteps(SimPar);
			cout << "MAXTIMESTEPS = " << SimPar->MAXTIMESTEPS << "\n";
			CSimulation Simulation(SimPar);
			CLattice lattice(SimPar->HEIGHT,SimPar->WIDTH, SimPar->COUPLING_CONSTANT);
			CObject *lattice2 = &lattice;
			Simulation.addObject(lattice2);

			Simulation.run(pos_stream,mag_stream, SimPar);
			pos_of.flush();
			mag_of.flush();
			en_of.flush();
			float energy = Simulation.calcEnergy();
			cout << "Final Energy: " << energy << "\n";
			en_stream << energy << "\n";
		}
	}
	delete SimPar;
    return 0;
}
