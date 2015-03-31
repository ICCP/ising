/*
 * SimulationParametersStruct.cpp
 *
 *  Created on: Mar 13, 2015
 *      Author: Mark Wittekind and Drew Murray
 */
#ifndef SIM_PARAMS_UNIQUE
#define SIM_PARAMS_UNIQUE

struct SimulationParameters
{
	const unsigned short WIDTH = 50;
	const unsigned short HEIGHT = 50;
	unsigned short MAXTIME = 50;
	float TIMESTEP = 1;
	unsigned int MAXTIMESTEPS = short(MAXTIME/TIMESTEP);
	float COUPLING_CONSTANT = 1;
	float temperature = 2;
	float BOLTZMAN_CONSTANT = 1;
	unsigned int PRINT_FREQ = 1;
	float MIN_TEMP = 1;
	float MAX_TEMP = 3;
	float TEMP_STEP = .25;
	bool CRIT_MODE = false;
	unsigned short CRIT_REPEATS = 10;
	bool var_timestep = true;

	//TODO
	// Critical exponents:
	// 		log(X) vs log(T)
	// 		log(mag) vs log(T)
};

#endif

