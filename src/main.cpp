#include <vector>
#include <cmath>
//#include <boost/tuple/tuple.hpp>
#include <random>
#include <iostream>
#include <math.h>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
//#include "gnuplot-iostream.h"

// TODO: Compare time avg to ensemble

////////////
// GLOBAL //
////////////
int genSeed = time(0);
std::default_random_engine engine(genSeed);
std::uniform_real_distribution<float> rDist(0, 1);

using std::endl;
using std::cout;
using std::vector;
using std::chrono::steady_clock;
using std::chrono::duration;

vector<vector<int >> fillArray(vector<vector<int >> lattice, int N, bool rand, float& magnetization) {
	std::uniform_int_distribution<int> dist(-1, 1);
	magnetization = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (rand == 1)
			{
				lattice[i][j] = dist(engine);
			}
			else {
				lattice[i][j] = 1;
			}
			magnetization += lattice[i][j];
		}
	}
	//cout << magnetization << endl;
	return lattice;
}

void metroTest(float energyDiff, vector<vector<int >>& lattice, int row, int col, float& magnetization, float currentTemp)
{
	const int boltz = 1;
	double beta;
	// R value random distribution
	float R;

	// Metropolis test
	if (energyDiff < 0)
	{
		// Flip state
		lattice[row][col] = -lattice[row][col];
		magnetization += 2 * lattice[row][col];

	}
	else {
		// Pick random value for R
		R = rDist(engine);
		beta = 1 / (boltz * currentTemp);

		if (exp(-beta * (energyDiff)) > R) {
			// Flip state
			lattice[row][col] = -lattice[row][col];
			magnetization += 2 * lattice[row][col];
		}
	}
}

float findEnergyDiff(int row, int col, vector<vector<int >> lattice, int N) {
	float energyDiff = 0;

	// Check left
	if (col == 0) {
		//energyDiff += lattice[row][col] * lattice[row][N-1];        
	}
	else {
		energyDiff += lattice[row][col] * lattice[row][col - 1];
	}

	// Check right
	if (col == N - 1) {
		//energyDiff += lattice[row][col] * lattice[row][0];

	}
	else {
		// Check col +1
		energyDiff += lattice[row][col] * lattice[row][col + 1];
	}

	// Check up
	if (row == N - 1) {
		// Check 0
		//energyDiff += lattice[row][col] * lattice[0][col];

	}
	else {
		// Check row +1
		energyDiff += lattice[row][col] * lattice[row + 1][col];
	}

	// Check down
	if (row == 0) {
		// Check N-1
		//energyDiff += lattice[row][col] * lattice[N-1][col];

	}
	else {
		// Check row -1
		energyDiff += lattice[row][col] * lattice[row - 1][col];

	}

	energyDiff = 2 * energyDiff;
	return energyDiff;
}

float averageArray(float myArray[], int size)
{
	float sum = 0;
	for (int i = 0; i < size; i++)
	{
		sum += myArray[i];
	}
	return sum / size;
}

float findStdDev(float avg, float myArray[], int size)
{
	float variance = 0;
	for (int i = 0; i < size; i++)
	{
		variance += pow((avg - myArray[i]), 2);
	}
	return sqrt(variance / size);
}

void generateMathematica(vector<vector<float>>& results)
{
	std::ofstream myfile;
	myfile.open("results.txt");
	myfile << "ErrorListPlot[{{";
	for (int i = 0; i < results.size(); i++)
	{
		myfile << results[i][1];
		if (i < results.size() - 1)
		{
			myfile << ",";
		}
		else {
			myfile << "},{";
		}
	}

	for (int i = 0; i < results.size(); i++)
	{
		if (i < results.size() - 1)
		{
			myfile << results[i][2] << ", ";
		}
		else {
			myfile << results[i][2] << "}} // Transpose, Epilog -> {Text[\"sine\", Scaled[{.2, .8}]]}, PlotRange -> {{0, 60}, {-.2, 1.2}}]";
		}
	}
	myfile.close();
}

void reportTime(int iter, steady_clock::time_point& currentTime, steady_clock::time_point& lastTime, int currentIteration, int totalIterations)
{
	currentTime = steady_clock::now();
	duration<double> time_span = std::chrono::duration_cast<duration<double>>(currentTime - lastTime);
	double remainingTime = (((time_span.count()) / iter) * (totalIterations - currentIteration));
	cout << "Iterations remaining: " << totalIterations - currentIteration << endl;

	if (remainingTime > 60)
	{
		if (remainingTime > 3600)
		{
			int hours = remainingTime / 3600;
			int minutes = fmod(remainingTime, 3600) / 60;
			int seconds = fmod(fmod(remainingTime, 60), minutes);
			cout << "Remaining Time: " << hours << " hours " << minutes << " minutes " << seconds << " seconds." << endl;
		}
		else
		{
			int minutes = fmod(remainingTime, 3600) / 60;
			int seconds = fmod(fmod(remainingTime, 60), minutes);
			cout << "Remaining Time: " << minutes << " minutes " << seconds << " seconds." << endl;
		}
	}
	else {
		int seconds = remainingTime;
		cout << "Remaining Time: " << seconds << " seconds." << endl;
	}
}

//void plotRealTime(Gnuplot latticePlot, vector<vector<int >> lattice)
//{
//	// Display state
//	latticePlot << "unset key\n";
//	latticePlot << "set pm3d\n";
//	latticePlot << "set hidden3d\n";
//	latticePlot << "set view map\n";
//	latticePlot << "set xrange [ 0 : " << N << " ] \n";
//	latticePlot << "set yrange [ 0 : " << N << " ] \n";
//	latticePlot << "splot '-'\n";
//	latticePlot.send2d(lattice);
//	latticePlot.flush(); 
//}

int main() {
	//Gnuplot latticePlot;
	//Gnuplot tempVmagnetization;
	typedef std::chrono::high_resolution_clock Clock;
	const int N = 10;
	vector<vector<int>> lattice(N, vector<int>(N, 1));
	int randRow;
	int randCol;
	float energyDiff;
	float currentTemp = .001; // Kelvin
	float maxTemp = 6;
	float tempStep = .05;
	float magnetization = 0;
	const int iteraLattice = 10 * N * N;
	bool rand = 1;
	const int iteraTemp = 30;
	float magnetAvgArray[iteraTemp];
	vector<vector<float>> results;
	float magnetAvg;
	float stdDev;
	int tempIndex = 0;
	int currentIteration = 0;
	int totalIterations = (maxTemp / tempStep) * iteraLattice * iteraTemp;
	steady_clock::time_point currentTime = steady_clock::now();
	steady_clock::time_point lastTime = steady_clock::now();

	// Array element selection random distribution
	std::uniform_int_distribution<int> latticeDist(0, N - 1);

	while (currentTemp < (maxTemp + tempStep)) {
		int i = 0;
		while (i < iteraTemp)
		{
			//Reset and calculate initial average magnetization
			lattice = fillArray(lattice, N, rand, magnetization);
			// Iterate over lattice
			for (int i = 0; i < iteraLattice; i++)
			{
				// Pick random array element
				randRow = latticeDist(engine);
				randCol = latticeDist(engine);

				// Calculate energy difference (E_new - E_old)
				energyDiff = findEnergyDiff(randRow, randCol, lattice, N);
				
				// Metropolis test
				metroTest(energyDiff, lattice, randRow, randCol, magnetization, currentTemp);

				int iter = 100000;
				if (currentIteration % iter == 0 && currentIteration > 5)
				{
					reportTime(iter, currentTime, lastTime, currentIteration, totalIterations);
					lastTime = currentTime;
				}
				if (currentIteration == 1)
				{
					cout << "Calculating remaining time..." << endl;
				}
				currentIteration++;
			}
			// Normalize magnetization
			magnetization = (magnetization / (N * N));

			magnetAvgArray[i] = magnetization;
			i++;

			//plotRealTime(latticePlot, lattice);
		}

		cout << "Current temperature: " << currentTemp << endl;
		// Find average
		magnetAvg = averageArray(magnetAvgArray, iteraTemp);
		// Find standard deviation
		stdDev = findStdDev(magnetAvg, magnetAvgArray, iteraTemp);
		// Store results
		vector<float> resultsVec = {currentTemp, magnetAvg, stdDev};
		results.insert(results.end(), resultsVec);
		// Increment temperature
		currentTemp += tempStep;
		tempIndex++;
	}

	cout << "Complete... " << endl;
	generateMathematica(results);
	std::cin.get();
}