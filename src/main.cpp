#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <random>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::queue;
using std::pair;
using std::make_pair;

void addNearestNeighbors(vector< vector<int> > &lattice, queue < pair<int,int> > &index_queue, int x, int y, int dimension, double J, double kb, double T, std::default_random_engine generator){

	std::uniform_real_distribution<double> distribution(0.0,1.0);
	generator.seed(rand());

	int left_coor  = x - 1;
	int right_coor = x + 1;
	int up_coor    = y + 1;
	int down_coor  = y - 1;

	if (left_coor == -1){
		left_coor = dimension - 1;
	}
	if (right_coor == dimension){
		right_coor = 0;
	}
	if (up_coor == dimension){
		up_coor = 0;
	}
	if (down_coor == -1){
		down_coor = dimension - 1;
	}

	int current = lattice[x][y];

	int up_neighbor = lattice[x][up_coor];
	int down_neighbor = lattice[x][down_coor];
	int right_neighbor = lattice[right_coor][y];
	int left_neighbor = lattice[left_coor][y];

	double R = distribution(generator);
	if (R < (1.0 - exp(-2.0*J/kb/T)) && (current != left_neighbor)){

		index_queue.push(make_pair(left_coor, y));
		lattice[left_coor][y] *= -1;

	}

	R = distribution(generator);
	if (R < (1.0 - exp(-2.0*J/kb/T)) && (current != right_neighbor)){

		index_queue.push(make_pair(right_coor,y));
		lattice[right_coor][y] *= -1;

	}

	R = distribution(generator);
	if (R < (1.0 - exp(-2.0*J/kb/T)) && (current != up_neighbor)){

		index_queue.push(make_pair(x, up_coor));
		lattice[x][up_coor] *= -1;

	}

	R = distribution(generator);
	if (R < (1.0 - exp(-2.0*J/kb/T)) && (current != down_neighbor)){

		index_queue.push(make_pair(x, down_coor));
		lattice[x][down_coor] *= -1;

	}

}

double calculateInternalEnergy(vector< vector<int> > lattice, int dimension, double J){

	double energy = 0;

	for (int i = 0; i < dimension; i++){
		for (int j = 0; j < dimension; j++){

			int up_coor = j + 1;
			int right_coor = i + 1;

			if (up_coor == dimension){
				up_coor = 0;
			}

			if (right_coor == dimension){
				right_coor = 0;
			}

			energy += lattice[i][up_coor] * lattice[right_coor][j];

		}
	}

	return -J * energy;

}

vector<double> calculateHeatCapacity(vector<double> interalEnergyList, vector<double> temperatureList){

	vector<double> heatCapacityList;
	int limit = interalEnergyList.size();
	limit = limit - 1;

	for (int i = 0; i < limit; i++){

		double deltaE = interalEnergyList[i+1] - interalEnergyList[i];
		double deltaT = temperatureList[i+1] - temperatureList[i];

		heatCapacityList.push_back(fabs(deltaE / deltaT));
	
	}

	return heatCapacityList;

}

vector<double> calculateSusceptibility(vector<double> magnetizationList, vector<double> temperatureList){
	
	vector<double> susceptibility;
	int limit = magnetizationList.size();
	limit = limit - 1;

	for (int i = 0; i < limit; i++){

		double deltaMag = magnetizationList[i+1] - magnetizationList[i];
		double deltaTemp = temperatureList[i+1] - temperatureList[i];

		susceptibility.push_back(-1.0 * fabs(deltaMag / deltaTemp));
	}

	return susceptibility;
}

double calculateCriticalTemperature(vector<double> susceptibility, vector<double> temperatureList){

	double maxS = 0;

	double criticalTemp;

	int sLimit = susceptibility.size();

	for (int i = 0; i < sLimit; i++){

		if (fabs(susceptibility[i]) > maxS){

			maxS = fabs(susceptibility[i]);
			criticalTemp = temperatureList[i];

		}

	}

	return criticalTemp;

}

void buildCluster(vector< vector<int> > &lattice, int dimension, int x, int y, double J, double kb, double T, std::default_random_engine generator){

	// cout << "In build cluster" << endl;

	std::uniform_real_distribution<double> distribution(0.0,1.0);
	generator.seed(rand());

	queue < pair <int, int> > index_queue;
	lattice[x][y] *= -1;

	addNearestNeighbors(lattice, index_queue, x, y, dimension, J, kb, T, generator);

	while (!index_queue.empty()){

		pair <int, int> current_index = index_queue.front();
		index_queue.pop();
		int new_x = current_index.first;
		int new_y = current_index.second;

		addNearestNeighbors(lattice, index_queue, new_x, new_y, dimension, J, kb, T, generator);

	}

	/*std::ostream_iterator<int> out_it(std::cerr, " ");
	for(auto row: lattice) {
		std::copy(row.begin(), row.end(), out_it);
		std::cerr << std::endl;
	}
	cout << endl;*/

}

double findAverage(vector< vector<int> > lattice, int dimension){

	double sum = 0.0;

	for (int i = 0; i < dimension; i++){

		for (int j = 0; j < dimension; j++){

			sum += lattice[i][j];
		}
	}

	return fabs(sum / (dimension * dimension));
}

double findStdDev(vector<double> stdDevCalcList, int dimension){

	double summation = 0.0;
	double z = 0;

	for (int i = 0; i < dimension; i++){

		z += stdDevCalcList[i];
	}

	double avg = z / dimension;

	for (int i = 0; i < dimension; i++){

			summation += pow((avg - stdDevCalcList[i]), 2);
	}

	double variance = summation / (dimension * dimension);
	return sqrt(variance);
}

int main () {

	srand(time(NULL));

	int dimension = 100;
	double start_T = 0.01; // Kelvin
	double end_T = 4.0; // Kelvin
	double iter_T = 0.01; // Kelvin
	const double kb = 1.0;
	const double J = 1.0;

	std::uniform_real_distribution<double> distribution(0.0,1.0);
	std::default_random_engine generator;
	generator.seed(rand());

	vector<double> magnetizationList;
	vector<double> temperatureList;
	vector<double> stdDevCalcList;
	vector<double> stdDevList;
	vector<double> interalEnergyList;

	/*std::ostream_iterator<int> out_it(std::cerr, " ");
	for(auto row: lattice) {
		std::copy(row.begin(), row.end(), out_it);
		std::cerr << std::endl;
	}*/

	double magnetization;
	
	int clusterN = 10;
	int stdDevRuns = 20;

	// Iterate through every temperature value
	for (double T = start_T; T <= end_T; T += iter_T){

		magnetization = 0;
		vector< vector<int> > lattice(dimension, vector<int> (dimension,1));
		
		// Uncomment below to randomize initial lattice

		/*for (int i = 0; i < dimension; i++){
			for (int j = 0; j < dimension; j++){

				double rand_float = distribution(generator);

				if (rand_float < 0.5){

					lattice[i][j] = -1;

				} else{

					lattice[i][j] = 1;

				}
			}
		}*/

		if (T > 2.0 && T <= 2.5){
			clusterN = 2 * dimension;
		}
		if (T > 2.5){
			clusterN = 2 * dimension * dimension;
		}

		double totMag = 0.0;
	
		for (int k = 0; k < stdDevRuns; k++){

			for (int i = 0; i <= clusterN; i++){

				int i_rand = rand() % dimension;
				int j_rand = rand() % dimension;
				//int i_rand = 0;
				//int j_rand = 0;
				buildCluster(lattice, dimension, i_rand, j_rand, J, kb, T, generator);


			} // End of Metropolis Test loop

			magnetization = findAverage(lattice, dimension);
			totMag += magnetization;
			stdDevCalcList.push_back(magnetization);

		}

		double avgMag = totMag / stdDevRuns;
		magnetizationList.push_back(avgMag);
		temperatureList.push_back(T);

		double stdDev = findStdDev(stdDevCalcList, stdDevCalcList.size());
		stdDevList.push_back(stdDev);

		double internalEnergy = calculateInternalEnergy(lattice, dimension, J);
		interalEnergyList.push_back(internalEnergy);

		cout << "For temperature " << T << ", magnetization: " << avgMag << ", Std Dev: " << stdDev << endl;

		// output << T << "\t" << "\t" << avgMag << "\t" << "\t" << stdDev << endl;

		
	} // End of temperature iteration loop

	ofstream output;
	output.open("Output.txt");
	output << "Temperature" << "\t" << "Magnetization" << "\t" << "Standard Deviation" << "\t" << "Internal Energy" << "\t" << "Susceptibility" << "\t" << "Heat Capacity" << endl;

	vector<double> susceptibility = calculateSusceptibility(magnetizationList, temperatureList);
	vector<double> heatCapacityList = calculateHeatCapacity(interalEnergyList, temperatureList);

	int arraySize = magnetizationList.size();

	for (int i = 0; i < arraySize; i++){

		if (i == arraySize - 1){

			output << temperatureList[i] << "\t" << "\t" << magnetizationList[i] << "\t" << "\t" << stdDevList[i] << "\t" << "\t" << interalEnergyList[i] << endl;

		}else{

			output << temperatureList[i] << "\t" << "\t" << magnetizationList[i] << "\t" << "\t" << stdDevList[i] << "\t" << "\t" << interalEnergyList[i] << "\t" << "\t" << susceptibility[i] << "\t" << "\t" << heatCapacityList[i] << endl;

		}
	}

	double criticalTemp = calculateCriticalTemperature(susceptibility, temperatureList);

	// Gamma calculation loop below

	double gamma;
	int gammaLimit = magnetizationList.size();
	double Tc = 2.21;
	for (int i = 0; i < gammaLimit; i++){

		double tempDiff = fabs(temperatureList[i] - Tc);
		gamma = log(magnetizationList[i]) / log(tempDiff);

		cout << gamma << endl;

	}

	cout << "Critical Temperature: " << criticalTemp << endl;

	output.close();
	cout << "Completed." << endl;
	cout << "Output.txt file created" << endl;
	return 0;

}