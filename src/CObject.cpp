/*
 * File:   CObject.cpp
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 6, 2015, 11:53 AM
 */

#include "CObject.h"

CObject::CObject() {
}


CObject::~CObject() {
}

void CObject::update(SimulationParameters* SimPar){
    cout << "Called CObject update (bad)\n";
	return;
}

void CObject::print_all(ostream& pos_stream, ostream& mag_stream, SimulationParameters* SimPar){
    cout << "Called CObject print_all (bad)\n";
	return;
}
float CObject::calcMag(){
	cout << "Called CObject calcMag (bad)\n";
	return 0.0;
}
float CObject::calcEnergy(){
	cout << "Called CObject calcEnergy (bad)\n";
	return 0.0;
}

