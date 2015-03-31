/*
 * File:   CObject.h
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 6, 2015, 11:53 AM
 */

#ifndef COBJECT_H
#define	COBJECT_H
#include <fstream>
#include <streambuf>
#include <iostream>
#include "SimulationParametersStruct.cpp"
using namespace std;

class CObject {
public:
    CObject();
    virtual ~CObject();
    virtual void update(SimulationParameters* SimPar);
    virtual void print_all(ostream& pos_stream, ostream& mag_stream, SimulationParameters* SimPar);
    virtual float calcMag();
    virtual float calcEnergy();
private:

};

#endif	/* COBJECT_H */
