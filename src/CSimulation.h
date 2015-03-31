/*
 * File:   CSimulation.h
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 6, 2015, 11:33 AM
 */

#ifndef CSIMULATION_H
#define	CSIMULATION_H
#include <vector>
#include "CObject.h"
#include "SimulationParametersStruct.cpp"
#include "CLattice.h"
using namespace std;


class CSimulation {
public:
    CSimulation(SimulationParameters *SimPar);
    virtual ~CSimulation();
    void run(ostream& pos_stream, ostream& mag_stream, SimulationParameters *SimPar);
    void addObject(CObject* obj);
    float calcMag();
    float calcEnergy();

private:
    SimulationParameters* m_simPar;
    vector<CObject*> m_objects;

};

#endif	/* CSIMULATION_H */
