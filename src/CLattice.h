/*
 * File:   CLattice.h
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 6, 2015, 11:51 AM
 */

#ifndef CLATTICE_H
#define	CLATTICE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include "SimulationParametersStruct.cpp"

#include "CObject.h"

using namespace std;

class CLattice : public CObject{
public:
    CLattice(const unsigned short h, const unsigned short w, const float j);
    ~CLattice();
    void print_all(ostream& pos_stream, ostream& mag_stream, SimulationParameters* SimPar);
    void update(SimulationParameters* SimPar);
    float calcMag();
    float calcEnergy();
private:
    vector<vector<bool> > m_lattice;
    unsigned short m_height;
    unsigned short m_width;
    float m_COUPLING_CONSTANT;
    queue<vector<unsigned short>> m_queue;

};

#endif	/* CLATTICE_H */
