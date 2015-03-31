/*
 * File:   CLattice.cpp
 * Author: Mark Wittekind and Drew Murray
 *
 * Created on March 6, 2015, 11:51 AM
 */

#include "CLattice.h"

#include <iostream>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;

CLattice::CLattice(const unsigned short h, const unsigned short w,
		const float j) :
		CObject() {
	m_height = h;
	m_width = w;
	m_COUPLING_CONSTANT = j;
	// Fill lattice with initial values
	for (unsigned short i = 0; i <= m_height; i++) {
		m_lattice.push_back(vector<bool>());
		for (unsigned short j = 0; j <= m_width; j++) {
			m_lattice[i].push_back(false);
		}
	}

}

CLattice::~CLattice() {
}

void CLattice::update(SimulationParameters* SimPar) {

	unsigned short y = rand() % m_height;
	unsigned short x = rand() % m_width;
	bool old_spin = m_lattice[y][x];
	m_lattice[y][x] = !old_spin;
	vector<unsigned short> temp = { y, x };
	m_queue.push(temp);
	while (!m_queue.empty()) {
		y = m_queue.front()[0];
		x = m_queue.front()[1];
		m_queue.pop();
		unsigned short newy, newx;
		for (short i = 0; i < 4; i++) {
			switch (i) {
			case 0:
				newy = y;
				if (x == 0)
					newx = m_width - 1;
				else
					newx = x - 1;
				break;
			case 1:
				newy = y;
				newx = (x + 1) % m_width;
				break;
			case 2:
				if (y == 0)
					newy = m_height - 1;
				else
					newy = y - 1;
				newx = x;
				break;
			case 3:
				newy = (y + 1) % m_height;
				newx = x;
				break;
			}
			if (m_lattice[newy][newx] == old_spin
					&& float(rand()) / float(RAND_MAX)
							< (1
									- exp(
											-2.0 * SimPar->COUPLING_CONSTANT
													/ (SimPar->BOLTZMAN_CONSTANT
															* SimPar->temperature)))) {
				m_lattice[newy][newx] = !old_spin;
				temp.clear();
				temp.push_back(newy);
				temp.push_back(newx);
				m_queue.push(temp);
			}
		}

	}
	return;
}

float CLattice::calcMag() {
	float s = 0;
	for (unsigned short i = 0; i < m_height; i++) {
		for (unsigned short j = 0; j < m_width; j++) {
			if (m_lattice[i][j])
				s++;
			else
				s--;
		}
	}
	s /= float(m_width * m_height);
	return abs(s);
}

float CLattice::calcEnergy() {
	signed short Energy = 0;
	for (unsigned short y = 0; y < m_height; y++) {
		for (unsigned short x = 0; x < m_width; x++) {
			for (short i = 0; i < 2; i++) {
				unsigned short newy = 0;
				unsigned short newx = 0;
				switch (i) {
					case 0:
						newy = y;
						if (x == 0)
							newx = m_width - 1;
						else
							newx = x - 1;
						break;
					case 1:
						if (y == 0)
							newy = m_height - 1;
						else
							newy = y - 1;
						newx = x;
						break;
				}
				signed short s[2] = { -1, -1 };
				if (m_lattice[newy][newx])
					s[0] = 1;
				if (m_lattice[y][x])
					s[1] = 1;
				Energy += (s[0] * s[1]);
				}
			}
		}
	return abs(2.0 * Energy * m_COUPLING_CONSTANT);
}

void CLattice::print_all(ostream& pos_stream, ostream& mag_stream, SimulationParameters* SimPar) {
	const string TRUE_CHAR = "1";
	const string FALSE_CHAR = "0";
	const bool labels = false;

	mag_stream << calcMag() << " ";

	for(unsigned short i = 1; i< m_height+1; i++) {
		string newstring = "";
		if(labels) {
			string adjust = "";

			if (i< 10)adjust = " ";
			string newstring = adjust + to_string(i);
		}
		string header("  ");
		if (i==1)
		{
			for(unsigned short j = 1; j<m_width+1; j++) {
				if(labels) {
					if(j<10)header+=" ";
					header+= to_string(j);
					header += " ";
				}
			}
			if(labels) {
				header += "\n";
				pos_stream << header;
			}
		}
		if(labels) {
			for(unsigned short j = 1; j<m_width+1; j++) {
				if(m_lattice[i][j])newstring+= " "+TRUE_CHAR+" ";
				else newstring+= " "+FALSE_CHAR+" ";
			}
		}
		else {
			for(unsigned short j = 0; j<m_width; j++) {
				if(m_lattice[i-1][j])newstring+= " "+TRUE_CHAR+" ";
				else newstring+= " "+FALSE_CHAR+" ";
			}
		}
		newstring += "\n";
		pos_stream << newstring;
	}
	pos_stream << "TIME!\n";
	return;
}

