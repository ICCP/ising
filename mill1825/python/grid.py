import math
from math import exp
import random
import numpy as np

random.seed()

def pause():
    garbage = raw_input("[PAUSE]")

class ising_grid(object):    
    def __init__(self, x_dim=10, y_dim=10):
        self.grid = []
        self.x_max = x_dim
        self.y_max = y_dim
        self.J = 1.0
        self.B = 1.0
        self.T = 1.0
        self.beta = 1/self.T   

    def initialize(self):
        self.grid = np.ones([self.y_max, self.x_max])

    def compute_boltzmann_factor(self, S):
        return (exp(2*(self.J*S + self.B)*self.beta), exp(-2*(self.J*S + self.B)*self.beta))

    def metropolis(self):

        self.initialize()

        boltzmann_factor = {}
        for S in [-4, -2, 0, 2, 4]: 
            boltzmann_factor[S] = self.compute_boltzmann_factor(S)

        M = []

        for k in range(self.xM*self.y_max*self.x_max):

            # ------------------------
            # Metropolis algorithm 
            # ------------------------

            i = random.randint(0,self.x_max-1)
            j = random.randint(0,self.y_max-1)
            val = self.spin(i,j)
            
            # if spin i is down, choose first column in precomputed boltzmann factors (corresponds to flipping the spin)
            if(val == -1):
                w = boltzmann_factor[self.get_S(i, j)][0]
            # if spin i is down, choose second column in precomputed boltzmann factors (corresponds to flipping the spin)
            elif(val == 1):
                w = boltzmann_factor[self.get_S(i, j)][1]
            if(w > 1):
                self.flip(i, j)
            else:
                if(w > random.random()):
                    self.flip(i, j)

            M.append(self.get_M())

        M = np.array(M)

        return np.mean(M), self.std_block(M)

    def std_block(self, data):
        std = []
        blocks = 50
        split = len(data)/blocks 
        for i in range(blocks):
            std.append(np.std(data[(i*split):((i+1)*split)]))
        return np.std(np.array(std))

    def set_T(self,T):
        self.T = float(T)
        self.beta = 1./self.T

    def spin(self, i, j):
        return self.grid[j][i]    

    def get_M(self):
        return np.sum(self.grid)/self.x_max/self.y_max

    def flip(self, i, j):
        self.grid[j][i] *= -1

    def get_E(self):
        E = 0.0
        for j in range(self.y_max):
            for i in range(self.x_max):
                E += self.get_S(i, j) 
        return E

    def get_S(self, i_in, j_in, p=False):
        i, j = i_in%self.x_max, j_in%self.y_max
        # This is where the periodic BCs are implemented
        ip1, jp1 = (i_in+1)%self.x_max, (j_in+1)%self.y_max
        im1, jm1 = (i_in-1)%self.x_max, (j_in-1)%self.y_max
        return (self.grid[jm1][i] + self.grid[jp1][i] + self.grid[j][ip1] + self.grid[j][im1])

