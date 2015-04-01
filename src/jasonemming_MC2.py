#Jason Emming
#Cluster Monte Carlo Method: Ising Model

from __future__ import division
import numpy as np
import random
import math
import matplotlib.pyplot as plt

global J, kb, N 	#Global Variables defined


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


J=1
kb=1 				#Constants are set to 1 (like any good physist)

graph = 1

start_T = 2.1
stop_T = 2.4 		#Stoping temperature
diff_T = 0.01 		#Step size of temperature
Temp = start_T		#Starting temperature
N = 10				#Size of lattice
repeat = 10			#Repeats run for each temperature Y times, so the end magnitization can be averaged


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#Initializes a 'list of lists' for the 2D lattice structure
#Fills each cell with spin up (1) or random spins
def initialize(N):
	lattice = []
	for i in range(N):
		sub_list = []
		for j in range(N):
			sub_list.append(1)
		lattice.append(sub_list)
	return lattice


#Prints out the NxN lattice of spin states
def print_lattice(lattice):
	for i in range(N):
		out = []
		for j in range(N):
			out.append(lattice[i][j])
		print out


#Returns two random numbers from 0 to N-1
def random_coord(N):
	x = random.randint(1, 1000)%N
	y = random.randint(1, 1000)%N	
	return x,y


#Returns bool on whether or not to grow cluster
def yay_or_nay(T):
	rand = random.random()
	expo = 1 - math.exp(-2.0*J/(kb*T))
	if expo > rand:
		yay = True
	else:
		yay = False
	return yay


#Checks if the neighbors can and should be added to the queue given the temperature
def neighbors(lattice,i,j,Temp,queue):
	if lattice[(i-1)%N][j] != lattice[i][j]:
		if yay_or_nay(Temp):
			queue.append( [(i-1)%N,j] )
			lattice[(i-1)%N][j] *= -1
	if lattice[(i+1)%N][j] != lattice[i][j]:
		if yay_or_nay(Temp):
			queue.append( [(i+1)%N,j] )
			lattice[(i+1)%N][j] *= -1
	if lattice[i][(j-1)%N] != lattice[i][j]:
		if yay_or_nay(Temp):
			queue.append( [i,(j-1)%N] )
			lattice[i][(j-1)%N] *= -1
	if lattice[i][(j+1)%N] != lattice[i][j]:
		if yay_or_nay(Temp):
			queue.append( [i,(j+1)%N] )
			lattice[i][(j+1)%N] *= -1


#Returns the TOTAL magnitization. Not average.
def magnitization(lattice):
	total = 0
	N1 = len(lattice)
	for i in range(N):
		N2 = len(lattice[i])
		for j in range(N2):
			total += lattice[i][j]		
	return total


#Returns the average magnitization for an individual temperature
def avg_mag_temp(last_mag):
	last_mag = np.array(last_mag)
	mean = np.mean(last_mag)/(N*N)
	stdev = np.std(last_mag)/(N*N)
	return mean, stdev


#Used for the determination of the change in magetization
def dmdt(last_mag,temp_list):
	out =[]
	for i in range(len(last_mag) - 1):
		dm = last_mag[i+1] - last_mag[i]
		dt = temp_list[i+1] - temp_list[i]
		out.append(-abs(dm/dt))
	return out


#Used for the determination of the change in internal energy
def dEdt(last_mag,temp_list):
	out =[]
	for i in range(len(last_mag) - 1):
		dm = last_mag[i+1] - last_mag[i]
		dt = temp_list[i+1] - temp_list[i]
		out.append(dm/dt)
	return out


#Determines the interal energy given a lattice of spins
def energy(lattice):
	sigma = 0
	for i in range(len(lattice)):
		for j in range(len(lattice[i])):
			sigma += lattice[i][j]*( lattice[i][(j+1)%N] + lattice[(i+1)%N][j] )
	out = sigma*(-J)
	return out


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


mag_temp = []
mag_temp_std = []

dmdt_temp = []

avg_energy = []
energy_std = []

temp_list = []

#Starts a new series of runs for a range of temperatures
while Temp < stop_T:
	Temp += diff_T

	if Temp < 2:				#Grows X amount of clusters for each temperature
		cluster = 10
	elif Temp < 2.3:
		cluster = 2*N
	elif Temp < 2.5:
		cluster = N*N//2
	elif 2.55< Temp < 2.65:
		cluster = 3*N*N
	else: 
		cluster = 2*N*N

	temp_list.append(Temp)
	print "Running Temperature: ", Temp

	last_mag = []				#Holds the last magetization for every run
	int_energy = []				#Holds the internal energy of all runs of a specific temperature

	run = 0
	while run < repeat:
		run += 1
		print "RUN: ", run

		mag = []				#Holds the magnitization after each cluster finishes
		queue = []				#Holds the members of each cluster that have neighbors to be checked

		lattice = initialize(N)

		#Runs through growing a cluster
		for num in range(1,cluster+1):
			r_i,r_j = random_coord(N)					#Chooses random coordinates for first spot to flip
			lattice[r_i][r_j] = lattice[r_i][r_j]*(-1)	#Flips the spin of location
			
			neighbors(lattice,r_i,r_j,Temp,queue)		#Neighbors with opposite spin are probabilistically added to queue with their spin flipped

			#Runs until cluster cannot expand any further
			while len(queue) != 0:
				new_i, new_j = queue[0][0], queue[0][1] #New spot chosen from next in queue
				queue.pop(0)							#Removes spot out of queue

				neighbors(lattice, new_i, new_j, Temp, queue)	#Runs same neighbor function

		int_energy.append(energy(lattice))				#Adds the internal energy of the lattice of X amount of cluster runs
		last_mag.append(abs(magnitization(lattice)))	#Adds the final magetization to a list after X amount of cluster runs
	
	avg_energy_mean, avg_energy_std = avg_mag_temp(int_energy)
	avg_energy.append(avg_energy_mean)					#Adds the averaged value of internal energy for X amount of runs
	energy_std.append(avg_energy_std)
	
	avg_mag_mean, avg_mag_std = avg_mag_temp(last_mag)
	mag_temp_std.append(avg_mag_std)
	mag_temp.append(avg_mag_mean)						#Adds the averaged value of magetization for each temperature


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# Open a file
fo = open("data1.txt", "wb")
fo.write('\n'+"Temperature"+'\n')
for t in range(len(temp_list)):
	fo.write(str(temp_list[t])+'\t')


fo.write('\n'+"Magnetization"+'\n')
for m in range(len(mag_temp)):
	fo.write(str(mag_temp[m])+'\t')

fo.write('\n'+"Magnetization STD"+'\n')
for s in range(len(mag_temp_std)):
	fo.write(str(mag_temp_std[s])+'\t')

dmdt_list = dmdt(mag_temp,temp_list)
fo.write('\n'+"dMdT"+'\n')
for d in range(len(dmdt_list)):
	fo.write(str(dmdt_list[d])+'\t')


fo.write('\n'+"Energy"+'\n')
for e in range(len(avg_energy)):
	fo.write(str(avg_energy[e])+'\t')

fo.write('\n'+"Energy STD"+'\n')
for T in range(len(energy_std)):
	fo.write(str(energy_std[T])+'\t')


dEdt_list = dEdt(avg_energy,temp_list)
fo.write('\n'+"dEdT"+'\n')
for D in range(len(dEdt_list)):
	fo.write(str(dEdt_list[D])+'\t')


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


plt.figure()

if graph == 1:
	plt.title("Magnetization vs Temperature")
	plt.scatter(temp_list,mag_temp)
	plt.errorbar(temp_list,mag_temp,yerr=mag_temp_std, linestyle="None")
	plt.ylim([-0.5,1.1])
	plt.xlabel("Temperature")
	plt.ylabel("Magnetization")
	plt.axhline(y=0,c='k')
	plt.show()	

	plt.title("Change in Magnetization vs Temperature")
	plt.scatter(temp_list[:-1], dmdt_list)
	plt.xlabel("Temperature")
	plt.ylabel("dM/dT")
	plt.axhline(y=0,c='k')
	plt.show()

	plt.title("Internal Energy vs Temperature")
	plt.scatter(temp_list,avg_energy)
	plt.errorbar(temp_list,avg_energy,yerr=avg_energy_std, linestyle="None")
	plt.xlabel("Temperature")
	plt.ylabel("Internal Energy")
	plt.axhline(y=0,c='k')
	plt.show()

	plt.title("Specific Heat vs Temperature")
	plt.scatter(temp_list[:-1],dEdt_list)
	plt.xlabel("Temperature")
	plt.ylabel("Specific Heat (dE/dT)")
	plt.axhline(y=0,c='k')
	plt.show()	

	Tc = 2.28					#Best estimate received from the peak of the susceptibility graph (dMdT v T)
	plt.title("ln(m) vs ln(|T-Tc|)")
	plt.scatter(np.log(np.abs(np.array(temp_list) - Tc)) , np.log(np.array(mag_temp)))
	plt.xlabel("ln(|T-Tc|)")
	plt.ylabel("ln(magetization)")
	plt.xlim([-5,0])	
	plt.axhline(y=0,c='k')
	plt.show()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
