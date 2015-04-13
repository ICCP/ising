


#PHY 480 Project 3 (Advanced Monte Carlo)
#Brandon Ewert 
#3/22/15

import random, math, time


########################################################################

def Setup(dim):
    '''Collects variables and seeds the RNG'''

    import random, math, time
    seed_num = time.time()
    random.seed(seed_num)
    setup_tup = dim
    return setup_tup

########################################################################

def Lattice_Initializer(dim,rand='no'):
    '''Initializes a lattice either uniformly spin up or randomly oriented'''

    lattice =[]
    
    #If unrandomize is chosen:
    if rand.lower() == 'no':
        
        #Create #dim rows
        for row_index in range(0,dim):
            row = []
            
            #Create #dim columns within row
            for col_index in range(0,dim):
                
                #All spins are 'up'
                row.append(int(1))
                
            lattice.append(row)
            
    elif rand.lower() == 'yes':
        for row_index in range(0,dim):
            row = []
            
            for col_index in range(0,dim):
                
                #All spins are 'random'
                spin_chooser = random.random()
                if spin_chooser < .5:
                    chosen_spin = int(1)
                elif spin_chooser >= .5:
                    chosen_spin = int(-1)
                    
                row.append(chosen_spin)
                
            lattice.append(row)
            
    return lattice

########################################################################

def Site_Chooser(dim, lattice):
    '''Chooses a random site in a lattice to begin growing a cluster'''

    #Find random site wuith RNG
    num = random.randint(0,((dim**2)-1))
    
    #Calculates Row/Col
    row = (num // dim)
    col = (num % dim)
    
    #Grabs spin
    spin = lattice[row][col]
    

    return spin, row, col

########################################################################

def Growth_Calc(temp):
    '''Calculates whether a new site will flip spin based on its
    Boltzmann Factor'''

    r = random.random()
    z =  1 - math.exp(-2/temp)
    if r < z:
        add = 'y'
    else:
        add = 'n'

    return add

########################################################################

def Cluster_NN(dim, temp, spin, row, col, lattice, queue, energy):
    '''Adds neighbors to cluster if they are of opposite sign and
        their growth_calc value is "yes". We use periodic boundary
        conditions to simulate infinite lattice.'''
    
    add_up = Growth_Calc(temp)
    if row-1 < 0:
        nn_up = [spin, dim-1, col]
    else:
        nn_up = [spin, row-1, col]

    if add_up == 'y' and lattice[row-1][col] == -spin:
        if row-1 < 0:
            lattice[dim-1][col] *= -1
        else:    
            lattice[row-1][col] *= -1 
        queue.append(nn_up)

                
    add_down = Growth_Calc(temp)
    nn_down = [spin, (row+1)%dim, col]
    
    if add_down == 'y' and lattice[(row+1)%dim][col] == -spin:
        lattice[(row+1)%dim][col] *= -1 
        queue.append(nn_down)


    add_left = Growth_Calc(temp)
    if col-1 < 0:
        nn_left = [spin, row, dim-1]
    else:
        nn_left = [spin, row, col-1]


    if add_left == 'y' and lattice[row][col-1] == -spin:
        if col-1 < 0:
            lattice[row][dim-1] *= -1
        else:
            lattice[row][col-1] *= -1 
        queue.append(nn_left)


    add_right = Growth_Calc(temp)
    nn_right = [spin, row, (col+1)%dim]
    
    if add_right == 'y' and lattice[row][(col+1)%dim] == -spin:
        lattice[row][(col+1)%dim] *= -1 
        queue.append(nn_right)
            

    energy += spin*((nn_up[0]) + (nn_down[0]) + (nn_left[0]) + (nn_right[0]))

    return queue, lattice, energy

########################################################################

def Initial_Queue(lattice, dim, temp):
    '''Builds the initial queue to start growing clusters'''
    queue = []
    energy = 0

    #Choose random site to start growing cluster
    spin, row, col = Site_Chooser(dim, lattice)
    
    #Flip spin of chosen site, add its nearest neighbors to the
    #queue
    lattice[row][col] *= -1
    spin = lattice[row][col]
   
    queue, lattice, energy = Cluster_NN(dim, temp, spin, row, col, \
                                     lattice, queue, energy)
        
    
    return queue, lattice, energy 

########################################################################

def Queue_Run(dim, temp, N):
    
    lattice = Lattice_Initializer(dim)
    for iteration in range(0,N):
    
        queue, lattice, energy = Initial_Queue(lattice, dim, temp)
        while len(queue) != 0:

            site = queue.pop(0)
            spin = site[0]
            row = site[1]
            col = site[2]

            queue, lattice, energy = Cluster_NN(dim, temp, spin, row, col, \
                                             lattice, queue, energy)
##    #Visualization (Crude)
##    for rows in lattice:
##        print()
##        for spin in rows:
##            if spin == 1:
##                print("-", end=' ')
##            else:
##                print("#", end=' ')
    
    
    mag_tot = 0
    for rows in lattice:

        for spins in rows:
            if spins == 1:
                mag_tot += 1
            else:
                mag_tot -= 1
    

    

    return mag_tot, lattice, energy 
    
    
########################################################################

def Temp_Run(dim, temp_min, temp_max, step):
    
    mag_list, en_list = [],[]
    final = int((temp_max - temp_min)/step)
    temp = temp_min
    
    for iteration in range(0,final):

        #Even number for low temp so the mag isn't negative
        if temp < 1.7:
            N = 8
        if 1.9 >= temp >= 1.7:
            N = 45
        if 2.5 > temp > 1.9:
            N = dim
        if temp >= 2.5:
            N = dim^2
        mag_tot, lattice, energy = Queue_Run(dim, temp, N)
        temp += step
        if temp < 2.5:
            mag_norm = abs(round(mag_tot/dim**2,5))
        if temp >= 2.5:
            mag_norm = round(mag_tot/dim**2,5)
        mag_list.append(mag_norm)
        
        
        en_list.append(energy)
        print('E:', energy, 'M:', mag_norm, 'T:', round(temp,3))
        
    return mag_list, en_list

########################################################################

def Temp_Run_Average(dim, temp_min, temp_max, step, ave_itera):
    
    ave_mag = [0]*int(((temp_max-temp_min)/step))
    ave_en = [0]*int(((temp_max-temp_min)/step))
    
    for iteration in range(0,ave_itera):
        
        mag_list, en_list = Temp_Run(dim, temp_min, temp_max, step)
        
        index = 0
        for s in mag_list:
            ave_mag[index] += s
            index += 1
        index = 0
        for e in en_list:
            ave_en[index] += e
            index += 1

        print()
        print('Done', ((iteration+1)/ave_itera)*100, '%')
    
        print()
        print('Ave_Mag', ave_mag)
        print()
        print('Ave_En', ave_en)
        print()
        print()

    index = 0
    for unnorm_mag in ave_mag:
        norm = (unnorm_mag) / (iteration+1)
        ave_mag[index] = round(norm,5)
        index += 1

    index = 0
    for unnorm_en in ave_en:
        norm = (unnorm_en) / (iteration+1)
        ave_en[index] = round(norm,5)
        index += 1
        
                
    print()
    print('Ave_Mag', ave_mag)
    print()
    print('Ave_En', ave_en)
    
    return ave_mag, ave_en

########################################################################

def Numerical_Derivative(a_list):

    deriv_list = []
    try:
        index = 0
        for num in a_list:
            deriv = round((a_list[index+1] - a_list[index]),5)
            deriv_list.append(abs(deriv))
            index += 1
            
    except IndexError: 
        pass

    print(deriv_list)

########################################################################
    
def Temp_Spread(dim, temp, time_i, time_f, time_step):
    
    for N in range(time_i, time_f, time_step):
        mag_norm = Temp_Run(dim, temp, N)
        print(mag_norm, ',', end = '')

########################################################################


Temp_Run_Average(100, 1, 3.5, .007, 10)

