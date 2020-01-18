import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import signal
import time
import ga

def extended_model_1(Z, t, alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = alpha / (1 + C ** n + (A / kA) ** nA) + alpha0 - a
    dbdt = alpha / (1 + A ** n + (B / kB) ** nB) + alpha0 - b
    dcdt = alpha / (1 + B ** n + (C / kC) ** nC) + alpha0 - c

    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)

    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]

    return dzdt


#Initial parameters
Z0 = [ 7.02028482, 1.57288111, 58.50876295, 8.69333201, 1.40155783, 53.49915358]
# number of time points
t_end = 100
n = t_end * 10
# time points
t = np.linspace(0, t_end, n)
alpha = 216
alpha0 = 0.001 * alpha
n = 2
beta = 5

#example of making model
#params_1 = (alpha, alpha0, beta, n, 2, 2, 2, 100, 100, 100)
#Z = odeint(extended_model_1, Z0, t, args=params_1)

#A = Z[:,3]
#B = Z[:,4]
#C = Z[:,5]


############################################
# Helper functions for genetic algorithm
############################################
def calculate_pop_fitness():
    #TODO
    return 1

#Number of weight we are looking to optimize
num_weights = 6

#How many solution will one population have
solutions_per_population = 8

pop_size = (solutions_per_population,3)

############################
#Creating initial population
############################

#for nA, nB, nC
inital_1 = np.random.randint(1,5,size=pop_size)

#for kA, kB, kC
inital_2 = np.random.randint(1,100000,size=pop_size)

#8 chromosomes with 6 genes, first three for nX, last three for kX
initial_population = np.concatenate((inital_1,inital_2), axis=1)

###############################
# Generation and parents mating
###############################

num_generations = 5
num_parents_mating = 4

for generation in range(num_generations)
    #Measuring the fitness of each chromosome in the population

