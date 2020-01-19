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
def signal_oscilates(a):
    peaks = signal.find_peaks(a)
    num_of_peaks = len(peaks[0])

    # Worst case scenario, model not even close to oscilating
    if num_of_peaks <= 3:
        return (False, 0, 0)

    peak_values = a[peaks[0]]

    peak_difference = abs(peak_values[-3] - peak_values[-2])

    if (abs(peak_values[-3] - peak_values[-2])) <= 0.001:
        peak_values = a[peaks[0]]
        return (True, 1, peak_values[-2])
    else:
        return (False, 1 - peak_difference, 0)
    

def calculate_gene_fitness(input):
    params = (alpha, alpha0, beta, n, input[0], input[1], input[2], input[3], input[4], input[5])
    result = odeint(extended_model_1, Z0, t, args=params)

    A = result[:,3]

    return signal_oscilates(result[:,3])

def select_best_parents(population, fitness, num_parents_mating, generation):
    indices = (-fitness).argsort()[:num_parents_mating]
    print("Best for generation: " + str(generation) + ", Fitness: " + str(fitness[indices[0]]))
    print(population[indices[0],:])
    return population[indices,:]

def perform_crossover(parents, offspring_size,num_parents):
    offspring = np.empty(offspring_size, dtype=int)
    for k in range(offspring_size[0]):
        parent_1 = parents[k % num_parents,:]
        parent_2 = parents[(k+1) % num_parents,:]
        offspring[k,:] = generate_from_parents(parent_1,parent_2)     
    return offspring
    

def generate_from_parents(parent_1,parent_2):
    result = np.empty(6,dtype=int)
    for i in range(6):
        parent_to_choose = np.random.randint(0,2)
        if parent_to_choose == 0:
            result[i] = parent_1[i]
        else:
            result[i] = parent_2[i]
    return result

def apply_mutation(offspring,num_of_offspring):
    #we change single weight randomly for each element in offspring
    for i in range(num_of_offspring):
        position_to_change = np.random.randint(0,6)
        if position_to_change < 3:
            offspring[i,position_to_change] = np.random.randint(1,5)
        else:
            offspring[i,position_to_change] = np.random.randint(0,100000)
    return offspring

def append_random(offspring,num_of_random):
    size_to_add = (num_of_random,3)

    paramN = np.random.randint(1,5,size=size_to_add)
    paramK = np.random.randint(1,100000,size=size_to_add)

    randomParam = np.concatenate((paramN,paramK), axis=1)
    return np.concatenate((offspring,randomParam), axis=0)

def show_solution(nA,nB,nC,kA,kB,kC):
    #testing outcome
    params = (alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC)
    result = odeint(extended_model_1, Z0, t, args=params)

    A = result[:,3]
    B = result[:,4]
    C = result[:,5]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t,A,'g:',label='A(t)')
    ax.plot(t,B,'b-',label='B(t)')
    ax.plot(t,C,'r--',label='C(t)')
    ax.set_ylabel('concentrations')
    ax.set_xlabel('time')
    ax.legend(loc='best')
    plt.show()


#Number of weight we are looking to optimize
num_weights = 6

#How many solution will one population have
solutions_per_population = 16

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

# Parameters
num_generations = 200
num_parents_mating = 8
perform_mutation = False
perform_genetic = False

if perform_genetic:
    for generation in range(num_generations):
        new_fitness = []
        for i in range(solutions_per_population):
            (oscilates, fitness_to_add, amplitude) = calculate_gene_fitness(initial_population[i,:])
            if oscilates:
                print("Found oscilating solution for parameters: " + str(initial_population[i,:]))
                print("Generation number: " + str(generation) + " Amplitude: " + str(amplitude))
            new_fitness.append(fitness_to_add)

        #select best parents based on fitness
        parents = select_best_parents(initial_population, np.array(new_fitness), num_parents_mating, generation)

        #perform crossover
        offspring = perform_crossover(parents,(num_parents_mating,6),num_parents_mating)

        #perform mutation
        if perform_mutation:
            offspring = apply_mutation(offspring,num_parents_mating)
        
        #append random genes 
        initial_population = append_random(offspring,solutions_per_population - num_parents_mating)


#show_solution(2, 1, 3, 33852, 548, 77207)
#solution that was found shown graphical