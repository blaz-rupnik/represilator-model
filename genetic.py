import numpy as np
import matplotlib.pyplot
from scipy.integrate import odeint
import libs
import time
import multiprocessing as mp
import sys


#The GA methods are taken from https://github.com/ahmedfgad/GeneticAlgorithmPython
def select_mating_pool(pop, fitness, num_parents):
    # Selecting the best individuals in the current generation as parents for producing the offspring of the next generation.
    parents = np.empty((num_parents, pop.shape[1]))
    for parent_num in range(num_parents):
        max_fitness_idx = np.where(fitness == np.max(fitness))
        max_fitness_idx = max_fitness_idx[0][0]
        parents[parent_num, :] = pop[max_fitness_idx, :]
        fitness[max_fitness_idx] = -99999999999
    return parents


def crossover(parents, offspring_size):
    offspring = np.empty(offspring_size)
    # The point at which crossover takes place between two parents. Usually, it is at the center.
    crossover_point = np.uint8(offspring_size[1]/2)

    for k in range(offspring_size[0]):
        # Index of the first parent to mate.
        parent1_idx = k%parents.shape[0]
        # Index of the second parent to mate.
        parent2_idx = (k+1)%parents.shape[0]
        # The new offspring will have its first half of its genes taken from the first parent.
        offspring[k, 0:crossover_point] = parents[parent1_idx, 0:crossover_point]
        # The new offspring will have its second half of its genes taken from the second parent.
        offspring[k, crossover_point:] = parents[parent2_idx, crossover_point:]
    return offspring


def mutation(offspring_crossover, num_mutations=1):
    mutations_counter = np.uint8(offspring_crossover.shape[1] / num_mutations)
    # Mutation changes a number of genes as defined by the num_mutations argument. The changes are random.
    for idx in range(offspring_crossover.shape[0]):
        gene_idx = mutations_counter - 1
        for mutation_num in range(num_mutations):
            # The random value to be added to the gene.
            random_value = np.random.randint(-1, 1, 1)
            offspring_crossover[idx, gene_idx] = offspring_crossover[idx, gene_idx] + random_value
            gene_idx = gene_idx + mutations_counter
    return offspring_crossover


def calculate_fitness_for_specimen(results):
    #Currently we count the number of combinations that oscilate
    flat = np.ravel(results)
    osc_count = np.count_nonzero(flat == 1)
    return osc_count

def calculate_results_for_specimen(model, nA, nB, nC):
    new = []
    new_amp = []
    new_per = []
    logspace_vector = np.logspace(0, 5, 1000).astype(int)
    for j in range(0, len(logspace_vector)):
        params = (alpha, alpha0, beta, n, nA, nB, nC, logspace_vector[j], logspace_vector[j], logspace_vector[j])
        Z = odeint(model, Z0, t, args=params)
        A = Z[:, 3]
        is_osc, period, amplitude = libs.osc_detect(A)
        new.append(is_osc)
        new_amp.append(amplitude)
        new_per.append(period)
    return new


def calculate_results_for_generation(model, population):
    pop_res = []
    for pop in population:
        specimen_result = calculate_results_for_specimen(model, pop[0], pop[1], pop[2])
        specimen_fitness = calculate_fitness_for_specimen(specimen_result)
        print("Fitness calculated:", specimen_fitness)
        pop_res.append(specimen_fitness)
    return pop_res


def calculate_results_for_specimen_mp(params):
    model, nA, nB, nC = params
    new = []
    new_amp = []
    new_per = []
    logspace_vector = np.logspace(0, 5, 1000).astype(int)
    for j in range(0, len(logspace_vector)):
        params = (alpha, alpha0, beta, n, nA, nB, nC, logspace_vector[j], logspace_vector[j], logspace_vector[j])
        Z = odeint(model, Z0, t, args=params)
        A = Z[:, 3]
        is_osc, period, amplitude = libs.osc_detect(A)
        new.append(is_osc)
        new_amp.append(amplitude)
        new_per.append(period)
    fitness = calculate_fitness_for_specimen(new)
    return fitness


def calculate_results_for_generation_mp(model, population):
    pop_res = []
    work = []
    for pop in population:
        work.append((model, pop[0], pop[1], pop[2]))
    p = mp.Pool(num_processes)
    mp_solutions = p.map(calculate_results_for_specimen_mp, work)
    return mp_solutions

# Params
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

# Number of the weights we are looking to optimize.
num_weights = 3
num_parents_mating = 4

if __name__ == "__main__":
    # Defining the population size.
    pop_size = (20, num_weights)
    # Creating the initial population.
    new_population = np.random.randint(low=1, high=5, size=pop_size)
    print(new_population)
    best_outputs = []
    num_generations = 10

    # Parameters for multi processing
    num_processes = 8


    for generation in range(num_generations):
        print("#### STARTED generation", generation, "####")
        start = time.time()
        population_fitness = calculate_results_for_generation_mp(libs.extended_model_1, new_population)
        print("Population fitness:")
        print(population_fitness)
        end = time.time()
        print("Generation execution time:", end - start)
        best_specimen = np.max(population_fitness)
        best_outputs.append(best_specimen)
        parents = select_mating_pool(new_population, population_fitness, num_parents_mating)
        print("Parents")
        print(parents)
        offspring_crossover = crossover(parents, offspring_size=(pop_size[0] - parents.shape[0], num_weights))
        print("Crossover")
        print(offspring_crossover)
        # Adding some variations to the offspring using mutation.
        offspring_mutation = mutation(offspring_crossover, num_mutations=2)
        print("Mutation")
        print(offspring_mutation)
        # Creating the new population based on the parents and offspring.
        new_population[0:parents.shape[0], :] = parents
        new_population[parents.shape[0]:, :] = offspring_mutation
        print("New population:", new_population)
        print("#### FINISHED generation", generation, "####")

    print("Best values:", best_outputs)

