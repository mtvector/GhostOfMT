
import array
import random

import numpy
import json
from deap import algorithms
from deap import base
from deap import creator
from deap import tools

x=True
guess = [1,2,3,4,5,6]
if x:
    creator.create("FitnessMax", base.Fitness, weights=(1.0,1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    
    def initIndividual(icls, content):
        return icls(content)
    
    def initPopulation(pcls, ind_init,num_ind):
        contents = [[x + random.randint(-9,10) for x in guess] for x in range(num_ind) ]
        return pcls(ind_init(c) for c in contents)
    
    toolbox = base.Toolbox()
    
    toolbox.register("individual", initIndividual, creator.Individual)
    toolbox.register("population", initPopulation, list, toolbox.individual)
    
    #pop = toolbox.population_guess()
    
    def mutAddInt(individual, lo, hi, indpb):
        for i in range(len(individual)-1):
            individual[i] += random.randint(lo, hi)
        return individual,
    
    def evalOneMax(individual):
        return tuple([sum(individual), min(individual)])
    
    toolbox.register("evaluate", evalOneMax)
    toolbox.register("mate", tools.cxOnePoint)
    #toolbox.register("mutate", tools.mutUniformInt,low= -1,up=1, indpb=0.05)
    toolbox.register("mutate", mutAddInt,lo= -2, hi = 2, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)
    NPOP, CXPB, MUTPB, NGEN = 100,0.1, 0.4, 40
    print("Start of evolution")
    pop = toolbox.population(NPOP)
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    print zip(pop, fitnesses)
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    
    print("  Evaluated %i individuals" % len(pop))
    
    # Begin the evolution
    for g in range(NGEN):
        print("-- Generation %i --" % g)
        
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
    
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):

            # cross two individuals with probability CXPB
            if random.random() < CXPB:
                toolbox.mate(child1, child2)

                # fitness values of the children
                # must be recalculated later
                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:

            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
    
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        
        print("  Evaluated %i individuals" % len(invalid_ind))
        
        # The population is entirely replaced by the offspring
        pop[:] = offspring
        
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
        best_ind = tools.selBest(pop, 5)
        for bi in best_ind:
            print( bi + list(bi.fitness.values))
    print("-- End of (successful) evolution --")
    
    best_ind = tools.selBest(pop, 1)[0]
print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))




'''
def main():
    random.seed(64)
    
    pop = toolbox.population(300)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)
    
    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=40, 
                                   stats=stats, halloffame=hof, verbose=True)
    
    return pop, log, hof
'''
#if __name__ == "__main__":
#    main()
    
