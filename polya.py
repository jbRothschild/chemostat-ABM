import numpy as np
from python_src import settings
from scipy.special import beta
import random
from math import comb, floor
import matplotlib.pyplot as plt
import multiprocessing

random.seed()
plt.style.use('python_src/custom.mplstyle')

import psutil
# print(psutil.cpu_count())
NUM_PROCESSES = int(np.floor( psutil.cpu_count() / 2))
psutil.cpu_times_percent(interval=1, percpu=False)


def polya_probability(vecStartAbundance, vecFinalAbundance):
    # As seen: https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution
    # for starting abundances, tells me the probability of the final
    # abundances
    
    # check to see that defined the same number of species at start and final
    if len(vecStartAbundance) != len(vecFinalAbundance):
        raise Exception("Incorrect number of species")
    
    # count number of trials that have occurred
    numDraws= np.sum(vecFinalAbundance - vecStartAbundance)
    
    # numerator of the Dirichlet-multinomial
    probability = numDraws * beta(np.sum(vecStartAbundance), numDraws)
    
    # find the count added to each species. Remove the zeros
    countAbundance = vecFinalAbundance - vecStartAbundance
    nonzeroStartAbundance = vecStartAbundance[np.where(countAbundance != 0)[0]]
    nonzeroCountAbundance = countAbundance[np.where(countAbundance != 0)[0]]
    
    # denominator of the Dirichlet-multinomial 
    for i, count in enumerate(nonzeroCountAbundance):
        probability /= count * beta(nonzeroStartAbundance[i], count)
        
    return probability


def test_polya_probability(numDraws, vecStartAbundance):
    # testing polya_probability with 2 species to check that it gives the same
    # result as polya_two_species_distribution() below
    # It seems to work!
    
    # how many individuals are there at the end of the drawing
    totalFinalAbundance = numDraws + np.sum(vecStartAbundance)
    
    # initial number of each species
    initSpecies0 = vecStartAbundance[0]
    initSpecies1 = vecStartAbundance[1]

    # the probability distirbution of the fraction of 0th species
    probabilityDistribution = np.zeros(totalFinalAbundance + 1)
    
    for k in range(numDraws + 1):
        vecFinalAbundance = vecStartAbundance \
                            + np.array([k, numDraws - k])
        probabilityDistribution[initSpecies0 + k] = \
            polya_probability(vecStartAbundance, vecFinalAbundance)
    
    return probabilityDistribution
    


def polya_two_species_distribution(numDraws, vecStartAbundance):
    # how many individuals are there at the end of the drawing
    totalFinalAbundance = numDraws + np.sum(vecStartAbundance)
    
    # initial number of each species
    initSpecies0 = vecStartAbundance[0]
    initSpecies1 = vecStartAbundance[1]

    # the probability distirbution of the fraction of 0th species
    probabilityDistribution = np.zeros(totalFinalAbundance + 1)
    
    # calculating the distribution
    print("Number of draws : ", numDraws)
    for k in range(numDraws + 1):
        # print(comb(numDraws, k))
        # print(beta(initSpecies0 + k, initSpecies1 + numDraws - k))
        # print(beta(initSpecies0, initSpecies1))
        probabilityDistribution[initSpecies0 + k] = comb(numDraws, k) * \
            beta(initSpecies0 + k, initSpecies1 + numDraws - k) / \
            beta(initSpecies0, initSpecies1)
            
    return probabilityDistribution


def polya_draws(numDraws, vecStartAbundance, vecAdvantage):
    # at each time step, one cell gets selected to duplicate. Some cells can
    # have an advantage.

    copyVecStart = np.copy(vecStartAbundance)
    for _ in range(numDraws):
        i = 0
        p_sum = 0.0
        propensity = copyVecStart * vecAdvantage
        propensity = propensity / np.sum(propensity)
        r = random.uniform(0, 1)
        
        # one way to cheeck which process happens
        while p_sum < r: # r isn't in the propensity bracket, go next species
            p_sum += propensity[i]
            i += 1
        
        copyVecStart[i - 1] += 1
    
    return copyVecStart / np.sum(copyVecStart)
    
    
def polya_trials(numTrials, numDraws, vecStartAbundance, vecAdvantage):
    # Given some initial condition, what is the distribution of species that
    # we get after a certain number of draws. Assume we're well sampling the
    # space with numTrials
    distribution = np.zeros((numTrials, len(vecStartAbundance)))
    args = [(numDraws, vecStartAbundance, vecAdvantage) for i in range(numTrials)]
    """
    for i in range(numTrials):
        distribution[i, :] = polya_draws(numDraws,
                                         vecStartAbundance,
                                         vecAdvantage)
    """
    pool = multiprocessing.Pool(NUM_PROCESSES)
    distribution = np.array(pool.starmap(polya_draws, args))
    return distribution

def testing_polya():
    initialSpecies = np.array([2, 2])
    advantage = np.array([1, 1])
    numDraws = 100
    numTrials = 100000
    
    for i in np.arange(1, 5):
        theory = polya_two_species_distribution(numDraws, np.array([i, i]))
        theoryTest = test_polya_probability(numDraws, np.array([i, i]))
        x = np.linspace(0, 1.0, len(theory))
        plt.plot(x, theory)
        plt.scatter(x, theoryTest)
        
    simula = polya_trials(numTrials,
                    numDraws,
                    initialSpecies,
                    advantage)
    weight = np.ones_like(simula[:, 0]) / len(simula[:,0])
    plt.hist(simula[:, 0], bins=numDraws+1, weights=weight, zorder=0, alpha=0.5)
    plt.xlim(0., 1.0)
    # plt.ylim(0.0, 0.5)
    
    plt.show()
    return 0

def initial_density_plot(numGenerations, chamberSize):
    coexistence = (1. / chamberSize, 1. - 1. / chamberSize)
    initDensity = [1, 5, 9, 13, 17, 21]
    probDensityLow = []
    probDensityHigh = []
    for i in initDensity:
        vecStartAbundance = np.array([i, i])
        initCells = np.sum(vecStartAbundance)
        numDraws = initCells * ( 2 ** numGenerations - 1)
        # numDraws = 256 - initCells
        #dist = polya_two_species_distribution(numDraws, vecStartAbundance)
        dist = test_polya_probability(numDraws, vecStartAbundance)      
        x = np.linspace(0, 1.0, len(dist))
        
        probDensityLow.append(np.sum(dist[x <= coexistence[0]]))
        probDensityHigh.append(np.sum(dist[ x >= coexistence[1]]))
        
        plt.plot(x, dist)
    
    plt.xlabel('final fraction')
    plt.ylabel('probability')
    plt.show()
    
    plt.plot(initDensity, probDensityLow, color=settings.colors['mCherry'])
    plt.plot(initDensity, probDensityHigh, color=settings.colors['eGFP'])
    plt.plot(initDensity,
             1. - np.array(probDensityLow) - np.array(probDensityHigh),
             color=settings.colors['coexistence'])
    plt.xlabel('initial abundance of the species')
    plt.ylabel('probability')
    plt.ylim(0.0, 1.0)
    plt.xlim(0, 21)
    plt.show()
        
    
    return 0
    
  
if __name__ == '__main__':
    # testing_polya(100000, 100, np.array([2, 2]), np.array([1, 1]))
    
    numGenerations = 12
    chamberSize = 3
    
    initial_density_plot(numGenerations, chamberSize)
    
    