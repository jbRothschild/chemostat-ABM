import numpy as np
from python_src import settings
from scipy.special import beta
import random
from math import comb, floor
import matplotlib.pyplot as plt
import multiprocessing
import matplotlib
import os

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import Normalize #color in log the maps

random.seed()
plt.style.use('python_src/custom.mplstyle')

import psutil
NUM_PROCESSES = int(np.floor( psutil.cpu_count() / 2))
psutil.cpu_times_percent(interval=1, percpu=False)

FIGDIR = 'polya'

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


def polya_probability_two_species(numDraws, vecStartAbundance):
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
    # checking whether 2 methods are identical (they are)
    initialSpecies = np.array([2, 2])
    advantage = np.array([1, 1])
    numDraws = 100
    numTrials = 100000
    
    for i in np.arange(1, 5):
        theory = polya_two_species_distribution(numDraws, np.array([i, i]))
        theoryTest = polya_probability_two_species(numDraws, np.array([i, i]))
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


def initial_density_plot_A(numGenerations, fixateFraction):
    # define regimes of fixation of coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    
    # density arrays for results
    initDensity = [1, 5, 9, 13, 17, 21]
    probDensityLow = []
    probDensityHigh = []
    cmapname = 'plasma'
    cmap = matplotlib.cm.get_cmap(cmapname)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(initDensity))]
    
    for i, density in enumerate(initDensity):
        vecStartAbundance = np.array([density, density])
        initCells = np.sum(vecStartAbundance)
        # numDraws = initCells * ( 2 ** numGenerations - 1)
        numDraws = numGenerations - initCells
        #dist = polya_two_species_distribution(numDraws, vecStartAbundance)
        dist = polya_probability_two_species(numDraws, vecStartAbundance)      
        x = np.linspace(0, 1.0, len(dist))
        
        plt.plot(x, dist, color=colors[i], label=f"{density}")
    plt.xlabel('Final fraction eGFP')
    plt.ylabel('Probability')
    #plt.axvline(x = coexistence[0], color=settings.colors['mCherry'], ls=':')
    #plt.axvline(x = coexistence[1], color=settings.colors['eGFP'], ls=':')
    #plt.text(0.4, 0.00011, 'Coexistence', color=settings.colors['coexistence'])
    #plt.text(0.025, 0.00011, 'mCherry fix.', color=settings.colors['mCherry'])
    #plt.text(0.825, 0.00011, 'eGFP fix.', color=settings.colors['eGFP'])
    plt.legend(title='Initial density', loc=(0.01, 0.1), framealpha=0.9)
    # ax1.set_ylim(0.0, 1.0)
    plt.xlim(0.0, 1.0)
    
    filename = 'polya_density_V1'
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
        
    return 0
    

def initial_density_plot(numGenerations, fixateFraction):
    # define regimes of fixation of coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    
    # density arrays for results
    initDensity = [1, 5, 9, 13, 17, 21]
    probDensityLow = []
    probDensityHigh = []
    cmapname = 'plasma'
    cmap = matplotlib.cm.get_cmap(cmapname)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(initDensity))]
    
    # define figure
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    
    for i, density in enumerate(initDensity):
        vecStartAbundance = np.array([density, density])
        initCells = np.sum(vecStartAbundance)
        numDraws = numGenerations - initCells
        #dist = polya_two_species_distribution(numDraws, vecStartAbundance)
        dist = polya_probability_two_species(numDraws, vecStartAbundance)      
        x = np.linspace(0, 1.0, len(dist))
        
        probDensityLow.append(np.sum(dist[x <= coexistence[0]]))
        probDensityHigh.append(np.sum(dist[ x >= coexistence[1]]))
        
        ax1.plot(x, dist, color=colors[i], label=f"{density}")
    ax1.set_xlabel('Final fraction eGFP')
    ax1.set_ylabel('Probability')
    ax1.axvline(x = coexistence[0], color=settings.colors['mCherry'], ls=':')
    ax1.axvline(x = coexistence[1], color=settings.colors['eGFP'], ls=':')
    ax1.text(0.4, 0.00011, 'Coexistence', color=settings.colors['coexistence'])
    ax1.text(0.025, 0.00011, 'mCherry fix.', color=settings.colors['mCherry'])
    ax1.text(0.825, 0.00011, 'eGFP fix.', color=settings.colors['eGFP'])
    ax1.legend(title='Initial density', loc=(0.01, 0.1), framealpha=0.9)
    # ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.0, 1.0)
    
    ax2.plot(initDensity, probDensityLow, color=settings.colors['mCherry'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['mCherry'],
             label=settings.labels['mCherry'])
    ax2.plot(initDensity, probDensityHigh, color=settings.colors['eGFP'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['eGFP'],
             label=settings.labels['eGFP'])
    ax2.plot(initDensity,
             1. - np.array(probDensityLow) - np.array(probDensityHigh),
             color=settings.colors['coexistence'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['coexistence'],
             label=settings.labels['coexistence'])
    ax2.set_xlabel('Initial abundance of the species')
    ax2.set_ylabel('Probability')
    ax2.legend()
    ax2.set_ylim(0.0, 1.0)
    ax2.set_xlim(0, 21)
    
    filename = 'polya_density'
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
        
    return 0


def initial_density_growth_plot(numTrials, numGenerations, fixateFraction,
                                advantages):
    # binning simulations
    nbins = 21
    bins = np.linspace(0.0, 1.0, nbins)
    center = (bins[:-1] + bins[1:]) / 2
    
    # defining regions of fixation and coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    maskL = center <= coexistence[0]
    maskH = center >= coexistence[1]
    maskC = (center > coexistence[0]) & (center < coexistence[1])
    
    # density arrays for results
    initDensity = [1, 5, 9, 13, 17, 21]
    probDensityLow = []
    probDensityHigh = []
    probDensityCoexist = []
    cmapname = 'plasma'
    cmap = matplotlib.cm.get_cmap(cmapname)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(initDensity))]
    
    # define figure
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    vecAdvantage = np.array(advantages)
    for i, density in enumerate(initDensity):
        # settings for simualtion
        vecStartAbundance = np.array([density, density])
        initCells = np.sum(vecStartAbundance)
        numDraws = numGenerations - initCells
        
        # simulation
        simula = polya_trials(numTrials, numDraws, vecStartAbundance,
                              vecAdvantage)

        # turning simulation results into a probability distribution
        weight = np.ones_like(simula[:, 0]) / len(simula[:,0])
        hist, _ = np.histogram(np.array(simula[:, 0]), bins, weights = weight)

        # appending to results arrays
        probDensityLow.append(np.sum(hist[maskL]))
        probDensityHigh.append(np.sum(hist[maskH]))
        probDensityCoexist.append(np.sum(hist[maskC]))
        
        # plotting distributions for each fitness advantage
        ax1.plot(center, hist, color=colors[i], label=f"{density}")
    
    # plot settings for probability of fractions
    ax1.set_xlabel('Final fraction eGFP')
    ax1.set_ylabel('Probability')
    ax1.axvline(x = coexistence[0], color=settings.colors['mCherry'], ls=':')
    ax1.axvline(x = coexistence[1], color=settings.colors['eGFP'], ls=':')
    ax1.text(0.4, 0.27, 'Coexistence', color='k')
    ax1.text(0.025, 0.27, 'mCherry', color='k')
    ax1.text(0.85, 0.27, 'eGFP', color='k')
    ax1.legend(title='Initial density', loc=(0.03, 0.1), framealpha=0.9)
    ax1.set_ylim(0.0, 0.3)
    ax1.set_xlim(0.0, 1.0)
    
    # plot for probability of certain densities
    ax2.plot(initDensity, probDensityLow, color=settings.colors['mCherry'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['mCherry'],
             label=settings.labels['mCherry'])
    ax2.plot(initDensity, probDensityHigh, color=settings.colors['eGFP'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['eGFP'],
             label=settings.labels['eGFP'])
    ax2.plot(initDensity,
             1. - np.array(probDensityLow) - np.array(probDensityHigh),
             color=settings.colors['coexistence'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['coexistence'],
             label=settings.labels['coexistence'])
    
    # plot settings
    ax2.set_xlabel('Initial abundance of the species')
    ax2.set_ylabel('Probability')
    ax2.legend()
    ax2.set_ylim(0.0, 1.0)
    ax2.set_xlim(0, 21)
    
    filename = 'polya_density_g_5'
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
        
    return 0


def selection_advantage_plot(numTrials, numGenerations, vecStartAbundance,
                             advantage, fixateFraction):
    # binning simulations
    nbins = 21
    bins = np.linspace(0.0, 1.0, nbins)
    center = (bins[:-1] + bins[1:]) / 2
    
    # defining regions of fixation and coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    maskL = center <= coexistence[0]
    maskH = center >= coexistence[1]
    maskC = (center > coexistence[0]) & (center < coexistence[1])
    
    # arrays for plotting results and coloring scheme
    probDensityLow = []
    probDensityHigh = []
    probDensityCoexist = []
    cmapname = 'viridis'
    cmap = matplotlib.cm.get_cmap(cmapname)
    range_cmap = (0.0, 0.8)
    colors = [cmap(nbr) for nbr in np.linspace(range_cmap[0], range_cmap[1],
                                               num=len(advantage))]
    short_range = cmap(np.linspace(range_cmap[0], range_cmap[1], cmap.N))
    cmap_short = matplotlib.colors.LinearSegmentedColormap.from_list('cut_map', short_range)
    
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    for i, fitness in enumerate(advantage):
        # simualtions
        vecAdvantage = np.array([fitness, 1.0])
        initCells = np.sum(vecStartAbundance)
        numDraws = numGenerations - initCells
        simula = polya_trials(numTrials, numDraws, vecStartAbundance,
                              vecAdvantage)

        # turning simulation results into a probability distribution
        weight = np.ones_like(simula[:, 0]) / len(simula[:,0])
        hist, _ = np.histogram(np.array(simula[:, 0]), bins, weights = weight)

        # appending to results arrays
        probDensityLow.append(np.sum(hist[maskL]))
        probDensityHigh.append(np.sum(hist[maskH]))
        probDensityCoexist.append(np.sum(hist[maskC]))
        
        # plotting distributions for each fitness advantage
        ax1.plot(center, hist, color=colors[i], label=f"{fitness:.2f}")
    ax1.set_xlabel('Final fraction eGFP')
    ax1.set_ylabel('Probability')
    ax1.axvline(x = coexistence[0], color=settings.colors['mCherry'], ls=':')
    ax1.axvline(x = coexistence[1], color=settings.colors['eGFP'], ls=':')
    ax1.text(0.4, 0.9, 'Coexistence', color='k')
    ax1.text(0.025, 0.9, 'mCherry', color='k')
    ax1.text(0.85, 0.9, 'eGFP', color='k')
    bbox = (0.35, 0.3, 0.3, 0.03)
    cbaxes = inset_axes(ax1, width='100%', height='100%', bbox_to_anchor=bbox,
                        bbox_transform=ax1.transAxes,
                        loc='center') 
    norm = Normalize(1.0, 2.0)
    #cbaxes = inset_axes(ax1, width="30%", height="3%", loc=cloc, bbox_to_anchor=bbox)  
    cax = f.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_short),
                     cax=cbaxes, orientation='horizontal',
                     label=r'Growth Rate Ratio')
    cbaxes.xaxis.set_label_position('top')
    cax.ax.tick_params(size=0)
    # ax1.legend(title='Growth rate ratio', loc=(0.25, 0.1), framealpha=0.1)
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.0, 1.0)
    
    # plotting fixations and coexistence fractions for each fitness ratio 
    ax2.plot(advantage, probDensityLow, color=settings.colors['mCherry'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['mCherry'],
             label=settings.labels['mCherry'])
    ax2.plot(advantage, probDensityHigh, color=settings.colors['eGFP'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['eGFP'],
             label=settings.labels['eGFP'])
    ax2.plot(advantage, probDensityCoexist,
             color=settings.colors['coexistence'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['coexistence'],
             label=settings.labels['coexistence'])
    ax2.legend()
    ax2.set_xlabel('Growth rate ratio (eGFP/mCherry)')
    ax2.set_ylabel('Probability')
    ax2.set_ylim(0.0, 1.0)
    # ax2.set_xlim(1.0, 2.0)
    
    filename = 'polya_growth'
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
    
    return 0


def multispecies_coexistence(nbrCellsEnd, fixateFraction):
    # define regimes of fixation of coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    
    # density arrays for results
    initDensity = np.arange(1, 22, 4)
    averageRichness = []
    cmapname = 'viridis'
    cmap = matplotlib.cm.get_cmap(cmapname)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(initDensity))]
    
    # define figure
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
    
    for i, density in enumerate(initDensity):
        vecStartAbundance = np.array([1, density])
        initCells = np.sum(vecStartAbundance)
        numDraws = nbrCellsEnd - initCells
        dist = polya_probability_two_species(numDraws, vecStartAbundance)      
        x = np.linspace(0, 1.0, len(dist))
        
        averageRichness.append(initCells * np.sum(dist[x >= coexistence[0]]))
        
        ax1.plot(x, dist, color=colors[i], label=f"{initCells}")
        
    ax1.set_xlabel('Final fraction of an initial species')
    ax1.set_ylabel('Probability')
    ax1.axvline(x = coexistence[0], color='k', ls=':')
    ax1.text(0.375, 0.0018, 'Presence', color='k')
    ax1.legend(title='Number of species', framealpha=0.9)
    ax1.set_xlim(0.0, 1.0)
    
    ax2.plot(initDensity + 1, averageRichness)

    ax2.set_xlabel('Initial number of the species')
    ax2.set_ylabel('Mean richness')

    filename = 'polya_richness'
    plt.tight_layout()
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
        
    return 0


def distribution_time(vecStartAbundance, fixateFraction):
    # define regimes of fixation of coexistence
    coexistence = (fixateFraction, 1. - fixateFraction)
    
    # density arrays for results
    nbrGenerations = [13, 130, 1300, 13000]
    probDensityLow = []
    probDensityHigh = []
    width = []
    mean = []
    cmapname = 'cividis'
    cmap = matplotlib.cm.get_cmap(cmapname)
    colors = [cmap(nbr) for nbr in np.linspace(0.0, 1.0, num=len(nbrGenerations))]
    
    # define figure
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
    
    for i, nbrGen in enumerate(nbrGenerations):
        initCells = np.sum(vecStartAbundance)
        numDraws = initCells * ( 2 ** numGenerations - 1)
        numDraws = nbrGen
        #dist = polya_two_species_distribution(numDraws, vecStartAbundance)
        dist = polya_probability_two_species(numDraws, vecStartAbundance)     
        x = np.linspace(0, 1.0, len(dist))
        
        probDensityLow.append(np.sum(dist[x <= coexistence[0]]))
        probDensityHigh.append(np.sum(dist[ x >= coexistence[1]]))
        sample = np.arange(0, len(dist))
        sample = x
        #width.append(np.sqrt(np.sum((sample ** 2 * dist )) / np.sum((sample *  dist)) ** 2 - 1) )
        mean.append(np.sum(sample *  dist))
        width.append(np.sqrt(np.sum( sample ** 2 * dist ) - np.sum((sample *  dist) ** 2)))
        
        dist /= np.max(dist)
        ax1.plot(x, dist, color=colors[i], label=f"{nbrGen}")
    ax1.set_xlabel('Final fraction eGFP')
    ax1.set_ylabel('Probability')
    ax1.axvline(x = coexistence[0], color=settings.colors['mCherry'], ls=':')
    ax1.axvline(x = coexistence[1], color=settings.colors['eGFP'], ls=':')
    ax1.text(0.4, 0.00011, 'Coexistence', color=settings.colors['coexistence'])
    ax1.text(0.025, 0.00011, 'mCherry fix.', color=settings.colors['mCherry'])
    ax1.text(0.825, 0.00011, 'eGFP fix.', color=settings.colors['eGFP'])
    ax1.legend(title='Number of draws', loc=(0.01, 0.1), framealpha=0.9)
    # ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.0, 1.0)
    axins = ax1.inset_axes([0.6, 0.6, 0.37, 0.37])
    axins.plot(nbrGenerations, width)
    #axins.plot(nbrGenerations, mean)
    axins.set_xscale('log')
    #axins.set_yscale('log')
    axins.set_xlabel('Number of draws')
    axins.set_ylabel('Variance')
    
    ax2.plot(nbrGenerations, probDensityLow, color=settings.colors['mCherry'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['mCherry'],
             label=settings.labels['mCherry'])
    ax2.plot(nbrGenerations, probDensityHigh, color=settings.colors['eGFP'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['eGFP'],
             label=settings.labels['eGFP'])
    ax2.plot(nbrGenerations,
             1. - np.array(probDensityLow) - np.array(probDensityHigh),
             color=settings.colors['coexistence'],
             ls=settings.linestyles['sim'],
             marker=settings.markers['coexistence'],
             label=settings.labels['coexistence'])
    ax2.set_xlabel('Number of draws')
    ax2.set_ylabel('Probability')
    ax2.legend()
    ax2.set_xscale('log')
    ax2.set_ylim(0.0, 1.0)
    # ax2.set_xlim(0, 21)
    
    filename = 'polya_time_distribution'
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
        
    return 0
    
  
if __name__ == '__main__':
    directory = os.getcwd() + os.sep + FIGDIR
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    # testing_polya(100000, 100, np.array([2, 2]), np.array([1, 1]))
    
    numGenerations = 2500
    fixateFraction = 1. / 5.
    
    # density
    # initial_density_plot(numGenerations, fixateFraction)
    
    # selective advantage
    initAbundance = np.array([1, 1])
    fitnessDiff = np.logspace(0.0, 0.30102999566, 11)
    advantage = [1.1, 1.0]
    selection_advantage_plot(20000, numGenerations, initAbundance, fitnessDiff,
                            fixateFraction)
    
    # density plot with selection
    initial_density_growth_plot(20000, numGenerations, fixateFraction, advantage)
    
    #distribution_time(initAbundance, fixateFraction)
    
    # multispecies_coexistence(numGenerations, fixateFraction)
    