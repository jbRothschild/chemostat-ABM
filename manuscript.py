from python_src import plotting, analysis, settings
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
import scipy.stats as stats

FIGDIR = 'figure'

plt.style.use('python_src/custom.mplstyle')

def timeDistributionPlot(filename, numberTrajectories):
    dataDir = 456
    expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{dataDir}'
    params = pd.read_csv(expDir + os.sep + 'params.txt')
    timestep = params.iloc[0][" save_time"] / 60.
    data = analysis.collect_data_array(expDir, numberTrajectories)
    # max_t = np.shape(data)[2] * timestep
    max_t = 16.
    nbr_species = np.shape(data)[1]
    extinctions = [[] for _ in range(nbr_species)]
    nbr_extinctions = np.zeros((nbr_species))
    for i in np.arange(0, numberTrajectories):
        for j in np.arange(0, nbr_species):
            zeros = np.where(data[i, j, :] == 0.0)[0]
            if zeros.size > 0:
                extinctions[j].append(zeros[0] * timestep)
                nbr_extinctions[j] += 1

    # figure for distribution of each species extinction times
    fig, ax = plt.subplots(1)
    """
    for j in np.arange(0, nbr_species):
        ax.hist(extinctions[j], num_bins,
                facecolor=settings.colors['coexistence'],
                alpha=0.5,
                density=True)
        ax.axvline(np.mean(extinctions[j]), color=BACT_COL[str(j)],
                   linestyle='dashed', linewidth=1)
    """
    
    num_bins = 20
    n, x, _ = ax.hist(extinctions, num_bins,
                      facecolor=settings.colors['single_fixate'],
                      # alpha=0.5,
                      density=True)

    density = stats.gaussian_kde(extinctions[0] + extinctions[1])
    plt.plot(x, density(x), c=settings.colors['single_fixate'],
             ls=settings.linestyles['sim'])
    ax.set_title(r"simulation")
    
    left, bottom, width, height = [0.15, 0.5, 0.4, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])
    totalExtinctions = np.sum(nbr_extinctions)
    ax2.pie(np.array([totalExtinctions, numberTrajectories - totalExtinctions]),
            labels=['fixation', 'coexistence'],
            colors=[settings.colors['single_fixate'],
                    settings.colors['coexistence']],
            #labeldistance=0.2,
            wedgeprops = {"linewidth": 1, "edgecolor": "white"}
            )
    ax2.axis('off')
    

    ax.set_xlim([0.0, max_t])
    #plt.yscale('log')
    ax.set_ylabel(r'probability')
    ax.set_xlabel(r'fixation time, $h$')
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
    # plt.show()
    return 0

def growthComparison(filename, numberTrajectories, bootstrap = False):
    listDirectories = [30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60]
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    listDirectories.reverse()
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(expDir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        
        """
        for timepoint in lastTimepoint:
            if timepoint[0] > 0.:
                if timepoint[1] > 0.:
                    numberCoexist += 1.
                else:
                    numberGfp += 1.
            elif timepoint[1] > 0.:
                numberMch += 1.
        gfpFixateFraction.append(numberGfp / numberTrajectories)
        mchFixateFraction.append(numberMch / numberTrajectories)
        coexistFraction.append(numberCoexist / numberTrajectories)             
        """
        randomIdx = np.random.permutation(len(lastTimepoint))
        bootCoexist = []
        bootGfp = []
        bootMch = []
        for splitIdx in np.split(randomIdx, 100):
            numGfp = 0
            numMch = 0
            numCoexist = 0
            for timepoint in lastTimepoint[splitIdx]:
                if timepoint[0] > 0.:
                    if timepoint[1] > 0.:
                        numCoexist += 1.
                    else:
                        numGfp += 1.
                elif timepoint[1] > 0.:
                    numMch += 1.
            bootGfp.append(numGfp / len(splitIdx))
            bootMch.append(numMch / len(splitIdx))
            bootCoexist.append(numCoexist / len(splitIdx))
        gfpFixateFraction.append(np.mean(bootGfp))
        mchFixateFraction.append(np.mean(bootMch))
        coexistFraction.append(np.mean(bootCoexist))
        if bootstrap:
            yerrGfp = np.std(bootGfp)
            yerrMch = np.std(bootGfp)
            yerrCoexist = np.std(bootCoexist)
        else:
            yerrGfp = 0
            yerrMch = 0
            yerrCoexist = 0

    growthRatio = [listDirectories[0] / float(x) for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.errorbar(growthRatio, coexistFraction, yerr=yerrCoexist,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['coexistence'],
                 marker=settings.markers['coexistence'],
                 label=settings.labels['coexistence'])
    ax1.errorbar(growthRatio, gfpFixateFraction, yerr=yerrGfp,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['eGFP'],
                 marker=settings.markers['eGFP'],
                 label=settings.labels['eGFP'])
    ax1.errorbar(growthRatio, mchFixateFraction, yerr=yerrMch,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['mCherry'],
                 marker=settings.markers['mCherry'],
                 label=settings.labels['mCherry'])
    ax1.set_title(r"growth rate comparison")
    ax1.set_xlabel(r'growth rate ratio (eGFP/mCherry)')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.9, 2.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + ".png")
    # plt.show()
    return 0

def geometryComparison(filename, numberTrajectories, bootstrap=False):
    listDirectories = [356, 376, 396, 416, 436, 456, 476, 496, 516, 536, 556]
    diff = 0
    add_title = r''
    # listDirectories = [656, 676, 696, 716, 736, 756, 776, 796, 816, 836, 856]
    # diff = 300
    # add_title = r'_growthdiff'
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(expDir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        """
        for timepoint in lastTimepoint:
            if timepoint[0] > 0.:
                if timepoint[1] > 0.:
                    numberCoexist += 1.
                else:
                    numberGfp += 1.
            elif timepoint[1] > 0.:
                numberMch += 1.
            
        gfpFixateFraction.append(numberGfp / numberTrajectories)
        mchFixateFraction.append(numberMch / numberTrajectories)
        coexistFraction.append(numberCoexist / numberTrajectories)
        """
        randomIdx = np.random.permutation(len(lastTimepoint))
        bootCoexist = []
        bootGfp = []
        bootMch = []
        for splitIdx in np.split(randomIdx, 100):
            numGfp = 0
            numMch = 0
            numCoexist = 0
            for timepoint in lastTimepoint[splitIdx]:
                if timepoint[0] > 0.:
                    if timepoint[1] > 0.:
                        numCoexist += 1.
                    else:
                        numGfp += 1.
                elif timepoint[1] > 0.:
                    numMch += 1.
            bootGfp.append(numGfp / len(splitIdx))
            bootMch.append(numMch / len(splitIdx))
            bootCoexist.append(numCoexist / len(splitIdx))
        gfpFixateFraction.append(np.mean(bootGfp))
        mchFixateFraction.append(np.mean(bootMch))
        coexistFraction.append(np.mean(bootCoexist))
        if bootstrap:
            yerrGfp = np.std(bootGfp)
            yerrMch = np.std(bootGfp)
            yerrCoexist = np.std(bootCoexist)
        else:
            yerrGfp = 0
            yerrMch = 0
            yerrCoexist = 0

    lengthRatio = [float(x - diff) / (listDirectories[5] - diff)
                   for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.errorbar(lengthRatio, coexistFraction, yerr=yerrCoexist,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['coexistence'],
                 marker=settings.markers['coexistence'],
                 label=settings.labels['coexistence'])
    ax1.errorbar(lengthRatio, gfpFixateFraction, yerr=yerrGfp,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['eGFP'],
                 marker=settings.markers['eGFP'],
                 label=settings.labels['eGFP'])
    ax1.errorbar(lengthRatio, mchFixateFraction, yerr=yerrMch,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['mCherry'],
                 marker=settings.markers['mCherry'],
                 label=settings.labels['mCherry'])
    ax1.set_title(r'simulation')
    ax1.set_xlabel(r'length ratio (eGFP/mCherry)')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(np.min(lengthRatio) - 0.1, np.max(lengthRatio) + 0.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                +".png")
    # plt.show()
    return 0


def densityComparison(filename, numberTrajectories, bootstrap=False):
    listDirectories = [1, 5, 9, 13, 17]
    listDirectories = [101, 105, 109, 113, 117, 121]
    diff = 100
    add_title = r'_growthdiff'
    
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(expDir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        """
        for timepoint in lastTimepoint:
            if timepoint[0] > 0.:
                if timepoint[1] > 0.:
                    numberCoexist += 1.
                else:
                    numberGfp += 1.
            elif timepoint[1] > 0.:
                numberMch += 1.         
        gfpFixateFraction.append(numberGfp / numberTrajectories)
        mchFixateFraction.append(numberMch / numberTrajectories)
        coexistFraction.append(numberCoexist / numberTrajectories)
        """
        randomIdx = np.random.permutation(len(lastTimepoint))
        bootCoexist = []
        bootGfp = []
        bootMch = []
        for splitIdx in np.split(randomIdx, 100):
            numGfp = 0
            numMch = 0
            numCoexist = 0
            for timepoint in lastTimepoint[splitIdx]:
                if timepoint[0] > 0.:
                    if timepoint[1] > 0.:
                        numCoexist += 1.
                    else:
                        numGfp += 1.
                elif timepoint[1] > 0.:
                    numMch += 1.
            bootGfp.append(numGfp / len(splitIdx))
            bootMch.append(numMch / len(splitIdx))
            bootCoexist.append(numCoexist / len(splitIdx))
        gfpFixateFraction.append(np.mean(bootGfp))
        mchFixateFraction.append(np.mean(bootMch))
        coexistFraction.append(np.mean(bootCoexist))
        if bootstrap:
            yerrGfp = np.std(bootGfp)
            yerrMch = np.std(bootGfp)
            yerrCoexist = np.std(bootCoexist)
        else:
            yerrGfp = 0
            yerrMch = 0
            yerrCoexist = 0
        

    density = [x - diff for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.errorbar(density, coexistFraction, yerr=yerrCoexist,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['coexistence'],
                 marker=settings.markers['coexistence'],
                 label=settings.labels['coexistence'])
    ax1.errorbar(density, gfpFixateFraction, yerr=yerrGfp,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['eGFP'],
                 marker=settings.markers['eGFP'],
                 label=settings.labels['eGFP'])
    ax1.errorbar(density, mchFixateFraction, yerr=yerrMch,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['mCherry'],
                 marker=settings.markers['mCherry'],
                 label=settings.labels['mCherry'])
    ax1.set_title(r"simulation")
    ax1.set_xlabel(r'initial abundance of each species')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0, 21)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                + ".png")
    # plt.show()
    return 0


def frictionComparison(filename, numberTrajectories, bootstrap=False):
    listDirectories = [201, 202, 205, 210, 220]
    diff = 200
    add_title = r''
    listDirectories = [401, 402, 405, 410, 420]
    diff = 400
    add_title = r'_growthdiff'
    #listDirectories = [901, 902, 905, 910, 920]
    #diff = 900
    #add_title = r'_growthdiff_coccus'
    
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(expDir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        """
        for timepoint in lastTimepoint:
            if timepoint[0] > 0.:
                if timepoint[1] > 0.:
                    numberCoexist += 1.
                else:
                    numberGfp += 1.
            elif timepoint[1] > 0.:
                numberMch += 1.  
        gfpFixateFraction.append(numberGfp / numberTrajectories)
        mchFixateFraction.append(numberMch / numberTrajectories)
        coexistFraction.append(numberCoexist / numberTrajectories)
        """
        randomIdx = np.random.permutation(len(lastTimepoint))
        bootCoexist = []
        bootGfp = []
        bootMch = []
        for splitIdx in np.split(randomIdx, 100):
            numGfp = 0
            numMch = 0
            numCoexist = 0
            for timepoint in lastTimepoint[splitIdx]:
                if timepoint[0] > 0.:
                    if timepoint[1] > 0.:
                        numCoexist += 1.
                    else:
                        numGfp += 1.
                elif timepoint[1] > 0.:
                    numMch += 1.
            bootGfp.append(numGfp / len(splitIdx))
            bootMch.append(numMch / len(splitIdx))
            bootCoexist.append(numCoexist / len(splitIdx))
        gfpFixateFraction.append(np.mean(bootGfp))
        mchFixateFraction.append(np.mean(bootMch))
        coexistFraction.append(np.mean(bootCoexist))
        if bootstrap:
            yerrGfp = np.std(bootGfp)
            yerrMch = np.std(bootGfp)
            yerrCoexist = np.std(bootCoexist)
        else:
            yerrGfp = 0
            yerrMch = 0
            yerrCoexist = 0

    lengthRatio = [ float(x) - diff / (listDirectories[0] - diff) 
                   for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.errorbar(lengthRatio, coexistFraction, yerr=yerrCoexist,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['coexistence'],
                 marker=settings.markers['coexistence'],
                 label=settings.labels['coexistence'])
    ax1.errorbar(lengthRatio, gfpFixateFraction, yerr=yerrGfp,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['eGFP'],
                 marker=settings.markers['eGFP'],
                 label=settings.labels['eGFP'])
    ax1.errorbar(lengthRatio, mchFixateFraction, yerr=yerrMch,
                 ls=settings.linestyles['sim'],
                 color=settings.colors['mCherry'],
                 marker=settings.markers['mCherry'],
                 label=settings.labels['mCherry'])
    ax1.set_title(r'simulation')
    ax1.set_xlabel(r'damping ratio (mCherry/eGFP)')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(np.min(lengthRatio) - 0.1, np.max(lengthRatio) + 0.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                + ".pdf")
    plt.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename + add_title
                + ".png")
    # plt.show()
    return 0

def distributionInitSpecies(filename, numberTrajectories):
    listDirectories = [72, 73, 74, 75, 76, 77, 78, 79]
    violinPlot = []

    # scatter plot
    cmap_name = 'viridis'
    cmap = matplotlib.cm.get_cmap(cmap_name)
    colors = [cmap(nbr) for nbr in
              np.linspace(0.0, 0.8, num=len(listDirectories))]
    custom_cycler = cycler(color=colors)
    fig1, ax1 = plt.subplots(1)
    ax1.set_prop_cycle(custom_cycler)
    for nbrSpecies, exp_nbr in enumerate(listDirectories):
        expDir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(expDir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        countRichness = [np.count_nonzero(timep) for timep in lastTimepoint]
        distRichness = np.bincount(countRichness) / numberTrajectories
        richness = np.arange(len(distRichness))
        violinPlot.append(countRichness)
        ax1.plot(richness, distRichness, label= 2 + nbrSpecies)
    
    ax1.set_title(r"diversity various initial seeding species")
    ax1.set_xlabel(r'number of species')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 0.8)
    ax1.set_xlim(0, 9)
    
    ax1.legend(title='\# initial species')
    fig1.tight_layout()
    fig1.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename
                 + "_dist.pdf")
    fig1.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename
                 + "_dist.png")
    
    # Violin plot
    fig2, ax2 = plt.subplots(1)
    plots = ax2.violinplot(violinPlot,[x - 70 for x in listDirectories],
                           showmeans=True, showextrema=False, points=1000)
    
    # Set the color of the violin patches
    for pc, color in zip(plots['bodies'], colors):
        pc.set_facecolor(color)
    
    ax2.set_title(r"diversity various initial seeding species")
    ax2.set_xlabel(r'number of initial species')
    ax2.set_ylabel(r'number final species')
    
    fig2.tight_layout()
    fig2.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename
                 + "_violin.pdf")
    fig2.savefig(os.getcwd() + os.sep + FIGDIR + os.sep + filename
                 + "_violin.png")    
    
    # SAD?
    return 0

    
if __name__ == '__main__':
    directory = os.getcwd() + os.sep + FIGDIR
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    timeDistributionPlot("time_distribution", 5000)
    growthComparison("compare_growth_rate", 5000, True)
    geometryComparison("compare_geometry", 5000, True)
    densityComparison("compare_density", 5000, True)
    frictionComparison("compare_friction", 5000, True)
    # distributionInitSpecies("lanes_formed", 5000)
    