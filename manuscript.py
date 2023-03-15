from python_src import plotting, analysis
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

def growthComparison(filename, numberTrajectories):
    listDirectories = [30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60]
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    listDirectories.reverse()
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(exp_dir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        
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

    growthRatio = [listDirectories[0] / float(x) for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.plot(growthRatio, gfpFixateFraction, color='forestgreen', marker='o',
             label='GFP fixates')
    ax1.plot(growthRatio, mchFixateFraction, color='firebrick', marker='o',
             label='mCherry fixates')
    ax1.plot(growthRatio, coexistFraction, color='gold', marker='o',
             label='coexistence')
    ax1.set_title(r"growth rate comparison")
    ax1.set_xlabel(r'growth rate GFP / growth rate mCherry')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.9, 2.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
    # plt.show()
    return 0

def geometryComparison(filename, numberTrajectories):
    listDirectories = [356, 376, 396, 416, 436, 456, 476, 496, 516, 536, 556]
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(exp_dir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        
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

    lengthRatio = [float(x) / listDirectories[5] for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.plot(lengthRatio, gfpFixateFraction, color='forestgreen', marker='o',
             label='GFP fixates')
    ax1.plot(lengthRatio, mchFixateFraction, color='firebrick', marker='o',
             label='mCherry fixates')
    ax1.plot(lengthRatio, coexistFraction, color='gold', marker='o',
             label='coexistence')
    ax1.set_title(r"length comparison")
    ax1.set_xlabel(r'length GFP / length mCherry')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(np.min(lengthRatio) - 0.1, np.max(lengthRatio) + 0.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
    # plt.show()
    return 0

def densityComparison(filename, numberTrajectories):
    listDirectories = [1, 5, 9, 13, 17]
    listDirectories = [101, 105, 109, 113, 117, 121]
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(exp_dir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        
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

    density = [x - 100 for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.plot(density, gfpFixateFraction, color='forestgreen', marker='o',
             label='GFP fixates')
    ax1.plot(density, mchFixateFraction, color='firebrick', marker='o',
             label='mCherry fixates')
    ax1.plot(density, coexistFraction, color='gold', marker='o',
             label='coexistence')
    ax1.set_title(r"initial density comparison")
    ax1.set_xlabel(r'initial number of each species')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0, 21)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
    # plt.show()
    return 0


def frictionComparison(filename, numberTrajectories):
    listDirectories = [201, 202, 205, 210, 220]
    gfpFixateFraction = []
    mchFixateFraction = []
    coexistFraction = []
    
    for exp_nbr in listDirectories:
        numberCoexist = 0.
        numberGfp = 0.
        numberMch = 0.
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(exp_dir, numberTrajectories)
        lastTimepoint = data[:, :, -1]
        
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

    lengthRatio = [ float(x) - 200 / (listDirectories[0] - 200) 
                   for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.plot(lengthRatio, gfpFixateFraction, color='forestgreen', marker='o',
             label='GFP fixates')
    ax1.plot(lengthRatio, mchFixateFraction, color='firebrick', marker='o',
             label='mCherry fixates')
    ax1.plot(lengthRatio, coexistFraction, color='gold', marker='o',
             label='coexistence')
    ax1.set_title(r"A22 increased damping")
    ax1.set_xlabel(r'mCherry damping / eGFP damping')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(np.min(lengthRatio) - 0.1, np.max(lengthRatio) + 0.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
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
        exp_dir = os.getcwd() + os.sep + 'data' + os.sep + f'c_exp_{exp_nbr}'
        data = analysis.collect_data_array(exp_dir, numberTrajectories)
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
    
    ax1.legend(title='# initial species')
    fig1.tight_layout()
    fig1.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename
                 + "_dist.pdf")
    fig1.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename
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
    fig2.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename
                 + "_violin.pdf")
    fig2.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename
                 + "_violin.png")    
    
    # SAD?
    return 0

    
if __name__ == '__main__':
    # growthComparison("compare_growth_rate", 5000)
    # geometryComparison("compare_geometry", 5000)
    # densityComparison("compare_density", 5000)
    # frictionComparison("compare_friction", 5000)
    distributionInitSpecies("lanes_formed", 5000)
    