from python_src import plotting, analysis
import os
import matplotlib.pyplot as plt
import  numpy as np

def growthComparison(filename, numberTrajectories):
    listDirectories = [30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60]
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

    growthRatio = [float(x) / listDirectories[-1] for x in listDirectories]
    fig, ax1 = plt.subplots(1)
    ax1.plot(growthRatio, gfpFixateFraction, color='forestgreen', marker='o',
             label='GFP fixates')
    ax1.plot(growthRatio, mchFixateFraction, color='firebrick', marker='o',
             label='mCherry fixates')
    ax1.plot(growthRatio, coexistFraction, color='gold', marker='o',
             label='coexistence')
    ax1.set_title(r"fixation and coexistence")
    ax1.set_xlabel(r'growth rate mCherry / growth rate GFP')
    ax1.set_ylabel(r'fraction')
    ax1.set_ylim(0.0, 1.0)
    ax1.set_xlim(0.4, 1.1)
    
    plt.legend()
    fig.tight_layout()
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".pdf")
    plt.savefig(os.getcwd() + os.sep + 'data' + os.sep + filename + ".png")
    # plt.show()
    return 0
    
if __name__ == '__main__':
    growthComparison("compare_growth_rate", 5000)