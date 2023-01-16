# chemostat-ABM

This is code written to perform an agent-based model of rod-shaped bacteria growing in an open-ended chemostat.
The purpose of this code is to model investigate the competition of populations of bacteria in a confined environment.
The bacteria can be from different species with different physical characteristics such as growth rate, division lengths, length, and radius (amongst others).
There is also the possibility to give some variability to the division lengths and the angles at which the cells divide. 

## Papers

This work was inspired by some very nice simulations developed in these works.

- [Geometry and Mechanics of Microdomains in Growing Bacterial Colonies, *You et al* (2018)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031065)
- [Confinement-induced self-organization in growing bacterial colonies, *You et al* (2021)](https://www.science.org/doi/full/10.1126/sciadv.abc8685)
- [Modeling mechanical interactions in growing populations of rod-shaped bacteria, *Winkle et al* (2017)](https://iopscience.iop.org/article/10.1088/1478-3975/aa7bae/meta)
- [Emergent spatiotemporal population dynamics with cell-length control of synthetic microbial consortia, *Winkle et al* (2021)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009381)

## Getting started

### Dependencies

The code is split in two parts: the ABM simulations are written and run in C++, whereas the analysis is in Python.
The environment I am currently using is

1. **Python**: 3.10.6
  1. numpy 1.21.5
  2. matplotlib 3.6.2
  3. imageio 2.22.1
  4. re 2.2.1
2. **C++**: g++ compiler version 11.3.0, C++17 Standard Library

### Installation

I recommend forking the main repo, and trying to pull often to keep updated with the newest version that is online.
If you make changes to the code that you think are beneficial, please feel free to make a Pull Request!

### Executing program

#### Locally

It's pretty straightforward to run simulation code (I'm assuming that you're running all these commands in console).
Once you've set up the geometry of your simulations with the bacteria you would like to test, first compile your code
```console
foo@bar:~$ g++ -std=c++17 des.cpp -o {executable-name}.o
```
Then run it
```console
foo@bar:~$ ./{executable-name}.o {simulation-name} {sim-number}
```
We have three variables
1. ```{executable-name}``` : string, doesn't matter what you choose so long as you make sure that the you name different types of simulations different names. You wouldn't want to be confused about what simulations are running. For example: ```my-first-des```
2. ```{simulation-name}``` : string, descriptive name of what you want your simulation to be. For example: ```three-species-same-length```
3. ```{executable-name}``` : integer, instance number of the simulation. For all the simulation with the similar ```{simulation-name}```, you should start with zero and then keep climbing up.
Make sure that you do not keep the parenthesis in the names.

#### Niagara

What is important here is to obtain statistics of the simulations.
One instance is not enough to come to conclusions of the competition in the chemostat.
You can run many simualtions on your own machine if you have the hardware for it, however if you have the luxury of access to ComputeCanada resources, then I've provided you with a short script to help run things in parallel.
If you are unfamiliar with niagara, you can look into the documentation [here](https://docs.alliancecan.ca/wiki/Niagara).
I'll be assuming as of now that you are familiar with Niagara enough to run the code below.

First load the correct compiler module in Niagara which (as of right now) is
```console
user@nia-login0X:~$ module load gcc
```
Again, compile your code
```console
user@nia-login0X:~$ g++ -std=c++17 -o {executable-name}.o des.cpp -lstdc++fs
```
Now the script ```niagara_submit.sh``` will help you run as many simulations in parallel as you'd like.
Make sure to request an appropriate amount of resources (nodes, cores, time) for the number of simualtions you're running.
You can do this by changing the values in lines 2-5, remember there are roughly 40 cores to each node on Niagara.
Then, you'll want to change line 22 in
```
parallel --joblog slurm-$SLURM_JOBID.log --sshdelay 0.1 --wd $PWD "./{executable-name}.o {simulation-number} {}" ::: {0..{number-simulations}}
```
Be sure to put in the appropriate names for ```{executable-name}```, ```{simulation-name}```, and ```{number-simulations}``` (without curly brackets).
Leave everything else as is, unless you are very familiar with what is going on.
You can submit a job on niagara by running
```console
user@nia-login0X:~$ sbatch niagara_submit.sh
```
and voila!

## Description

TBD

## Authors

- [Jeremy Rothschild](https://github.com/jbRothschild/)
