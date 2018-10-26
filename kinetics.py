"""Python module for kinetic reaction simulation

Supports text file input of elementar reactions and initial conditions
Example input:
    A=5,B=0,C=2,D=0
    A,B|C|10
    C|D|0.5
where first line specifies starting concentrations and subsequent
lines define elementary reactions, so that reaction
A + B -> C with rate 10 becomes A,B|C|10

Method supported:
    readFile(file): parses the input file
    reactions(Reags, Rxns): parses the output of readFile
    run(Reags, Rxns, steps, nspace=1): runs the simulation
    next_euler(Reags, reacts, rlist): calculates concentrations after timestep
"""


import matplotlib.pyplot as plt
import os
import time
from functools import reduce

# default values for simulation steps
h = 1e-3  # simulation timestep
steps = 10000
nspace = 10  # after how many steps to save concentrations for plot


class Rxn:
    """Class that describes each elementary reaction

    Initialised with:
        rate - float defining reaction rate
        reag - list of names of reagents in reaction
        prod - list of names of products in reaction

    Contains two functions:
        nReag(string r) returns the order of r in reaction
        nProd(string r) : returns number of molecules of r produced in reaction
    """
    def __init__(self, rate, reag, prod):
        self.rate = rate
        self.reag = reag
        self.prod = prod

    def nReag(self, rg):
        return self.reag.count(rg)

    def nProd(self, rg):
        return self.prod.count(rg)


def readFile(file):
    """Reads the file given

    Returns a dictionary of reagents with their concentrations - Reags

    Returns a list of Rxn's (see class)
    """
    f = open(file, 'r')

    line = f.readline().split(",")
    Reags = {}
    for reag in line:  # reads in names and concentrations
        reag = reag.split("=")
        Reags[reag[0]] = float(reag[1])

    Rxns = []
    for line in f:  # reads in reactions and initialises rxns
        line = line.split("|")
        line[0] = line[0].split(",")
        line[1] = line[1].split(",")

        reag = [r for r in line[0]]
        prod = [p for p in line[1]]
        rxn = Rxn(float(line[2]), reag, prod)
        Rxns.append(rxn)

    f.close()

    return Reags, Rxns


def reactions(Reags, Rxns):
    """Given Reags and Rxns forms diff eq forms for each reagent

    Returns dictionary of diff eq's in form:
        {Reagent1: [(rate const wrt reagent, index of reaction),
                     ...]],
         Reagent2: ...,
         ...}
    and list of reactions as indexed in the dictionary:
        [[reag1, reag2], ...]

    intended use is to replace the reags by name with their concentrations,
    then multiply together and with rate corected for each reagent to yield
    the respective change
    """
    reactions_dict = {}
    reactions_list = []
    i = 0
    for r in Reags:
        reactions_dict[r] = []
        for rxn in Rxns:
            a = rxn.nReag(r)
            if a:
                if rxn.reag in reactions_list:
                    reactions_dict[r].append((-a*rxn.rate,
                                              reactions_list.index(rxn.reag)))
                else:
                    reactions_list.append(rxn.reag)
                    reactions_dict[r].append((-a*rxn.rate, i))
                    i = i + 1
            a = rxn.nProd(r)
            if a:
                if rxn.reag in reactions_list:
                    reactions_dict[r].append((a*rxn.rate,
                                              reactions_list.index(rxn.reag)))
                else:
                    reactions_list.append(rxn.reag)
                    reactions_dict[r].append((a*rxn.rate, i))
                    i = i + 1
    return (reactions_dict, reactions_list)


def next_euler(Reags, reacts, rlist):
    """Calls the next time step

    As tested with the Oregonator, the small concentration changes are very
    importent, so small time steps are required, therefore when tested with
    Runge-Kutta 4th order integrator the performance was worse.

    Takes in dictionary of reagents Reags and arguments reacts, rlist
    prepared by reactions method.
    """
    vals = [reduce(
            lambda x, y: x*y, [Reags[i] for i in reags]) for reags in rlist]

    for r in reacts:
        delta = 0
        for rate, indx in reacts[r]:
            delta = delta + rate*vals[indx]
        Reags[r] = Reags[r] + h*delta
    return


def run(Reags, Rxns, steps, nspace=1):
    """Simulates the system for specified steps and returns concentrations
    at each specified time.
    """
    (reacts, rlist) = reactions(Reags, Rxns)
    points, i = [], 0
    for r in Reags:
        points.append([])

    for j in range(steps):
        if j % nspace:
            next_euler(Reags, reacts, rlist)
        else:
            for k in Reags:
                points[i].append(Reags[k])
                i = i + 1
            i = 0
    return points


# used for testing the kinetic model used
def main():
    t1 = time.time()
    Reags, Rxns = readFile("reaction.txt")
    points = run(Reags, Rxns, steps, nspace)

    ax = plt.axes()
    i = 0
    for r in Reags:
        ax.plot([nspace*x*h for x in range(int(steps/nspace))], points[i],
                label=r)
        i = i + 1

    ax.legend()
#    ax.set_yscale('log')
    ax.set_xlabel("time / s")
    ax.set_ylabel("concentration / M")
    plt.show()
    t2 = time.time()
    print(t2-t1)


if __name__ == '__main__':
    main()
