"""To be added

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from functools import reduce

h = 1e-6
steps = 1000000
defFile = "oregonator.txt"
nspace = 10


class Rxn:
    """What am I doing??
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


def next_euler(Reags, reacts, rlist):  # will be called steps times
    vals = [reduce(
            lambda x, y: x*y, [Reags[i] for i in reags]) for reags in rlist]

    for r in reacts:
        delta = 0
        for rate, indx in reacts[r]:
            delta = delta + rate*vals[indx]
        Reags[r] = Reags[r] + h*delta
    return 0


def run(Reags, Rxns, steps, nspace=1):
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


def main():
    t1 = time.time()
    Reags, Rxns = readFile("oregonator.txt")
    points = run(Reags, Rxns, steps, nspace)
    ax = plt.axes()
    ax.plot([nspace*x*h for x in range(int(steps/nspace))], points[4],
            label='X')
    ax.plot([nspace*x*h for x in range(int(steps/nspace))], points[5],
            label='Y')
    ax.plot([nspace*x*h for x in range(int(steps/nspace))], points[6],
            label='Z')
    ax.legend()
    ax.set_yscale('log')
    ax.set_xlabel("time / s")
    ax.set_ylabel("concentration / M")
    plt.savefig("ayy.png")
    ax.plot()
    t2 = time.time()
    print(t2-t1)


if __name__ == '__main__':
    main()
