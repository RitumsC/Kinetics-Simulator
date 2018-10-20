"""To be added

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import time

h = 10e-6


class Rxn:
    """What am I doing??
    """
    def __init__(self, rate, reag, prod):
        self.rate = rate
        self.reag = reag
        self.prod = prod


    def isReag(self, rg):
        return rg in self.reag

    def isProd(self, rg):
        return rg in self.prod



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


def next(Reags, Rxns):  # will be called steps times
    reactions = {}
    for r in Reags:
        reactions[r] = []
        for rxn in Rxns:
            if rxn.isReag(r):
                reactions[r].append((-rxn.rate, rxn.reag))
            elif rxn.isProd(r):
                reactions[r].append((rxn.rate, rxn.reag))
    # now we have a list  of reactions to work with
    for diff in reactions:
        delta = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] for i in reags])
            delta = delta + rate*np.prod(reags)
        Reags[diff] = Reags[diff] + h*delta
    # print([Reags[i] for i in Reags])
    return Reags, Rxns


def rk4(f, y, x):
    a = f(y, x)
    b = f(y, x + h/2*a)
    c = f(y, x + h/2*b)
    d = f(y, x + h*c)
    return x + h/6*(a+2*b+2*c+d)


def euler(f, y, x):
    a = f(y, x+h)
    return x + h*a


def run(Reags, Rxns, steps):
    points, i = [], 0
    for r in Reags:
        points.append([])

    for j in range(steps):
        for r in Reags:
            points[i].append(Reags[r])
            i = i + 1
        i = 0
        Reags, Rxns = next(Reags, Rxns)

    ax = plt.axes()
    for r in Reags:
        ax.plot(points[i])
        i = i+1
    # ax.set_yscale('log')
    ax.set_xlabel("time / s")
    ax.set_ylabel("concentration / M")
    i = 0
    ret = [0, 0, 0]
    for r in points:
        ret[i % 3] = r[-1]
        i = i + 1
    plt.show()
    return ret



def main():
    t1 = time.time()
#    Reags, Rxns = readFile("reaction.txt")
#    print(run(Reags, Rxns, 100000))
    plot2 = []
    for file in sorted(os.listdir("./protein_input/")):
        Reags, Rxns = readFile("./protein_input/"+file)
        plot2.append(run(Reags, Rxns, 100000))
    plt.plot(plot2)
    t2 = time.time()
    print(t2-t1)

if __name__ == '__main__':
    main()


















