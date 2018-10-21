"""To be added

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import time

h = 10e-7


class Rxn:
    """What am I doing??
    """
    def __init__(self, rate, reag, prod):
        self.rate = rate
        self.reag = reag
        self.prod = prod

    def isReag(self, rg):
        return rg in self.reag

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
    reactions = {}
    for r in Reags:
        reactions[r] = []
        for rxn in Rxns:
            if rxn.isReag(r):
                reactions[r].append((-rxn.rate, rxn.reag))
            a = rxn.nProd(r)
            if a:
                reactions[r].append((a*rxn.rate, rxn.reag))
    return reactions


def next_euler(Reags, reactions):  # will be called steps times
    """
    """
    for diff in reactions:
        delta = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] for i in reags])
            delta = delta + rate*np.prod(reags)
        Reags[diff] = Reags[diff] + h*delta
    return Reags


def next_rk4(Reags, reactions):  # will be called steps times
    a, b,  c = {}, {}, {}
    for diff in reactions:
        j = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] for i in reags])
            j = j + rate*np.prod(reags)
        a[diff] = j
    for diff in reactions:
        j = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] + 0.5*h*a[i] for i in reags])
            j = j + rate*np.prod(reags)
        b[diff] = j
    for diff in reactions:
        j = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] + 0.5*h*b[i] for i in reags])
            j = j + rate*np.prod(reags)
        c[diff] = j
    for diff in reactions:
        j = 0
        for i in reactions[diff]:
            rate, reags = i
            reags = np.array([Reags[i] + h*c[i] for i in reags])
            j = j + rate*np.prod(reags)
#            if j < 1: print(j)
        Reags[diff] = Reags[diff] + h*(a[diff]+2*b[diff]+2*c[diff]+j)/6

    return Reags


def euler(f, y, x):
    a = f(y, x)
    return x + h*a


def run(Reags, Rxns, steps, plot=True):
    reacts = reactions(Reags, Rxns)
    points, i = [], 0
    for r in Reags:
        points.append([])

    for j in range(steps):
        for r in Reags:
            points[i].append(Reags[r])
            i = i + 1
        i = 0
        Reags = next_euler(Reags, reacts)
    if plot:
        ax = plt.axes()
        ax.plot([x*h for x in range(steps)], points[4])
        ax.plot([x*h for x in range(steps)], points[5])
        ax.plot([x*h for x in range(steps)], points[6])
#        a = np.array(points[0])
#        b = np.array(points[1])
#        c = np.array(points[2])
#        ax.plot(a+b+c)
        ax.set_yscale('log')
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
    Reags, Rxns = readFile("oregonator.txt")
    print(run(Reags, Rxns, 100000000))
#    urea = []
#    plot2 = []
#    for file in os.listdir("./protein_input/"):
#        urea.append(float(file.split('.txt')[0][7:]))
#        Reags, Rxns = readFile("./protein_input/"+file)
#        plot2.append(run(Reags, Rxns, 10000))
#    plt.plot(urea, plot2)
    t2 = time.time()
    print(t2-t1)


if __name__ == '__main__':
    main()
