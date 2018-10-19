"""To be added

"""
import numpy as np
import matplotlib.pyplot as plt

h = 1e-6


def readFile(file):
    f = open(file, 'r')
    line = f.readline().split(",")
    Reags = {}
    for reag in line:  # reads in concentrations and initialises reags
        reag = reag.split("=")
        Reags[reag[0]] = Reag(float(reag[1]), list())

    Rxns = []
    for line in f:  # reads in reactions and initialises rxns
        line = line.split("|")
        line[0] = line[0].split(",")
        line[1] = line[1].split(",")

        reag, prod = [], []
        for r in line[0]: reag.append(Reags[r])
        for p in line[1]: prod.append(Reags[p])  
        rxn = Rxn(float(line[2]), reag, prod)
        Rxns.append(rxn)
    f.close()

    for react in Rxns:  # adds reactions to reagents
        for r in react.list_r():
            r.add_rxn(react)
    return Reags, Rxns


def next(Reags, Rxns):
    for r in Reags.values():
        delta = 0
        for rxn in r.rxns:
            delta = delta + rxn.change(r)

        r.c = r.c + delta*h
    return Reags, Rxns


class Reag:
    """To be tested
    """
    def __init__(self, conc, rxns):
        self.c = conc
        self.rxns = rxns

    def add_rxn(self, rxn):
        self.rxns.append(rxn)

    def c_n(self, conc):
        self.c = conc


class Rxn:
    """What am I doing??
    """
    def __init__(self, rate, reag, prod):
        self.rate = rate
        self.reag = reag
        self.prod = prod

    def list_r(self):
        return self.reag + self.prod

    def isReag(self, rg):
        return rg in self.reag

    def change(self, rg):
        x = 1
        for r in self.reag: x = r.c*x

        if self.isReag(rg):
            return -1*self.rate*x
        else:
            return self.rate*x


def main():
    file = "./reaction.txt"
    Reags, Rxns = readFile(file)
    points, i = [], 0
    ax = plt.axes()
    for r in Reags.values():
        points.append([])
        
    for j in range(10000):
        for r in Reags.values():
            points[i].append(r.c)
            i = i + 1
        i = 0
        Reags, Rxns = next(Reags, Rxns)
    for r in points: ax.plot(r)
    ax.set_xlabel("time / 1e-3s")
    ax.set_ylabel("concentration / M")
    
    for r in points:
        print(r[-1])
    plt.show()

if __name__ == '__main__':
    main()
