"""File that uses kinetic.py to generate the protein output

takes up to 5 minutes to run
"""
import time
import os
import kinetics as kin
import matplotlib.pyplot as plt

# changing default constants
kin.h = 1e-6


t1 = time.time()

urea = []
plot2 = []

# takes starting reagents and concentrations
for file in os.listdir("./protein_input/"):
    rgs = kin.readFile("./protein_input/"+file)[0]

# for each file calculates equilibrium concs given constants
for file in os.listdir("./protein_input/"):
    urea.append(float(file.split('.txt')[0][7:]))

    Rxns = kin.readFile("./protein_input/"+file)[1]
    reacts, rlist = kin.reactions(rgs, Rxns)
    maxCh = 1
    points = dict(rgs)
    while maxCh > 1e-9:
        i = 0
        kin.next_euler(rgs, reacts, rlist)
        if not i % 1000:
            maxCh = max(abs(points[r] - rgs[r]) for r in rgs)
        i = i + 1
        points = dict(rgs)
        # print(maxCh)

    pts = [rgs[r] for r in rgs]
    plot2.append(pts)
    print(file.split('.txt')[0] + " done")

# plots the results
D = [pts[0] for pts in plot2]
I = [pts[1] for pts in plot2]
N = [pts[2] for pts in plot2]

ax = plt.axes()
ax.plot(urea, D, label='D', marker='x')
ax.plot(urea, I, label='I', marker='o')
ax.plot(urea, N, label='N', marker='+')
ax.legend()
ax.set_xlabel("Urea concentration / M")
ax.set_ylabel("molar fraction")
ax.set_title("Protein folding as concentration of urea")
plt.savefig("prot.png")
plt.show()

t2 = time.time()
print("Time taken for calculations: ", t2-t1, " seconds")
