"""File that uses kinetic.py to generate the oregonator output

takes up to 10 minutes to run
"""
import time
import kinetics as kin
import matplotlib.pyplot as plt

nspace = 1000
steps = 45000000
kin.h = 20e-7


t1 = time.time()
Reags, Rxns = kin.readFile("oregonator.txt")
points = kin.run(Reags, Rxns, steps, nspace)
ax = plt.axes()
ax.plot([nspace*x*kin.h for x in range(int(steps/nspace))], points[4],
        label='X')
ax.plot([nspace*x*kin.h for x in range(int(steps/nspace))], points[5],
        label='Y')
ax.plot([nspace*x*kin.h for x in range(int(steps/nspace))], points[6],
        label='Z')
ax.legend()
ax.set_title('Oregonator relevant species concentrations with time')
ax.set_yscale('log')
ax.set_xlabel("time / s")
ax.set_ylabel("concentration / M")
plt.savefig("test.png")
ax.plot()
t2 = time.time()
print(t2-t1)
