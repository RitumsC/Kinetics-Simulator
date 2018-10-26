"""
Script to generate input files at different urea concentrations

requires premade directory ./protein_input/
"""
import math

c_urea = []
for i in range(21):
    c_urea.append(8/20*i)

for i in c_urea:
    k1 = 26000*math.exp(-1.68*i)
    k2 = 6.0e-2*math.exp(+0.95*i)
    k3 = 730*math.exp(-1.72*i)
    k4 = 7.5e-4*math.exp(+1.20*i)
    f = open("./protein_input/protein"+str(i)+".txt", 'w+')
    print("Wrote: ", "./protein_input/protein:"+str(i)+".txt")
    f.write("D=1,I=0,N=0\n")
    f.write("D|I|"+str(k1)+"\n")
    f.write("I|D|"+str(k2)+"\n")
    f.write("I|N|"+str(k3)+"\n")
    f.write("N|I|"+str(k4)+"\n")
    f.close()
