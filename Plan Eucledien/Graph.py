import numpy as np
import matplotlib.pyplot as plt
import os

data = {}
for d in os.walk("Data"):
    for f in d[2]:
        file = open(d[0] + "/" + f)
        for line in file.readlines():
            N, density, fraction = line.split()
            N = int(N)
            density = float(density)
            fraction = float(fraction)
            if density not in data.keys():
                data[density] = {}
            
            if N not in data[density].keys():
                data[density][N] = []
                
            data[density][N].append(fraction)
    
for density, Ns in data.items():
    sub_data = []
    for N, fractions in Ns.items():
        mean = np.mean(fractions)
        if len(fractions) == 1: err = 0
        else: err = np.std(fractions) / np.sqrt(len(fractions))
        sub_data.append((N, mean, err))
    sub_data.sort(key=lambda e: e[0])
    data[density] = sub_data

sorted_data = []
for density, v in data.items():
    sorted_data.append((density, v))
sorted_data.sort(key=lambda e: e[0])

for density, v in sorted_data:
    N, mean, err = zip(*v)
    line = plt.errorbar(N, mean, yerr=err, capsize=3, capthick=1)
    line.set_label(str(density))

plt.legend()
plt.xscale("log")
plt.ylim(0, 1)
plt.show()
