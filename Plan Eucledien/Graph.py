import numpy as np
import scipy.optimize
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

asymptotes = []

for density, v in sorted_data:
    N, mean, err = zip(*v)
    N = np.asarray(N)
    mean = np.asarray(mean)
    err = np.asarray(err)
    match density:
        case 1.4365 | 1.437 | 1.4375 | 1.438:
            asymptotes.append([density, mean[-2], err[-2]])
        case 1.439 | 1.440:
            asymptotes.append([density, mean[-1], err[-1]])
    if density == 1.43635:
        line = plt.errorbar(N, mean, yerr=err, capsize=3, capthick=1)
        line.set_label(str(density))
    N = np.log(N)
    err /= mean
    mean = np.log(mean)
    a = (mean[1] - mean[0]) / (N[1] - N[0])
    err_a = (err[0]**2 + err[1]**2)**0.5 / (N[1] - N[0])
    print(f"{density}: {2 + 2*a:.5f} ± {2*err_a:.5f}")

plt.legend()
plt.xscale("log")
plt.yscale("log")
#plt.gca().set_aspect('equal')
#plt.ylim(0, 1)
plt.show()

density, mean, err = zip(*asymptotes)
line = plt.errorbar(density, mean, yerr=err, capsize=3, capthick=1, elinewidth = 1, linewidth=0, zorder=2.1, fmt='o', markersize=3, color="firebrick")
line.set_label(str(density))
def f(rho, rho_c, p, c):
    return c*abs(rho - rho_c)**p
popt, pcov = scipy.optimize.curve_fit(f,
                                      density,
                                      mean,
                                      sigma=err,
                                      absolute_sigma=True,
                                      bounds=((1.436, 0, 0), (1.4365, 1, 2)))

print(f"Seuil critique: {popt[0]:.5f} ± {pcov[0,0]**0.5:.5f}, exposant: 1/{1/popt[1]:.5f} ± {pcov[1,1]**0.5:.5f}")
rho = np.linspace(popt[0], max(density), 1000)
plt.plot(rho, f(rho, *popt), '--', color="tab:orange")
plt.show()
