import matplotlib.pyplot as plt
import json
import numpy as np

datafile = open("resvec.json")
data = json.load(datafile)

names = ["Jacobi", "GaussSeidel", "Relax"]


#resvec
for n in names:
    plt.plot(data[n]["resvec"])

plt.yscale('log')
plt.legend(["Jacobi", "Gauss-Seidel", "Relaxation"])
plt.title("Historiques de convergence")
plt.show()


#######time
plt.figure(2)

for n in names:
    plt.plot(data[n]["sizes"], data[n]["duration"])

plt.legend(names)
plt.ylabel("$\Delta t$ en millisecondes")
plt.xlabel("$N$")
plt.title("durées des méthodes en fonction de la dimension pour $A_N$")


plt.show()



##### number of iterations
plt.figure(3)

iterfile = open("iter.json")
niter = json.load(iterfile)
for n in names:
    plt.plot(niter["sizes"], niter[n]["niter"])

plt.legend(names)
plt.ylabel("nombre d'itérations")
plt.xlabel("$N$")
plt.title("Nombre d'itérations en fonction de la taille de $A_N$")
plt.xscale('log')
plt.yscale('log')

plt.show()
