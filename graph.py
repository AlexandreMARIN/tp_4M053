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
plt.title("Historiques de convergence pour $A_{"+str(data["Jacobi"]["N"])+"}$ et $b=(1,\ \dots\ , 1)$")
plt.xlabel("$k$")
plt.ylabel("$||r_{k}||_{\infty}/||b||_{\infty}$", fontsize='x-large')
plt.show()


#######time
plt.figure(2)

for n in names:
    plt.plot(data[n]["sizes"], data[n]["duration"])

plt.plot(data[n]["sizes"], np.array(data[n]["sizes"])**3, data[n]["sizes"], np.array(data[n]["sizes"])**4)
plt.legend(names+["$N^3$", "$N^4$"])
plt.ylabel("$\Delta t$ en millisecondes")
plt.xlabel("$N$")
plt.title("durées des méthodes en fonction de la dimension pour $A_N$")
plt.xscale('log')
plt.yscale('log')

plt.show()



##### number of iterations
plt.figure(3)

iterfile = open("iter.json")
niter = json.load(iterfile)
for n in names:
    plt.plot(niter["sizes"], niter[n]["niter"])

plt.plot(niter["sizes"], niter["sizes"], niter["sizes"], np.array(niter["sizes"])**2)
plt.legend(names+["$N$", "$N^2$"])
plt.ylabel("nombre d'itérations")
plt.xlabel("$N$")
plt.title("Nombre d'itérations en fonction de la taille de $A_N$")
plt.xscale('log')
plt.yscale('log')

plt.show()

####LU -- time
plt.figure(4)
LU_file = open("LU.json")
LU_data = json.load(LU_file)
plt.plot(LU_data["sizes"], LU_data["LU_dur"])
plt.plot(LU_data["sizes"], LU_data["Cholesky_dur"])
plt.plot(LU_data["sizes"], np.array(LU_data["sizes"])**3)
plt.legend(["LU", "Cholesky", "$N^3$"])
plt.title("Durées des décompositions pour $A_N$ en fonction de $N$")
plt.xlabel("$N$")
plt.ylabel("durées (en microsecondes)")
plt.xscale('log')
plt.yscale('log')
plt.show()
