import matplotlib.pyplot as plt
import json
import numpy as np

datafile = open("resvec.json")
data = json.load(datafile)


plt.plot(data["Jacobi"]["resvec"])
plt.plot(data["GaussSeidel"]["resvec"])
plt.plot(data["Relax"]["resvec"])
plt.yscale('log')
plt.legend(["Jacobi", "Gauss-Seidel", "Relaxation"])
plt.show();
