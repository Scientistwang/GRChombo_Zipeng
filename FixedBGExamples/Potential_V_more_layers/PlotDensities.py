import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("ProcaDensities.dat")
t_list = [row[0] for row in data]
rho_list= [row[1] for row in data]
plt.plot(t_list,rho_list)

plt.yscale("log")
plt.xlabel("time")
plt.ylabel("rho")
filename = "rho_vs_t.png"
plt.savefig(filename)
print("...saved the plot at", filename)

