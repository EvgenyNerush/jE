import subprocess
import time
import numpy             as np
import numpy.ctypeslib   as ctl
import ctypes            as ct
import matplotlib.pyplot as plt

n = 100000 # number of particles
nb = 300 # number of bins, bins from -1 to 1

subprocess.run(["make"])
lib = ctl.load_library("lib.so", "./")

lib.metropolis.restype = ct.POINTER(ct.c_double)
t1 = time.perf_counter()
x_pointer = lib.metropolis(n)
t2 = time.perf_counter()
print(str(t2 - t1) + " s")
x = ctl.as_array(x_pointer, shape = (n,))

p = np.zeros(nb)
for i in range(len(x)):
    j = int(round((x[i] + 1) * nb / 2))
    if j >= 0 and j < nb:
        p[j] += 1
p = p / np.sum(p)

xth = np.arange(nb) * 2 / nb - 1
#pth = np.exp(-100 * (xth - 0.5) * (xth - 0.5)) + (1 - xth)**10 * (1 + np.tanh(50 * xth));
pth = np.exp(-500 * (xth + 0.5) * (xth + 0.5)) + np.exp(-500 * (xth - 0.5) * (xth - 0.5));
pth = pth / np.sum(pth)

plt.figure(figsize = (4,9))

plt.subplot(211)
plt.plot(xth, p, 'b-', label = "metropolis")
plt.plot(xth, pth, 'y--', label = "target function")
plt.xlim(-1, 1)
plt.xlabel("x")
plt.ylabel("dN/dx")
plt.legend()

plt.subplot(212)
plt.plot(x, -np.arange(len(x)), '.', ms = 0.7, label = "metropolis")
plt.xlim(-1, 1)
plt.xlabel("x[i]")
plt.ylabel("-i")
plt.legend()
plt.savefig("metropolis_example.png")
