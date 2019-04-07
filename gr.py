import subprocess
import numpy as np
import numpy.ctypeslib as ctl
import ctypes as ct
import matplotlib.pyplot as plt
import statistics
import random
def radp(x): return np.sin(x)**2
n = 100000
subprocess.run(["make"])
main = ctl.load_library("main.so","./")
#main.metropolis_radpattern.restype = ct.POINTER(ct.c_double)
main.metropolis_spectrum.restype = ct.POINTER(ct.c_double)
#x_p = main.metropolis_radpattern(n)
x_p = main.metropolis_spectrum(n)
x = ctl.as_array(x_p, shape = (n,))
#y = ctl.as_array(y_p, shape = (n,))

plt.hist(x,bins=180,normed=True,color='r')
#y = np.pi * np.array(range(n)) / n - np.pi /2.0
#plt.plot(y,2*radp(y)/np.pi,color='b') for dipole
#print(np.std(x))
plt.show()
