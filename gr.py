import subprocess
import numpy as np
import numpy.ctypeslib as ctl
import ctypes as ct
import matplotlib.pyplot as plt
import statistics
import random
def radp(x): return np.sin(x)**2
n = 1000000
subprocess.run(["make"])
main = ctl.load_library("main.so","./")
main.metropolis_spectrum.restype = ct.POINTER(ct.c_double)
x_p = main.metropolis_spectrum(n)
x = ctl.as_array(x_p, shape = (n,))
plt.hist(x,bins=60,normed=True,color='r')

plt.show()
