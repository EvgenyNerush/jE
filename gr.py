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
main.metropolis_spectrum_sy.restype = ct.POINTER(ct.c_double)
x_p = main.metropolis_spectrum_sy(n)
x = ctl.as_array(x_p, shape = (n,))
omega = []
for i in range(int(n/2)):
    if i % 2 ==0:
        omega.append(x[i])
#weights=x,
plt.hist(x,bins=90,normed=True,weights=x,color='r')
plt.xlabel('omega, mc^2/h')
plt.ylabel('f(omega), a.u.')
plt.show()
