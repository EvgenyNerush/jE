import subprocess
import time
import numpy             as np
import numpy.ctypeslib   as ctl
import ctypes            as ct
import matplotlib.pyplot as plt

from scipy.special import kv

# Излучение с окружности в точку, лежащую в плоскости окружности.  Алгоритмом Метрополиса
# вычисляется спектр излучения (без интегрирования по углам), и сравнивается с формулой (14.83),
# стр. 484, c $\theta = 0$ из [Jackson J.D., Classical electrodynamics, Wiley, 1962].

n = 5000 # number of random numbers to generate
radius = 30;
g = 30;
omega_c = g * g * g / radius;
omega_max = 3 * omega_c
nb = 150 # number of bins, bins from 0 to omega_max
domega = omega_max / nb

subprocess.run(["make", "circletestlib.so"])
lib = ctl.load_library("circletestlib.so", "./")

lib.rnd.restype = ct.POINTER(ct.c_double)
t1 = time.perf_counter()
x_pointer = lib.rnd(n)
t2 = time.perf_counter()
print(str(t2 - t1) + " s")
x = ctl.as_array(x_pointer, shape = (n,))

p = np.zeros(nb)
for i in range(len(x)):
    j = int(x[i] / domega)
    if j >= 0 and j < nb:
        p[j] += 1

omega = domega * (0.5 + np.arange(nb))
xi = omega / (3 * omega_c) # Jackson, Eq. (14.80)
pth = omega * omega * kv(2/3, xi)**2
pth = pth / np.sum(pth)
p = p * omega
p = p / np.sum(p)

plt.figure(figsize = (4,9))

plt.subplot(211)
plt.plot(omega, p, 'r-', label = "metropolis")
plt.plot(omega, pth, 'c-', label = "Jackson (14.83)")
plt.plot([omega_c, omega_c], [0, 0.01], ':', color = 'gray')
plt.xlim(0, omega_max)
plt.xlabel("omega")
plt.ylabel("dN/dx")
plt.legend()

plt.subplot(212)
plt.plot(x, -np.arange(len(x)), '.', ms = 0.7, label = "metropolis")
#plt.xlim(0, omega_max)
plt.xlabel("x[i]")
plt.ylabel("-i")
plt.legend()
plt.savefig("circletest.png")
