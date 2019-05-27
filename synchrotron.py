import subprocess
import numpy as np
import numpy.ctypeslib as ctl
import ctypes as ct
import matplotlib.pyplot as plt
import statistics
import seaborn as sns
import random
import scipy.special as sp
def bk(x):
    return x**2 * sp.kv(2/3,x*0.99994999875/3)**2/11
n = 5000000
subprocess.run(["make"])
main = ctl.load_library("main.so","./")
main.metropolis_spectrum_sy.restype = ct.POINTER(ct.c_double)
x_p = main.metropolis_spectrum_sy(n)
x = ctl.as_array(x_p, shape = (n,))
M = 1000
y = 10*np.array(range(0,M+1))/M

#plt.hist(x,bins=80,normed=True,weights=x,color='r',alpha=0.5,edgecolor='k')
plt.hist(x,bins=50,normed=True,weights=x,color='r')
plt.plot(y,bk(y),color='g')

plt.xlabel(r'$\omega$, plank units.')
plt.ylabel(r'f($\omega$,$\theta = 0$), a.u.')
plt.xlim(0,10)
#plt.text(-0.023,30,r'$\sigma = %.5f$'%(np.std(x)),fontsize=15,color='k')
#plt.text(0,0.6,r'$\sin(\theta)^{2}$',fontsize=15,color='blue')
#plt.xlabel('omega, mc^2/h')
#plt.ylabel('f(omega), a.u.')
#print(np.std(x))
plt.tick_params(axis='both',direction='in')
plt.show()
