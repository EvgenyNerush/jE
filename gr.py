import subprocess
import numpy as np
import numpy.ctypeslib as ctl
import ctypes as ct
import matplotlib.pyplot as plt
import statistics
import seaborn as sns
import random
import scipy.special as sp
gamma = 100
beta = np.sqrt(1.0 - 1.0/gamma**2.0)
h = 0.0001
omega0 = h/gamma
def bk(x):
    return x*x * sp.kv(2/3,(x/omega0)*(beta/gamma**3.0)*(1.0/3.0) )**2*0.075# omega^2!!!
n = 100000
subprocess.run(["make"])
main = ctl.load_library("main.so","./")
main.metropolis_spectrum_sy.restype = ct.POINTER(ct.c_double)
x_p = main.metropolis_spectrum_sy(n)
x = ctl.as_array(x_p, shape = (n,))
M = 10000
y = 13*np.array(range(0,M+1))/M

#plt.hist(x,bins=80,normed=True,weights=x,color='r',alpha=0.5,edgecolor='k')
plt.hist(x,bins=100,density=True,color='b',alpha=0.5,edgecolor='k',label='Metropolis')
#plt.plot(y,bk(y),color='b',label=r'$\frac{dI}{d\Omega}(\omega,\theta = 0)$',linewidth=3,linestyle='--')
plt.legend(loc = 'best',prop={'size': 18})
plt.xlabel(r'$\omega$, MeV/$\hbar$.',fontsize = 18)
plt.ylabel(r'f($\omega$,$\theta = 0$), a.u.',fontsize = 18)
#plt.xlim(0,10)
#plt.text(-0.023,30,r'$\sigma = %.5f$'%(np.std(x)),fontsize=15,color='k')
#plt.text(0,0.6,r'$\sin(\theta)^{2}$',fontsize=15,color='blue')
#plt.xlabel('omega, mc^2/h')
#plt.ylabel('f(omega), a.u.')
#print(np.std(x))
#plt.xticks(np.arange(0,20001,5000),(0,r'$2500$',r'$5000$',r'$7500$',r'$10000$'))
#plt.xticks(np.arange(1050,1301,50),(525,550,575,600,625,650))
plt.tick_params(axis='both',direction='in',labelsize = 18)
plt.show()
