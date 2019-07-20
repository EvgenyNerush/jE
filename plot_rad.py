import subprocess
import numpy as np
import numpy.ctypeslib as ctl
import ctypes as ct
import matplotlib.pyplot as plt
import statistics
import seaborn as sns
import random
import scipy.special as sp
import scipy.stats as st
gamma = 100
beta = np.sqrt(1.0 - 1.0/gamma**2.0)
h = 0.0001
omega0 = h/gamma
def bk(x):
    return x*x * sp.kv(2/3,(x/omega0)*(beta/gamma**3.0)*(1.0/3.0) )**2*0.075# omega^2!!!
n = 50000
subprocess.run(["make"])
main = ctl.load_library("main.so","./")
main.m1.restype = ct.POINTER(ct.c_double)
x_p = main.m1(n)
x = ctl.as_array(x_p, shape = (n,))

main.m2.restype = ct.POINTER(ct.c_double)
y_p = main.m2(n)
y = ctl.as_array(y_p, shape = (n,))


#plt.hist(x,bins=80,normed=True,weights=x,color='r',alpha=0.5,edgecolor='k',label='Numerical results')
xx, gg = np.histogram(x,bins=100)
yy, ggg = np.histogram(y,bins=100)
tmp = 8000*np.array(range(0,len(xx)))/len(xx)
plt.plot(tmp,2*xx,color='b',label=r'$vacuum$',linewidth=1.5)
plt.plot(tmp,yy,color='r',label=r'$vacuum$',linewidth=1.5)

#plt.hist(x,bins=100,density=True,color='g',alpha=0.5,edgecolor='k',label=r'$\epsilon = 100GeV, H = 0.1H_{cr}$')

#plt.plot(y,bk(y),color='b',label='Theory',linewidth=3,linestyle='--')
plt.legend(loc = 'best',prop={'size': 18})
plt.xlabel(r'$\omega$, MeV/$\hbar$.',fontsize = 18)
plt.ylabel(r'f($\omega$,$\theta = 0$), a.u.',fontsize = 18)
#plt.xlim(0,10)
#plt.text(-0.023,30,r'$\sigma = %.5f$'%(np.std(x)),fontsize=15,color='k')
#plt.text(0,0.6,r'$\sin(\theta)^{2}$',fontsize=15,color='blue')
#plt.xlabel('omega, mc^2/h')
#plt.ylabel('f(omega), a.u.')
#print(np.std(x))
#plt.xticks(np.arange(0,9000,2000),(0,r'$1000$',r'$2000$',r'$3000$',r'$4000$'))
#plt.xticks(np.arange(0,18000,2500),(0,r'$1250$',r'$2500$',r'$3750$',r'$5000$',r'$6250$',r'$7500$',r'$8750$'))
#plt.xticks(np.arange(1050,1301,50),(525,550,575,600,625,650))
plt.tick_params(axis='both',direction='in',labelsize = 18)

plt.show()
