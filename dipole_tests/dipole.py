import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import matplotlib.colors as colors
import matplotlib.cm as mcm
import copy
import scipy.integrate as integrate
from matplotlib.colors import LogNorm

subprocess.run(["make", "dipole"])

params = open("dipole_data/parameters.txt").readlines();
b           = float(params[0]) # magnetic field normalized in accoradance with radiation.hpp
m           =   int(params[1]) # number of dots in theta range
n           =   int(params[2]) # number of dots in omega range
gamma       = float(params[3]) # particle energy
omega_left  = float(params[4])
omega_right = float(params[5])
width_L     = float(params[6])
omega_L     = float(params[7])
theta_right = float(params[8])
theta_left  = float(params[9])
domega = (omega_right - omega_left) / n
chi = gamma * b

omega = np.arange(omega_left, omega_right, domega)
# theta = np.arange(theta_right, theta_left, m)

im1 = np.loadtxt("dipole_data/vacuum.txt") # data with ri = 1
im1 = np.reshape(im1, (n, m))
im2 = np.loadtxt("dipole_data/ref_ind.txt") # data with ri = 1
im2 = np.reshape(im1, (n, m))
ch_angle = np.loadtxt("dipole_data/ch_angle.txt") # cherenkov angle as a function of omega

cmap = 'hot'

fig = plt.figure(figsize = (6, 4.2), constrained_layout = True)
gs = fig.add_gridspec(2, 2)

#================= ax1, left plots ======================
ax1 = fig.add_subplot(gs[0, 0])
vmax = np.max(im1)
im1 = im1 / vmax
im_1 = ax1.imshow( im1, origin = 'lower'
                 , cmap = cmap, aspect = 'auto'
                 , extent = [omega_left/omega_L, omega_right/omega_L, -theta_left, theta_right]
                 , vmin = 0
                 , vmax = 0.1 * np.max(im1)
                 #, norm = LogNorm(vmin = 0.001, vmax = 0.01)
                 )

ax1.set_xlabel("ω/ω_L")
ax1.set_ylabel("θ (rad)")
#================= ax2, right plots ====================
ax2 = fig.add_subplot(gs[0,1])
ax2.plot(omega, ch_angle) 
vmax = np.max(im2)
im2 = im2 / vmax
im_2 = ax2.imshow( im2, origin = 'lower'
                 , cmap = cmap, aspect = 'auto'
                 , extent = [omega_left/omega_L, omega_right/omega_L, -theta_left, theta_right]
                 , vmin = 0
                 , vmax = 0.1 * np.max(im2)
                 )
ax2.set_xlabel("ω/ω_L")
#================= ax3, field float ====================
ax3 = fig.add_subplot(gs[1,0])
sum_im1 = np.sum(im1, axis = 1)
sum_im1 = sum_im1 / np.max(sum_im1)
sum_im2 = np.sum(im2, axis = 1) 
sum_im2 = sum_im2 / np.max(sum_im2)
ax3.plot(omega, sum_im1, "r", label="no-ri")
ax3.plot(omega, sum_im2, "g", label="ri")
ax3.legend(loc="upper left", fontsize = 15)

plt.show()
# plt.savefig('sinus.pdf')
