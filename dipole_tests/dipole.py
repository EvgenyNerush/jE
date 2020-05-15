import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import matplotlib.colors as colors
import matplotlib.cm as mcm
import copy
import scipy.integrate as integrate

subprocess.run(["make", "sinus"])

params = open("sinus_data/parameters.txt").readlines();
b           = float(params[0]) # magnetic field normalized in accoradance with radiation.hpp
theta_b     = float(params[0]) # theta boundery
m           =   int(params[2]) # number of dots in theta range
n           =   int(params[3]) # number of dots in omega range
gamma       = float(params[4]) # particle energy
omega_left  = float(params[5])
omega_right = float(params[6])
omega_L     = float(params[7])
width_L     = float(params[8])
domega = (omega_right - omega_left) / n
chi = gamma * b
lw = 1
ms = 3

im1 = np.loadtxt("dipole_data/vacuum.txt") # data with ri = 1
im1 = np.reshape(im1, (n, m))

plt.style.use("plasma.mplstyle")
cmap = 'Greys'

fig = plt.figure(figsize = (6, 4.2), constrained_layout = True)
gs = fig.add_gridspec(1, 1)

#================= ax1, left plots ======================
ax1 = fig.add_subplot(gs[:, 0])
vmax = np.max(im1)
im12 = np.empty((2 * m, n))
im12[0:m] = im1.transpose()
im12[m:(2 * m)] = im2.transpose()
im12 = im12 / vmax
im_1 = ax1.imshow( im12, origin = 'lower'
                 , cmap = cmap, aspect = 'auto'
                 , extent = [omega_left * b / gamma, omega_right * b / gamma, -theta_b_mrad,
                     theta_b_mrad]
                 , vmin = 0
                 , vmax = 1
                 )

ax1.set_xlabel("ℏω/ε")
ax1.set_xlim([0, 1])
ax1.set_xticks([0, 0.5, 1])
ax1.set_xticklabels(["0", "0.5", "1"])
ax1.set_ylabel("θ (rad)")
ax1.set_ylim([0, theta_b])
#cb1 = fig.colorbar( im_1
#                  , orientation = 'horizontal'
#                  , aspect = 9
#                  , ticks = [0, 0.5, 1]
#                  , format = '%.1g'
#                  )
#cb1.set_label("d²I/dωdθ (a.u.)")

plt.savefig('sinus.pdf')
