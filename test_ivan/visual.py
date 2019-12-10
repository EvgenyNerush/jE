import numpy as np
import matplotlib.pyplot as plt

im1 = np.loadtxt("vv.txt", delimiter = '_')
im2 = np.loadtxt("vp.txt", delimiter = '_')
im3 = np.loadtxt("mp.txt", delimiter = '_')
im4 = np.loadtxt("mv.txt", delimiter = '_')
xy = np.loadtxt("od.txt", delimiter = '_')
theta_list = xy[0]
omega_list = xy[1]

fig, axs = plt.subplots(2,2)
img1 = axs[0,0].imshow(im1, cmap='hot', extent = [-1,1,-1,1])
axs[0,0].set_title(r'$\alpha \chi^{2/3} = 1; n = 1$')
axs[0,0].set_xticks(np.linspace(-1,1,np.size(theta_list)))
axs[0,0].set_yticks(np.linspace(-1,1,np.size(omega_list)))
axs[0,0].set_xticklabels(theta_list)
axs[0,0].set_yticklabels(omega_list)

img2 = axs[0,1].imshow(im2, cmap='hot', extent = [-1,1,-1,1])
axs[0,1].set_title(r'$\alpha \chi^{2/3} \approx 0.05; n = 1$')
axs[0,1].set_xticks(np.linspace(-1,1,np.size(theta_list)))
axs[0,1].set_yticks(np.linspace(-1,1,np.size(omega_list)))
axs[0,1].set_xticklabels(theta_list)
axs[0,1].set_yticklabels(omega_list)

img3 = axs[1,0].imshow(im3, cmap='hot', extent = [-1,1,-1,1])
axs[1,0].set_title(r'$\alpha \chi^{2/3} \approx 0.05; n \neq 1$')
axs[1,0].set_xticks(np.linspace(-1,1,np.size(theta_list)))
axs[1,0].set_yticks(np.linspace(-1,1,np.size(omega_list)))
axs[1,0].set_xticklabels(theta_list)
axs[1,0].set_yticklabels(omega_list)

img4 = axs[1,1].imshow(im4, cmap='hot', extent = [-1,1,-1,1])
axs[1,1].set_title(r'$\alpha \chi^{2/3} = 1; n \neq 1$')
axs[1,1].set_xticks(np.linspace(-1,1,np.size(theta_list)))
axs[1,1].set_yticks(np.linspace(-1,1,np.size(omega_list)))
axs[1,1].set_xticklabels(theta_list)
axs[1,1].set_yticklabels(omega_list)

plt.show()
