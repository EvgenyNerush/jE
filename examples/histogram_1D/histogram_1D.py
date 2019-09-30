import subprocess
import numpy             as np
import numpy.ctypeslib   as ctl
import ctypes            as ct
import matplotlib.pyplot as plt

# This example demonstrates that it is preferable to plot distribution functions with bins of
# non-constant width. Here 10 bins in both histogram_1D and plt.hist functions are used. Note that
# approximately 200-500 bins are requared by plt.hist to reproduce the picture given by
# histogram_1D with 10 bins. However, if one uses in plt.hist 200 bins of equal width, the
# resulting dN/dx would be noizy.

# Point generation #

subprocess.run(["g++", "-O3", "-fPIC", "-shared", "-Wall", "-Wpedantic", "lib.cpp", "-o",
    "lib.so"])
lib = ctl.load_library("lib.so", "./")

lib.get_n.restype  = ct.c_size_t
lib.get_nb.restype = ct.c_size_t
lib.points.restype = ct.POINTER(ct.c_double)
lib.bins.restype   = ct.POINTER(ct.c_double)

n  = lib.get_n()  # number of points
nb = lib.get_nb() # number of bins
xs_pointer = lib.points()
bs_pointer = lib.bins()

xs = ctl.as_array(xs_pointer, shape = (n,))  # coordinates of the points generated with Metropolis
                                            # algorithm
bs = ctl.as_array(bs_pointer, shape = (nb,)) # coordinates of the bin boundaries

# Distribution function computation #

hist_bins = 10 # number of bins for plt.hist which plots the number of particles per bin

x_vals = 0.5 * (bs[1:] + bs[:-1])
f_vals = n / ((nb - 1) * hist_bins * (bs[1:] - bs[:-1]))

z = (np.arange(500) + 1) / 500 # x values to plot the target distribution

plt.hist(xs, hist_bins, color = 'lightgreen', label = 'plt.hist')
plt.plot( z, 0.5 * n / hist_bins / np.sqrt(z) , '-', color = 'darkorange'
        , label = r'$1 / \sqrt{x}$')
plt.plot(x_vals, f_vals, '.', color = 'royalblue', label = 'histogram_1D')
plt.xlim(0, 1)
plt.ylim(0, 12.5 * n / hist_bins)
plt.legend()
plt.xlabel(r'$x$')
plt.ylabel(r'$dN/dx$, a.u.')
plt.savefig("out.png")
