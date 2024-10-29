import matplotlib.pyplot as plt
# import cmocean.cm as cmo
from matplotlib import colormaps
import numpy as np
import subprocess


# Reconstruct the free energy from HILLS using the plumed tool "sum_hills" :
p = subprocess.Popen("plumed sum_hills --hills HILLS --outfile fes.dat" +
                     " --bin 300,300 --min 0.3,-0.35 --max 1.2,1.56",
                     shell=True, stdout=subprocess.PIPE)
p.wait()

# Import free energy and reshape with the number of bins defined in the
# reconstruction process.
scm = np.loadtxt('fes.dat', usecols=0).reshape(301, 301)
tcm = np.loadtxt('fes.dat', usecols=1).reshape(301, 301)
fes = np.loadtxt('fes.dat', usecols=2).reshape(301, 301)

# Plot
fig, ax = plt.subplots(figsize=(10, 9))

# Plot free energy surface
im = ax.contourf(scm, tcm, fes, 10, cmap=colormaps['Blues_r'])# cmo.tempo_r)
cp = ax.contour(scm, tcm, fes, 10,
                linestyles='-', colors='darkgray', linewidths=1.2)

# Plot parameters
ax.set_xlabel('SCM', fontsize=40)
ax.set_ylabel('TCM', fontsize=40)
ax.tick_params(axis='y', labelsize=25)
ax.tick_params(axis='x', labelsize=25)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(label=r'FES[$\epsilon$]', fontsize=40)
cbar.ax.tick_params(labelsize=32)

plt.tight_layout()
plt.show()
