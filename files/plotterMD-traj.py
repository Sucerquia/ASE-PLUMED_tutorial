import matplotlib.pyplot as plt
import numpy as np


# Import free energy and reshape with the number of bins defined in the
# reconstruction process.
scm = np.loadtxt('COLVAR', usecols=1)
tcm = np.loadtxt('COLVAR', usecols=2)

# Plot
fig, ax = plt.subplots(figsize=(10, 9))

# Plot free energy surface
plt.plot(scm, tcm, 'o')
# Plot parameters
ax.set_xlabel('SCM', fontsize=40)
ax.set_ylabel('TCM', fontsize=40)
ax.tick_params(axis='y', labelsize=25)
ax.tick_params(axis='x', labelsize=25)

plt.tight_layout()
plt.show()
