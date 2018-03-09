# Generate critical configurations and measure the size and 
# duration of the cascades that appear when the lattice is 
# perturbed locally.

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sandpiles import Sandpiles

# Parameters
L = 20  # Lattice length
perturbation = 5  # Local perturbation
nsim = 2500  # Number of experiments

# Vectors to store the data
durations = np.zeros(nsim*L*L)
cascadesizes = np.zeros(nsim*L*L)

for j_sim in range(nsim):
    print("Iteration: {0}/{1}".format(j_sim, nsim-1))
    # Generate the lattice randomly with heights bigger or
    # equal to maxheight
    latt = Sandpiles(L, L, pbc=False)
    latt.randomfill(minval=latt.maxheight, maxval=latt.maxheight*5)

    # Let it relax
    relaxstep = 100
    latt.relax(stepsize=relaxstep)

    # Store the relaxed lattice in order to recover it after
    # each cascade.
    relaxedlatt = np.copy(latt.latt)

    # Perturbe each cell of the lattice to measure cascade 
    # sizes and durations
    for cell_idx, height in enumerate(latt.latt.flat):
        latt.latt = relaxedlatt
        latt.latt.flat[cell_idx] += perturbation
        duration, cascadesize = latt.measurecascade()
        durations[j_sim*(latt.size-1) + cell_idx] = duration
        cascadesizes[j_sim*(latt.size-1) + cell_idx] = cascadesize
 
# Save results
np.savetxt("cascadedurationsL{0}N{1}.dat".format(L, nsim), durations)
np.savetxt("cascadesizesL{0}N{1}.dat".format(L, nsim), cascadesizes)

# Plot the results
# Duration plot
size = 3.5
fig, ax = plt.subplots(figsize=(1.62*size, size))
histo, bin_edges = np.histogram(durations, bins=30)
bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
plt.loglog(bin_centers, histo)
ax.set_xlabel("Cascade duration")
ax.set_ylabel("Frequency")
fig.tight_layout()
fig.savefig("cascadedurationL{0}N{1}.png".format(L, nsim))

# Size plot
size = 3.5
fig, ax = plt.subplots(figsize=(1.62*size, size))
histo, bin_edges = np.histogram(cascadesizes, bins=30)
bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
plt.loglog(bin_centers, histo)
ax.set_xlabel("Cascade size")
ax.set_ylabel("Frequency")
fig.tight_layout()
fig.savefig("cascadesizeL{0}N{1}.png".format(L, nsim))
