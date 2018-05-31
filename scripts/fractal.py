# Create fractal figure

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sandpiles import Sandpiles

# Parameters
L = 50
mass_filledregion = 100000
size_filledregion = 2 

# Define the filled region according to the parameters
filledregion = np.zeros((L,L), dtype=bool)
filledregion[
        L//2 - size_filledregion//2 : L//2 + size_filledregion//2,
        L//2 - size_filledregion//2 : L//2 + size_filledregion//2] = True

# Create and initialize the lattice
latt = Sandpiles(L, L)
latt.latt[filledregion] = mass_filledregion

# Let the lattice relax
latt.relax(maxtime=10*mass_filledregion, stepsize=100)

# Plot the result
size = 10
fig, ax = plt.subplots(figsize=(size,size))
im = ax.imshow(latt.latt, cmap=latt.cmap, interpolation=None)
cbar = fig.colorbar(im, ax=ax)

# Save the plot to a file
fig.savefig("fractalL{0}M{1}s{2}.png".format(L, mass_filledregion, size_filledregion))

# Show the plot
#plt.show()

