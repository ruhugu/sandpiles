# Save an animation of the evolution of a perturbation in lattice
# with a critical configuration.

import os
import sys
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sandpiles import Sandpiles

# Parameters
L = 50
perturbation = 5
fps = 10
# Perturbation position
y_pert = 25
x_pert = 25

# Generate the lattice configuration
latt = Sandpiles(L, L)
latt.randomfill(minval=latt.maxheight, maxval=latt.maxheight*5)
latt.relax()
latt.latt[y_pert,x_pert] += perturbation

# Create the animation
anim = latt.animatecascade()

# Show the animation
plt.show()

# Save animation to file (this does not work if plt.show() has been used before)
#anim.save("clusterevolutionL{0}.mp4".format(L), dpi=300, fps=fps,
#          extra_args=['-vcodec', 'libx264'])
