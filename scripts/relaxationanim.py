# Generate an animation of the relaxation process to a
# critical configuration.

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sandpiles import Sandpiles

# Parameters 
L = 50  # Lattice length

# Generate the lattice randomly with heights between maxheight and
# maxrandomval.
latt = Sandpiles(L, L)
maxrandomval = 5*latt.maxheight
latt.randomfill(minval=latt.maxheight, maxval=5*latt.maxheight)

# Find the relaxation time in order to know the number
# of frames in the animation.
inilatt = np.copy(latt.latt)
relaxtime = latt.relax()
latt.latt = inilatt

# Create a new colormap to take into account values bigger
# than maxheight
latt.vmaxcolor = maxrandomval + 0.5
latt.cmap = plt.cm.get_cmap('Reds', maxrandomval + 1)

# Generate the animation
anim = latt.animate(relaxtime)

# Show the animation
plt.show()

# Save animation to file (this does not work if plt.show() has been used before)
#anim.save("relaxationL{0}.mp4".format(L), dpi=300, fps=20,
#          extra_args=['-vcodec', 'libx264'])

