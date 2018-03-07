# Find the find the distribution of limit cycle periods for different 
# total masses in a lattice with periodic boundary conditions 

import sys
import os 
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sandpiles import Sandpiles

# Parameters
sims_per_m = 100  # Number of simulations per lattice mass
L = 20  # Lattice size
maxheight = 4  # Maximun height before collapse

# Lattice masses that will be measured
msat = L*L*maxheight
#ms = np.linspace(msat*0.75, msat*1, 50, dtype=int)
ms = np.array((50, 70))

# Empty arrays to store the results
periods = np.zeros((ms.size, sims_per_m), dtype=int)
relaxtimes = np.zeros((ms.size, sims_per_m), dtype=int)

for m_idx, m in enumerate(ms):
    print("Mass: {0} ({1}/{2})".format(m, m_idx, ms.size))
    for j_sim in range(sims_per_m):
        # Generate random lattice
        latt = Sandpiles(L, L, maxheight=maxheight, pbc=True)
        latt.randomfill_mass(m)

        # Measure and store the limit cycle
        period, relaxtime = latt.findlimitcycle(maxtime=500)
        periods[m_idx, j_sim] = period
        relaxtimes[m_idx, j_sim] = relaxtime

# Save data
np.savetxt("lcycle_periodsL{0}N{1}.dat".format(L, sims_per_m), periods)
np.savetxt("lcycle_relaxtimesL{0}N{1}.dat".format(L, sims_per_m), relaxtimes)
np.savetxt("lcycle_msL{0}N{1}.dat".format(L, sims_per_m), ms)

# Analyse data 
periods_mean = np.mean(periods, 1)
periods_std = np.std(periods, 1)

periods_q25 = np.percentile(periods, 25, axis=1)
periods_q50 = np.percentile(periods, 50, axis=1)
periods_q75 = np.percentile(periods, 75, axis=1)

# Plot results
size = 3.5
fig, ax = plt.subplots(figsize=(1.62*size, size))

ax.plot(ms, periods_q25, "--", color="gray", label="25th percentile")
ax.plot(ms, periods_q50, "-", color="blue", label="50th percentile")
ax.plot(ms, periods_q75, "-.", color="gray", label="75th percentile")

ax.set_xlabel("Mass")
ax.set_ylabel("Limit cycle period")

lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), fancybox=True, 
          shadow=True)

fig.savefig("limitcycledistL{0}N{1}.pdf".format(L, sims_per_m),
            bbox_extra_artists=(lgd,), bbox_inches='tight')
