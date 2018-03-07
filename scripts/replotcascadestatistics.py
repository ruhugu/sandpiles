# Import the results from different runs of "cascadestatistics.py"
# and plot then together.

import numpy as np
import itertools 
from matplotlib import pyplot as plt

# Parameters of the data to import
# Lattice lengths
Ls= (20, 50, 100)
# Number of experiments
nsims = (1250, 200, 50)

# Generate marker vector
markers1 = itertools.cycle(("o", "s", "^", "v"))
markers2 = itertools.cycle(("o", "s", "^", "v"))

# Plots
size = 3.5
fig_time, ax_time = plt.subplots(figsize=(1.62*size, size))
fig_size, ax_size = plt.subplots(figsize=(1.62*size, size))

for L, nsim in zip(Ls, nsims):
    # Load data
    relaxtimes = np.loadtxt("cascadedurationsL{0}N{1}.dat".format(L, nsim)) + 1
    sizes = np.loadtxt("cascadesizesL{0}N{1}.dat".format(L, nsim)) + 1

    # Calculate number of cascades
    nexp = nsim*L*L

    # Create logarithmic bins
    maxtime_log = np.log10(np.amax(relaxtimes))
    mintime_log = np.log10(np.amin(relaxtimes))
    maxsize_log = np.log10(np.amax(sizes))
    minsize_log = np.log10(np.amin(sizes))
    
    bin_edges_time = np.power(10, np.linspace(1, maxtime_log, 30))
    bin_edges_size = np.power(10, np.linspace(1, maxsize_log, 30))

    bin_step_time = (bin_edges_time[1:] - bin_edges_time[:-1])
    bin_step_size = (bin_edges_size[1:] - bin_edges_size[:-1])

    # Create histograms
    histo_time, bin_edges_time = np.histogram(relaxtimes,
            bins=bin_edges_time)
    bin_centers_time = (bin_edges_time[1:] + bin_edges_time[:-1])/2

    histo_size, bin_edges_size = np.histogram(sizes,
            bins=bin_edges_size)
    bin_centers_size = (bin_edges_size[1:] + bin_edges_size[:-1])/2

    # Plot the data
    ax_time.loglog(bin_centers_time, histo_time/(bin_step_time*float(nexp)),
                   linestyle="", marker=markers1.next(),
                   label="L{0}, N{1}".format(L, nsim))
    ax_size.loglog(bin_centers_size, histo_size/(bin_step_size*float(nexp)),
                   linestyle="", marker=markers2.next(),
                   label="L{0}, N{1}".format(L, nsim))


# Set axis labels and legends
ax_time.set_xlabel("Cascade duration")
ax_time.set_ylabel("Rel. frequency")
ax_size.set_xlabel("Cascade size")
ax_size.set_ylabel("Rel. frequency")
fig_time.legend()
fig_time.tight_layout()
fig_size.legend()
fig_size.tight_layout()

# Save plots to file
fig_time.savefig("cascadeduration.pdf")
fig_size.savefig("cascadesize.pdf")
