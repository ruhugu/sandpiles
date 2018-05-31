#-*- coding: utf-8 -*-
from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

from cellularautomata2d import CellAutomata2D
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors
from matplotlib import animation

class Sandpiles(CellAutomata2D):  

    def __init__(self, xlen, ylen, pbc=False, maxheight=4, cmap="Reds"):
        CellAutomata2D.__init__(self, xlen, ylen, pbc=pbc, dtype=int)
        self.maxheight = int(maxheight)

        # Create auxiliar lattices
        # Stores the modifications in the last step
        self._auxlatt = np.zeros((self._ylen_bc, self._xlen_bc),
                                 dtype=self.dtype)
        # Stores the cells bigger than maxheight in the last step
        self._collapse = np.zeros((self._ylen_bc, self._xlen_bc),
                                 dtype=self.dtype)
        # Create colormap for plots
        self.vmincolor = 0 - 0.5
        self.vmaxcolor = 2*maxheight + 0.5
        bounds = np.arange(0, 2*maxheight, 2*maxheight + 1) - 0.5
        self.cnorm = colors.BoundaryNorm(bounds, 256)
        self.cmap = plt.cm.get_cmap(cmap, 9)


    def randomfill(self, minval=0, maxval=None):
        """Fill the lattice randomly with values in the given range.

        """
        if maxval == None:
            maxval = self.maxheight

        for idx, height in np.ndenumerate(self.latt):
            self.latt[idx] = np.random.randint(minval, maxval + 1)

        return

    def randomfill_mass(self, mass):
        """Fill the lattice randomly up to a given total mass.

        """
        for i in range(mass):
            flat_idx = np.random.randint(0, self.size-1)
            self.latt.flat[flat_idx] += 1

        return


#    def mass(self):
#        """Return the value of the total mass of the system.
#
#        """
#        lattmass = self.latt.sum()
#        return lattmass


    def _evolvestep(self):
        """Evolve the system one step.

        Returns
        -------
            is_active : bool
                True if the lattice have moved and False otherwise.

        """
        self._auxlatt.fill(0)
        self._collapse[self._latt_idx] = (self.latt 
                                          > self.maxheight).astype(int)
        self._auxlatt -= self.maxheight*self._collapse
        self._auxlatt += np.roll(self._collapse, 1, axis=0)
        self._auxlatt += np.roll(self._collapse, -1, axis=0)
        self._auxlatt += np.roll(self._collapse, 1, axis=1)
        self._auxlatt += np.roll(self._collapse, -1, axis=1)
        self.latt += self._auxlatt[self._latt_idx]
        is_active = self._collapse.any()

        return is_active
        

#    def evolve(self, nsteps=0): 
#        """Evolve the system in nsteps timesteps.
#            
#        Parameters
#        ----------
#            nsteps : int
#                Number of steps the system will be evolved.
#
#            affectedcells : bool 
#                If True, a bool array with the affected cells will
#                be returned.
#                
#        Returns
#        -------
#            is_active : bool
#                True if the lattice have moved and False otherwise.
#
#        """
#        is_active = False
#        #collapse = np.zeros((self._ylen_bc, self._xlen_bc), dtype=self.dtype)
#        for i in range(nsteps):
#            is_active = self._evolvestep()
#
#            for (y, x), height in np.ndenumerate(self.latt):
#                if height > self.maxheight:
#                    self._auxlatt[self._bc(y, x)] -= 4
#                    self._auxlatt[self._bc(y, x+1)] += 1 
#                    self._auxlatt[self._bc(y, x-1)] += 1
#                    self._auxlatt[self._bc(y+1, x)] += 1
#                    self._auxlatt[self._bc(y-1, x)] += 1
#                    is_active = True
#
#        return is_active


    def measurecascade(self, maxtime=10000, retarray=False):
        """Measure the duration and number of affected cells of a cascade.

        Parameters
        ----------
            maxtime : int
                Maximum number of steps.

            retarray : bool
                If True, return a bool array with True values in the 
                cells affected by the cascade.

        Returns
        -------
            duration : int
                Duration of the cascade in steps. If the cascade does
                not end before mastime is reached, -1 is returned.

            cascadesize : int
                Number of cells affected by the cascade.

            cascadecells : bool array, optional
                Array with True in the cells affected by the cascade.

        """

        duration = 0
        # Array of cells in the cascade
        cascadecells = np.zeros((self._ylen_bc, self._xlen_bc), dtype=bool)
        active = True
        while active and (duration < maxtime):
            active = self._evolvestep()
            np.logical_or(cascadecells, (self._auxlatt != 0), cascadecells)
            duration += 1*int(active)

        cascadesize = (cascadecells.astype(int)).sum()

        # If system have not relaxed return -1 as duration
        if active:
            duration = -1

        output = duration, cascadesize

        if retarray == True:
            output.append(cascadecells)
            
        return output
        

    def animatecascade(self, maxtime=100000, steps_per_frame=1,
                       frame_interval=300):
        """Animate the evolution of a cascade.
        
        Parameters
        ----------
            nframes : int
                Number of frames in the animation.

            steps_per_frame : int
                Number of time steps advanced in each frame.

            frame_interval : float
                Interval between frames in miliseconds.

            colorsteps : int
                Number of steps between color changes.

        Returns
        -------
            anim

        """

        # Measure the cascade duration (up to maxtime)
        inilatt = np.copy(self.latt)
        duration, size = self.measurecascade(maxtime=maxtime)
        self.latt = inilatt

        # Array where True cells are those who belong to the cascade
        cascadearray = np.zeros((self.ylen, self.xlen), dtype=bool)
        # Array which stores the time at which cells where
        # affected by the cascade. If -1, the cell has not been
        # reached by the cascade yet.
        cascadetime = np.full((self.ylen, self.xlen), -1, dtype=int)

        # Create the color norm according to colorsteps
        # Since arange does not include the endpoint,
        # sum colorsteps to duration
#        colorboundaries = np.arange(0., duration + colorsteps, colorsteps)
#        ncolors = colorboundaries.size - 1
#        colornorm = colors.BoundaryNorm(colorboundaries, ncolors)

        # Create the colormap
        originalcmap = "jet"
        alpha = 0.6  # alpha of the colormap
        cmap_clist = (plt.get_cmap(originalcmap))(np.linspace(0, 1, 256))
        cmap_clist[:,-1][:] = alpha  # Set the alpha of all colors
        cmap_alpha = colors.ListedColormap(cmap_clist)
        cmap_alpha.set_under(alpha=0.)  # Set the alpha of values under 0

        # Plot configuration
        fig, ax = plt.subplots(figsize=(5,3))
        im = ax.imshow(self.latt, cmap="Greys", vmin=self.vmincolor, 
                       vmax=self.vmaxcolor, interpolation=None)
        im_cluster = ax.imshow(cascadetime, cmap=cmap_alpha,
                               vmin=-0.1, vmax=duration,
                               interpolation=None)
        # Color bars
        cbar_time = fig.colorbar(im_cluster, ax=ax)  # Time colorbar
        cbar_time.set_label("Time")
        cbar_heigth = fig.colorbar(im, ax=ax)  # Height colorbar
        cbar_heigth.set_label("Height")

        fig.tight_layout()

        nframes = duration
        
        # Update function
        def update(i, cascadearray, cascadetime, self, im, im_cluster):
            self._evolvestep()
            cascadetime += (i + 1)*np.logical_and(
                                np.logical_not(cascadearray),
                                self._auxlatt[self._latt_idx] != 0)
            np.logical_or(cascadearray, (self._auxlatt[self._latt_idx] != 0), 
                          cascadearray)
            im.set_array(self.latt)
            #im_cluster.set_array(cascadearray)
            im_cluster.set_array(cascadetime)
            #print(cascadearray)
            return im, im_cluster

        anim = animation.FuncAnimation(fig, update, frames=nframes, 
                                       blit=False, fargs=(cascadearray,
                                       cascadetime, self, im, im_cluster))
        return anim


    def scaleinvariance(self, maxtime=10000):
        # Array where True cells are those who belong to the cascade
        cascadearray = np.zeros((self.ylen, self.xlen), dtype=bool)
        # Array which stores the time at which cells where
        # affected by the cascade. If -1, the cell has not been
        # reached by the cascade yet.
        cascadetime = np.full((self.ylen, self.xlen), -1, dtype=int)

        time = 0
        active = True
        while active and (time < maxtime):
            time += 1
            active = self._evolvestep()
            cascadetime += (time + 1)*np.logical_and(
                                np.logical_not(cascadearray),
                                self._auxlatt[self._latt_idx] != 0)
            np.logical_or(cascadearray, (self._auxlatt[self._latt_idx] != 0), 
                          cascadearray)

        return cascadetime

    def plotcascadetime(self, cascadetime, region=None):

        if region == None:
            region = self._latt_idx

        # Create the colormap
        originalcmap = "jet"
        alpha = 0.6  # alpha of the colormap
        cmap_clist = (plt.get_cmap(originalcmap))(np.linspace(0, 1, 256))
        cmap_clist[:,-1][:] = alpha  # Set the alpha of all colors
        cmap_alpha = colors.ListedColormap(cmap_clist)
        cmap_alpha.set_under(alpha=0.)  # Set the alpha of values under 0

        # Plot configuration
        fig, ax = plt.subplots(figsize=(5,3))
        im = ax.imshow(self.latt[region], cmap="Greys", vmin=self.vmincolor, 
                       vmax=self.vmaxcolor, interpolation=None)
        im_cluster = ax.imshow(cascadetime[region], cmap=cmap_alpha,
                               vmin=-0.1, #vmax=duration,
                               interpolation=None)
        # Color bars
        cbar_time = fig.colorbar(im_cluster, ax=ax)  # Time colorbar
        cbar_time.set_label("Time")
        cbar_heigth = fig.colorbar(im, ax=ax)  # Height colorbar
        cbar_heigth.set_label("Height")

        fig.tight_layout()



#    def findlimitcycle(self, maxtime=50):
#        history = np.zeros((maxtime+1, self.ylen, self.xlen))
#
#        foundcycle = False
#        relaxtime = -1
#        period = -1  # This will be returned if no cycle is found
#
#        j_step = 0
#        while (not foundcycle) and j_step <= maxtime:
#            history[j_step] = self.latt 
#            self.evolve(1)
#            j_step += 1
#
#            for i in reversed(range(j_step)):
#                if (history[i] == self.latt).all(): 
#                    foundcycle = True
#                    period = j_step - i
#                    relaxtime = j_step
#                    break
#
#        return period, relaxtime
