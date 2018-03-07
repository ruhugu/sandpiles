# -*- coding: utf-8 -*-
from __future__ import (print_function, division, 
                        absolute_import, unicode_literals)

# TODO Cambiar línea para no usar *
from celullarautomata2d import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors

#TODO:add custom thresold
class Sandpiles(CelAutomata2D):  
    maxheight = 4
    def __init__(self, xlen, ylen, pbc=False, maxheight=4):
        CelAutomata2D.__init__(self, xlen, ylen, pbc=pbc, dtype=int)
        self.maxheight = int(maxheight)

        # Create auxiliar lattices
        # Stores the modifications in the last step
        self._auxlatt = np.zeros((self._xlen_bc, self._ylen_bc),
                                 dtype=self.dtype)
        # Stores the cells bigger than maxheight in the last step
        self._collapse = np.zeros((self._xlen_bc, self._ylen_bc),
                                 dtype=self.dtype)
        # Create colormap
        self.vmincolor = 0 - 0.5
        self.vmaxcolor = 2*maxheight + 0.5
        bounds = np.arange(0, 2*maxheight, 2*maxheight + 1) - 0.5
        self.cnorm = colors.BoundaryNorm(bounds, 256)
        self.cmap = plt.cm.get_cmap('Reds', 9)


    def randomfill(self, minval=0, maxval=None):
        """Fill the lattice randomly with values in the given range.

        """
        if maxval == None:
            maxval = self.maxheight

        for idx, height in np.ndenumerate(self.latt):
            self.latt[idx] = np.random.randint(self.maxheight, self.maxheight*5)

        return


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
        

    def evolve(self, nsteps=1, ret_activecells=True): 
        """Evolve the system in nsteps timesteps.
            
        Parameters
        ----------
            nsteps : int
                Number of steps the system will be evolved.

            affectedcells : bool 
                If True, a bool array with the affected cells will
                be returned.
                
        Returns
        -------
            is_active : bool
                True if the lattice have moved and False otherwise.

        """
        if ret_activecells:
            activecells = np.zeros((self._xlen_bc, self._ylen_bc), dtype=bool)

        is_active = False
        #collapse = np.zeros((self._xlen_bc, self._ylen_bc), dtype=self.dtype)
        for i in range(nsteps):
            is_active = self._evolvestep()

#            for (x, y), height in np.ndenumerate(self.latt):
#                if height > self.maxheight:
#                    self._auxlatt[self._bc(x, y)] -= 4
#                    self._auxlatt[self._bc(x, y+1)] += 1 
#                    self._auxlatt[self._bc(x, y-1)] += 1
#                    self._auxlatt[self._bc(x+1, y)] += 1
#                    self._auxlatt[self._bc(x-1, y)] += 1
#                    is_active = True

        return is_active


    def measurecascade(self, maxtime=10000):
        """Measure the duration and number of affected cells of a cascade.

        Parameters
        ----------
            maxtime : int
                Maximum number of steps.

        Returns
        -------
            duration : int
                Duration of the cascade in steps. If the cascade does
                not end before mastime is reached, -1 is returned.

            cascadesize : int
                Number of cells affected by the cascade.

            cascadecells : bool array
                Array with True in the cells affected by the cascade.
        """

        duration = 0
        # Array of cells in the cascade
        cascadecells = np.zeros((self._xlen_bc, self._ylen_bc), dtype=bool)
        active = True
        while active and (duration < maxtime):
            active = self._evolvestep()
            np.logical_or(cascadecells, (self._auxlatt != 0), cascadecells)
            duration += 1*int(active)

        cascadesize = (cascadecells.astype(int)).sum()

        # If system have not relaxed return -1 as duration
        if active:
            duration = -1

        return duration, cascadesize, cascadecells
        

    def animatecascade(self, maxtime=10000, steps_per_frame=1, frame_interval=300):
        """Animate the evolution of a cascade.
        
        Parameters
        ----------
            nframes : int
                Number of frames in the animation.

            steps_per_frame : int
                Number of time steps advanced in each frame.

            frame_interval : float
                Interval between frames in miliseconds.

        Returns
        -------
            anim

        """

        # Measure the cascade duration (up to maxtime)
        inilatt = np.copy(self.latt)
        duration, size, cascadearray = self.measurecascade(maxtime=maxtime)
        self.latt = inilatt

        # Initiliaze array of cells in the cascade
        cascadearray = np.zeros((self.xlen, self.ylen), dtype=bool)

        # Update function
        def update(i, cascadearray, self, im, im_cluster):
            self._evolvestep()
            np.logical_or(cascadearray, (self._auxlatt[self._latt_idx] != 0), 
                          cascadearray)
            im.set_array(self.latt)
            im_cluster.set_array(cascadearray)
            return im

        # Plot configuration
        fig, ax = plt.subplots()
        
        ncolors = 2
        cmap_alpha = (plt.get_cmap("Wistia_r"))(np.arange(ncolors))
        cmap_alpha[:,-1] = np.array((0., 0.6))
        cmap_alpha = colors.ListedColormap(cmap_alpha)
        im = ax.imshow(self.latt, cmap="Greys", vmin=self.vmincolor, 
                       vmax=self.vmaxcolor, interpolation=None)
        im_cluster = ax.imshow(cascadearray, cmap=cmap_alpha,
                               vmin=0, vmax=1, interpolation=None)
        cbar = fig.colorbar(im, ax=ax)

        nframes = duration
        anim = animation.FuncAnimation(fig, update, frames=nframes, 
                                       blit=False, fargs=(cascadearray, self,
                                       im, im_cluster))
        return anim


#    def findlimitcycle(self, maxtime=50):
#        history = np.zeros((maxtime+1, self.xlen, self.ylen))
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
