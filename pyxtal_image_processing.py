#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:35:41 2018

@author: mtrawick
"""

import numpy as np
import scipy.spatial as spat
import matplotlib.pyplot as plt


def do_triangulation(v):
    v.tri = spat.Delaunay(v.locations)
    v.triang = v.ax.triplot(v.locations[:,0], v.locations[:,1],
                           v.tri.simplices.copy(), #why the copy?
                         color='blue', zorder=1.5)
    v.imgCanvas.draw()
    #Still to do: figure out how to do this for periodic boundary conditions.
    
def do_circle_plot(v):
    v.circles = v.ax.scatter(v.locations[:,0], v.locations[:,1], 
                         color='green', zorder=2)
    v.imgCanvas.draw()

def do_raw_image(v):
    xsize, ysize = v.imgshape[0], v.imgshape[1]
    v.rawimg = v.ax.imshow(v.image, 
                              extent=[0, xsize, 0, ysize],
                              zorder=0,
                              cmap="gray")
    v.imgCanvas.draw()

