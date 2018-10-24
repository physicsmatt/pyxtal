#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:35:41 2018

@author: mtrawick
"""

import numpy as np
import scipy.spatial as spat
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib


def do_raw_image(v):
    xsize, ysize = v.imgshape[0], v.imgshape[1]
    v.rawimg = v.ax.imshow(v.image, 
                              extent=[0, xsize, 0, ysize],
                              zorder=0,
                              cmap="gray")
    v.imgCanvas.draw()
    
    
def do_circle_plot(v):
    v.circles = v.ax.scatter(v.locations[:,0], v.locations[:,1], 
                         color='green', zorder=3, facecolors='none')
    v.imgCanvas.draw()


def do_triangulation(v):
    #Performs Delaunay triangulation and creates a plot of it.
    #Also finds outermost vertices, which we'll want later.
    v.tri = spat.Delaunay(v.locations, qhull_options="Qj")
#    triangulation = v.ax.triplot(v.locations[:,0], v.locations[:,1],
#                           v.tri.simplices.copy(), #why the copy?
#                         color='blue', zorder=4)
    
#    triangulation.set_visible(0)
#    v.triang = triangulation
    v.imgCanvas.draw()
    
    #Now find outermost vertices.  
    #tri.neighbors gives the neighbors of each triangle ("simplex").
    #Index of -1 means no neighboring simplex, so the other two vertices
    #in that triangle must be on the outside.
    where_outer_tri = np.where(v.tri.neighbors == -1)[0] #rows that contain -1
    outer_tri = v.tri.simplices[where_outer_tri,:].reshape(-1)
    outer_neighbors = v.tri.neighbors[where_outer_tri,:].reshape(-1)
    v.tri.outer_vertices = np.unique(outer_tri[np.where(outer_neighbors != -1)])
    #above are indices of outer vertices
    
#not sure what I was thinking here.    
#    v.tri.outermost = np.unique(outer_tri.reshape(-1)) #indices of outer vertices


#    v.plt = SimpleNamespace()
#I would like to groupd plots togeher, like v.plt.outermost, but I don't know how.
    v.outermost = v.ax.scatter(v.locations[v.tri.outer_vertices,0], 
                               v.locations[v.tri.outer_vertices,1], 
                         color='red', zorder=5)
    v.imgCanvas.draw()
   
    
    #Still to do: figure out how to do this for periodic boundary conditions.

    
#def do_disclinations(v):
#Not sure if it's best to do these as a single plot with multiple colors,
#or multiple plots for each different flavor of disclination.  Wait to see how
#hard it is to handle zooming first.    

def do_angle_field(v):
    #This function ultimately produces a color image representing the local
    #orientation of the crystal.

    #First, get the x,y positions and angles of each "bond" between vertices
    num_tri = len(v.tri.simplices)
    num_outer = len(v.tri.outer_vertices)
    num_bonds = int((num_tri * 3 + num_outer ) / 2)
    v.tri.bondsx = np.zeros(num_bonds)
    v.tri.bondsy = np.zeros(num_bonds)
    v.tri.bondsl = np.zeros(num_bonds)
    v.tri.bondsangle = np.zeros(num_bonds)
    #There's probably a cute way to do this using np.where, but since there's
    #only order N of these to go through, I'm going to iterate instead.
    indptr = v.tri.vertex_neighbor_vertices[0]
    indices = v.tri.vertex_neighbor_vertices[1]
    bondi = 0
    #See https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.vertex_neighbor_vertices.html#scipy.spatial.Delaunay.vertex_neighbor_vertices
    for v1i in range(0,len(v.tri.points)):        
        for v2i in indices[indptr[v1i]:indptr[v1i + 1]]:
            if v2i > v1i: #This eliminates double counting of bonds
                x1 = v.tri.points[v1i,0]
                x2 = v.tri.points[v2i,0]
                y1 = v.tri.points[v1i,1]
                y2 = v.tri.points[v2i,1]
                v.tri.bondsx[bondi] = (x1 + x2) / 2
                v.tri.bondsy[bondi] = (y1 + y2) / 2
                v.tri.bondsl[bondi] =  np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                angle = np.pi + np.arctan2(y2-y1, x2-x1) #from 0 to 2pi
                #turn this into 0 to 60 degrees:
                v.tri.bondsangle[bondi] = angle % (np.pi / 3)
                bondi += 1

    #Now, create the orientation field:
    xsize, ysize = v.imgshape[0], v.imgshape[1]
    grid_x, grid_y = np.mgrid[0:xsize:1, 0:ysize:1]
    cosangle = np.cos(6 * v.tri.bondsangle)
    sinangle = np.sin(6 * v.tri.bondsangle)
    cosarr = scipy.interpolate.griddata((v.tri.bondsx, v.tri.bondsy), 
                                          cosangle, (grid_x, grid_y),method='linear')
    sinarr = scipy.interpolate.griddata((v.tri.bondsx, v.tri.bondsy), 
                                          sinangle, (grid_x, grid_y),method='linear')
    anglearr = (np.arctan2(sinarr,cosarr) + np.pi ) / (2 * np.pi)


    hsvimg = np.ones((xsize, ysize, 3))
    hsvimg[:,:,0] = anglearr
    v.rgbimg = matplotlib.colors.hsv_to_rgb(hsvimg)
    v.angleimg = v.ax.imshow(v.rgbimg, 
                              extent=[0, xsize, 0, ysize],
                              zorder=1, alpha = 0.7)
    v.imgCanvas.draw()
    
    
    



if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtalmain
    pyxtalmain.vp_start_gui()
