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


def do_inverted_images(v):
    v.inv_image = np.max(v.image) - v.image
    
    
def do_circle_plot(v):
    #Plot circles around the locations of the spheres.
    #I've used a collection of circle "patches", so that the radius can
    #scale with the axes as the figure is zoomed.
    #I'll still have to rescale the linewidth manually, however.
    

    radius = int(v.pmw.sphereSize[0]*0.7)
    patches = [matplotlib.patches.CirclePolygon(xy, radius) 
                    for xy in v.locations]
    coll = matplotlib.collections.PatchCollection(patches, 
                            edgecolor=('#B4FF64'), facecolor='None', zorder=2)
    #use RGB color (180,255,100), or #B4FF64.  The hex value seems much faster!
    v.circles = v.ax.add_collection(coll)

    #Note that this is how you set linewidths:
    #v.circles.set_linewidth(2)
    v.imgCanvas.draw()


def do_triangulation(v):
    #Performs Delaunay triangulation and creates a plot of it.
    #Also finds outermost vertices, which we'll want later.
    #Also gets orientation of each bond, which we'll use later for angle field.

    v.tri = spat.Delaunay(v.locations, qhull_options="Qj")
    #v.triang = v.ax.triplot(v.locations[:,0], v.locations[:,1],
    #                       v.tri.simplices.copy(), #why the copy?
    #                     color='blue', zorder=4)
    #The line above worked, but the object it created included some extra 
    #unneeded stuff, and the object couldn't be used with set_visible method.
    
    #Now find outermost vertices of the triangulation.  
    #tri.neighbors gives the neighbors of each triangle ("simplex").
    #Index of -1 means no neighboring simplex, so the other two vertices
    #in that triangle must be on the outside.
    where_outer_tri = np.where(v.tri.neighbors == -1)[0] #rows that contain -1
    outer_tri = v.tri.simplices[where_outer_tri,:].reshape(-1)
    outer_neighbors = v.tri.neighbors[where_outer_tri,:].reshape(-1)
    v.tri.outer_vertices = np.unique(outer_tri[np.where(outer_neighbors != -1)])
    #above are the indices of the outer vertices
#    v.outermost = v.ax.scatter(v.locations[v.tri.outer_vertices,0], 
#                               v.locations[v.tri.outer_vertices,1], 
#                         color='red', zorder=5)
    #I would like to figure out how to group all of the plots together, 
    #like v.plt.outermost and v.plt.circles but I don't know how.
    v.imgCanvas.draw()
    #Still to do: figure out how to do this for periodic boundary conditions.

    #Now step through and get bond information.
    #First, get the x,y positions and angles of each "bond" between vertices
    num_tri = len(v.tri.simplices)
    num_outer = len(v.tri.outer_vertices)
    num_bonds = int((num_tri * 3 + num_outer ) / 2)
    v.tri.bondsx = np.zeros(num_bonds)
    v.tri.bondsy = np.zeros(num_bonds)
    v.tri.bondsl = np.zeros(num_bonds)
    v.tri.bondsangle = np.zeros(num_bonds)
    v.tri.cnum = np.zeros(len(v.tri.points)) #coordination number of each vertex
    segs = np.zeros((num_bonds,2,2))
    #See https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.vertex_neighbor_vertices.html#scipy.spatial.Delaunay.vertex_neighbor_vertices
    indptr = v.tri.vertex_neighbor_vertices[0]
    indices = v.tri.vertex_neighbor_vertices[1]
    #There's probably a cute way to do the next part using np.where, but since 
    #there's only order N of these to go through, I'll just iterate instead.
    bondi = 0
    for v1i in range(0,len(v.tri.points)):
        v.tri.cnum[v1i] = indptr[v1i+1] - indptr[v1i]
        for v2i in indices[indptr[v1i]:indptr[v1i + 1]]:
            if v2i > v1i: #This eliminates double counting of bonds
                x1 = v.tri.points[v1i,0]
                x2 = v.tri.points[v2i,0]
                y1 = v.tri.points[v1i,1]
                y2 = v.tri.points[v2i,1]
                segs[bondi] = np.array([[x1,y1],[x2,y2]])
                v.tri.bondsx[bondi] = (x1 + x2) / 2
                v.tri.bondsy[bondi] = (y1 + y2) / 2
                v.tri.bondsl[bondi] =  np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                angle = np.pi + np.arctan2(y2-y1, x2-x1) #from 0 to 2pi
                #turn this into 0 to 60 degrees:
                v.tri.bondsangle[bondi] = angle % (np.pi / 3)
                bondi += 1
    line_coll = matplotlib.collections.LineCollection(segs, color='blue', zorder=4)
    v.triang = v.ax.add_collection(line_coll)
    v.imgCanvas.draw()

def disc_color(cnum):
    if cnum <= 4:
        return('#FF00FF') #magenta
    if cnum == 5:
        return('#FF0000') #red
    if cnum == 7:
        return('#00FF00') #green
    if cnum >=8:
        return('#00FFFF') #cyan

    
def do_disclinations(v):

    #Note that for some reason, the np.where returns a tuple with one item,
    #so I need to access it with [0].
    v.disc = np.where(v.tri.cnum != 6)[0]
    radius = int(v.pmw.sphereSize[0]*0.5)

    #The line below should be appreciated, as it is particularly pythonic.
    patches = [matplotlib.patches.Circle(v.locations[i], radius, 
                                         facecolor=disc_color(v.tri.cnum[i])) 
                                    for i in v.disc]

    coll = matplotlib.collections.PatchCollection(patches, match_original=True, zorder=5)
    v.pltdisc = v.ax.add_collection(coll)
    v.imgCanvas.draw()



def do_angle_field(v):
    #This function ultimately produces a color image representing the local
    #orientation of the crystal.

#    v.bondsplot = v.ax.scatter(v.tri.bondsx, 
#                               v.tri.bondsy, 
#                         color='red', zorder=15)

    #Now, create the orientation field:
    xsize, ysize = v.imgshape[0], v.imgshape[1]
    grid_x, grid_y = np.mgrid[0:xsize:1, 0:ysize:1]
    cosangle = np.cos(6 * (v.tri.bondsangle - np.pi/6))
    sinangle = np.sin(6 * (v.tri.bondsangle - np.pi/6))
    
    #Note below that I have switched grid_x and grid_y, contrary to the example
    #shown in the numpy reference manual.  I did so because x and y were clearly
    #reversed in the resulting rgb image.
    cosarr = scipy.interpolate.griddata((v.tri.bondsx, v.tri.bondsy), 
                                          cosangle, (grid_y, grid_x),method='linear')
    sinarr = scipy.interpolate.griddata((v.tri.bondsx, v.tri.bondsy), 
                                          sinangle, (grid_y, grid_x),method='linear')
    anglearr = (np.arctan2(sinarr,cosarr) + np.pi ) / (2 * np.pi)


    hsvimg = np.ones((xsize, ysize, 3))
    hsvimg[:,:,0] = anglearr
    v.rgbimg = matplotlib.colors.hsv_to_rgb(hsvimg)
#    v.rgbimg = np.flip(v.rgbimg, axis=1)
    v.rgbimg = np.flip(v.rgbimg, axis=0)
    rgbvar = v.rgbimg

    v.angleimg = v.ax.imshow(v.rgbimg, 
                              extent=[0, xsize, 0, ysize],
                              zorder=1, alpha = 0.3)
    v.imgCanvas.draw()
    
    
    



if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()
