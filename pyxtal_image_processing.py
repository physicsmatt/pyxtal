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
import pyxtal_dislocations as dislocs


colordict =	{
    "particles": "0.7",
    "background": "0.3",
    "triangulation": "blue",
    "circles": "#B4FF64", #RGB color (180,255,100)
    "disc4": '#FF00FF', #magenta
    "disc5": '#FF0000', #red
    "disc7": '#00FF00', #green
    "disc8": '#00FFFF', #cyan
    "dislocations": '#FFFF00', #yellow
}

def do_raw_image(v):
    if v.pmw.inFileType.get() in ["image", "assemblies"]:
        xsize, ysize = v.imgshape[0], v.imgshape[1]
        v.plt_image = v.ax.imshow(v.image, 
                                  extent=[0, xsize, 0, ysize],
                                  zorder=0,
                                  cmap="gist_gray")
    elif v.pmw.inFileType.get() == "particles":
        radius = v.pmw.sphereSize[0] / 2
        patches = [matplotlib.patches.CirclePolygon(xy, radius,
                        facecolor=colordict["particles"]) 
                                        for xy in v.locations]
        background = matplotlib.patches.Rectangle((0,0),
                        width=v.imgshape[0],height=v.imgshape[1],
                        facecolor=colordict["background"])
        patches.insert(0, background)
        coll = matplotlib.collections.PatchCollection(patches, 
                                zorder=0,match_original=True)
        v.plt_image = v.ax.add_collection(coll)
    v.imgCanvas.draw()


def do_filtered_image(v):
    import scipy.ndimage
    if v.pmw.inFileType.get() in ["image", "assemblies"]:
        sigma = v.pmw.sphereSize[0]/3
        input_ = np.fft.fft2(v.image)
        output_ = scipy.ndimage.fourier_gaussian(input_, sigma)
        v.filtered_image = np.fft.ifft2(output_).real
        #now rescale the output to match input
        v.filtered_image *= np.max(v.image)/np.max(v.filtered_image)
    
    
def do_circle_plot(v):
    #Plot circles around the locations of the spheres.
    #I've used a collection of circle "patches", so that the radius can
    #scale with the axes as the figure is zoomed.
    #I'll still have to rescale the linewidth manually, however.
    

    radius = int(v.pmw.sphereSize[0]*0.7)
    patches = [matplotlib.patches.CirclePolygon(xy, radius) 
                    for xy in v.locations]
    coll = matplotlib.collections.PatchCollection(patches, 
                            edgecolor=(colordict["circles"]), facecolor='None', zorder=2)
    v.plt_circles = v.ax.add_collection(coll)

    v.imgCanvas.draw()


def do_triangulation(v):
    #Performs Delaunay triangulation and creates a plot of it.
    #Also finds outermost vertices, which we'll want later.
    #Also gets orientation of each bond, which we'll use later for angle field.

    v.tri = spat.Delaunay(v.locations, qhull_options="Qj")
    #v.plt_triang = v.ax.triplot(v.locations[:,0], v.locations[:,1],
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
    #above are the indices of the outer vertices.  Note that there are many
    #others that are CLOSE to the outside besides these.  These are in no 
    #particular order.

#For debugging purposes, these lines plot the vertices on the perimeter:
#    v.outermost = v.ax.scatter(v.locations[v.tri.outer_vertices,0], 
#                               v.locations[v.tri.outer_vertices,1], 
#                         color='red', zorder=5)
#    print(v.locations[v.tri.outer_vertices])
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
    v.tri.indptr = v.tri.vertex_neighbor_vertices[0]
    v.tri.indices = v.tri.vertex_neighbor_vertices[1]
    #There's probably a cute way to do the next part using np.where, but since 
    #there's only order N of these to go through, I'll just iterate instead.
    bondi = 0
    for v1i in range(0,len(v.tri.points)):
        v.tri.cnum[v1i] = v.tri.indptr[v1i+1] - v.tri.indptr[v1i]
        for v2i in v.tri.indices[v.tri.indptr[v1i]:v.tri.indptr[v1i + 1]] :
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
    line_coll = matplotlib.collections.LineCollection(segs, 
                                  color=colordict["triangulation"], zorder=4)
    v.plt_triang = v.ax.add_collection(line_coll)
    v.imgCanvas.draw()

def disc_color(cnum):
    # These are the rgb colors used for various types of disclinations.
    # I should figure out how to implement this with a python dictionary,
    # but I don't know how to handle the cases <=4 or >=8.
    if cnum <= 4:
        return(colordict["disc4"])
    if cnum == 5:
        return(colordict["disc5"])
    if cnum == 7:
        return(colordict["disc7"])
    if cnum >=8:
        return(colordict["disc8"])

    
def do_disclinations(v):
    #Note that for some reason, the np.where returns a tuple with one item,
    #so I need to access it with [0].
    v.tri.disc = np.where(v.tri.cnum != 6)[0]
    radius = int(v.pmw.sphereSize[0]*0.45)

    #The line below should be appreciated, as it is particularly pythonic.
    patches = [matplotlib.patches.Circle(v.locations[i], radius, 
                                         facecolor=disc_color(v.tri.cnum[i])) 
                                    for i in v.tri.disc]

    coll = matplotlib.collections.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_disc = v.ax.add_collection(coll)
    v.imgCanvas.draw()


def do_label_points(v):
    for i in range(0,len(v.tri.points)):
        label = str(i)
        v.ax.annotate(label,v.tri.points[i],zorder=10,color="#FFFFFF",size=20)
    
def do_dislocations(v):
    dislocs.find_dislocations(v)
    #Now, step through the lists as we did in do_disclinations to compile
    #list of line segments to plot:
    num_dislocs = int(np.sum(v.tri.is_dislocation) / 2)
    segs = np.zeros((num_dislocs,2,2))
    disloci = 0
    for v1 in range(0,len(v.tri.points)):
        for index in range(v.tri.indptr[v1], v.tri.indptr[v1 + 1]):
            v2 = v.tri.indices[index]
            if v.tri.is_dislocation[index] > 0 and v2 > v1:
                x1 = v.tri.points[v1,0]
                x2 = v.tri.points[v2,0]
                y1 = v.tri.points[v1,1]
                y2 = v.tri.points[v2,1]
                #print(disloci,v1,v2)
                segs[disloci] = np.array([[x1,y1],[x2,y2]])
                disloci += 1
    line_coll = matplotlib.collections.LineCollection(segs, 
                                  color=colordict["dislocations"], zorder=6)
    v.plt_disloc = v.ax.add_collection(line_coll)
    v.imgCanvas.draw()
    
def do_unbound_discs(v):
    v.tri.unbound_discs = np.where(v.tri.unboundness != 0)[0]
    radius = int(v.pmw.sphereSize[0]*0.7)

    #The line below should be appreciated, as it is particularly pythonic.
    patches = [matplotlib.patches.Circle(v.locations[i], radius, 
                         edgecolor=disc_color(v.tri.unboundness[i]+6),
                         facecolor='None') 
                                    for i in v.tri.unbound_discs]

    coll = matplotlib.collections.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_unbound = v.ax.add_collection(coll)
    v.imgCanvas.draw()
  



def do_angle_field(v):
    #This function ultimately produces a color image representing the local
    #orientation of the crystal.

#For debugging purposes, these lines plot the centers of each bond:
#    v.bondsplot = v.ax.scatter(v.tri.bondsx, 
#                               v.tri.bondsy, 
#                         color='red', zorder=15)

    #Now, create the orientation field:
    xsize, ysize = v.imgshape[0], v.imgshape[1]
#    grid_x, grid_y = np.mgrid[0:xsize:1, 0:ysize:1]
    grid_x, grid_y = np.mgrid[0:ysize:1, 0:xsize:1]
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

    hsvimg = np.ones((ysize, xsize, 3))
    hsvimg[:,:,0] = anglearr
    w = np.where(np.isnan(hsvimg[:,:,0]))
    hsvimg[w] = np.array([0,0,1] )
    v.rgbimg = matplotlib.colors.hsv_to_rgb(hsvimg)
    v.rgbimg = np.flip(v.rgbimg, axis=0)

    v.plt_angleimg = v.ax.imshow(v.rgbimg, 
                              extent=[0, xsize, 0, ysize],
                              zorder=1, alpha = 0.3)
    v.imgCanvas.draw()
    
    
def do_stats(v):
     v.pmw.stats.viewer=v
     st =  v.pmw.statsText
     st.delete(1.0, "end")
     formstr = "{:28}{:10}{:10}\n"
     st.insert("end", formstr.format(" ","     Total","     Inbounds"))
     st.insert("end", formstr.format(
             "Number of spheres found:",
               len(v.tri.points), np.sum(v.tri.inbounds,dtype=int) ))
     st.insert("end", formstr.format(
             "    with <= 4 neighbors:",
               len(np.where(v.tri.cnum <= 4)[0]), 
               len(np.where((v.tri.cnum <= 4) * v.tri.inbounds)[0]) ))
     st.insert("end", formstr.format(
             "    with 5 neighbors:",
               len(np.where(v.tri.cnum == 5)[0]), 
               len(np.where(v.tri.cnum * v.tri.inbounds == 4)[0]) ))
     st.insert("end", formstr.format(
             "    with 6 neighbors:",
               len(np.where(v.tri.cnum == 6)[0]), 
               len(np.where(v.tri.cnum * v.tri.inbounds == 6)[0]) ))
     st.insert("end", formstr.format(
             "    with 7 neighbors:",
               len(np.where(v.tri.cnum == 7)[0]), 
               len(np.where(v.tri.cnum * v.tri.inbounds == 7)[0]) ))
     st.insert("end", formstr.format(
             "    with >= 8 neighbors:",
               len(np.where(v.tri.cnum >= 8)[0]), 
               len(np.where((v.tri.cnum >= 8) * v.tri.inbounds)[0]) ))
     st.insert("end", formstr.format(
             "    Total != 6 neighbors:",
               len(np.where(v.tri.cnum != 6)[0]), 
               len(np.where((v.tri.cnum != 6) * v.tri.inbounds)[0]) ))
     
     st.insert("end","\n")
     st.insert("end","Number of dislocations: "  + 
                   str(int(len(np.where(v.tri.is_dislocation != 0)[0]) / 2)) + "\n")
     st.insert("end","Number of unbound disclinations (inbounds only): "  + 
                   str(len(np.where((v.tri.unboundness != 0) & v.tri.inbounds)[0])) )


def do_output_files(v):
    import os
    splitpoint = v.filename.rfind('.')
    base = os.path.join(v.pmw.path, v.filename[0:splitpoint])
    if v.pmw.inFileType.get() != "image":
        base += "_f" + str(v.framenum)
    if v.pmw.outCircles.get():
        v.plt_circles.set_visible(True)
        v.plt_triang.set_visible(False)
        v.plt_angleimg.set_visible(False)
        v.plt_disc.set_visible(False)
        v.plt_disloc.set_visible(False)
        v.plt_unbound.set_visible(False)
        v.fig.tight_layout(pad=-1.08)
        v.fig.savefig(base + "_circ.tif")#, bbox_inches='tight',pad_inches=0)
    if v.pmw.outAll.get():
        v.plt_circles.set_visible(True)
        v.plt_triang.set_visible(False)
        v.plt_angleimg.set_visible(True)
        v.plt_disc.set_visible(True)
        v.plt_disloc.set_visible(True)
        v.plt_unbound.set_visible(True)
        
        v.fig.savefig(base + "_all.tif", bbox_inches='tight',pad_inches=0)
        X = np.array(v.fig.canvas.renderer._renderer)
    if v.pmw.outTriang.get():
        v.plt_circles.set_visible(False)
        v.plt_triang.set_visible(True)
        v.plt_angleimg.set_visible(False)
        v.plt_disc.set_visible(False)
        v.plt_disloc.set_visible(False)
        v.plt_unbound.set_visible(False)
        v.fig.savefig(base + "_triang.tif", bbox_inches='tight',pad_inches=0)

    #print(X.shape)
    #Use X as a frame for an mpeg output.
    
#    fig2 = plt.figure()
#    ax2 = fig2.add_subplot(111, frameon=False)
#    ax2.imshow(X)
#    plt.show()
#    plt.imsave("testout.tif", X)


if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()
