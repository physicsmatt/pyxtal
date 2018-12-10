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
    "border": '#FF00FF', #magenta
}

def pad_locations(locs, width, size, with_indices=False):
    #This function is analogous to the numpy.pad function in "wrap" mode.
    #The input locs is an array of (x,y) coordinates.
    #This function makes additional coordinates of the points just outside
    #the range, mimicking periodic boundaries.
    #If with_indices==True, then it
    #also adds a third column, which is a reference back to the index of
    #the original vertex.
    
    #First, create the third column
    indices = np.arange(0,len(locs))
    locs = np.stack((locs[:,0], locs[:,1], indices), axis=1)

    #Now make copies of vertices within range of the edges
    #use them to make new points for padding
    x1 = locs[np.where(locs[:,0] < width)] + [size[0],0,0]
    x2 = locs[np.where(locs[:,0] > size[0] - width)] - [size[0],0,0]
    locs = np.concatenate((locs, x1, x2))
    y1 = locs[np.where(locs[:,1] < width)] + [0,size[1],0]
    y2 = locs[np.where(locs[:,1] > size[1] - width)] - [0,size[1],0]
    locs = np.concatenate((locs, y1, y2))
    if with_indices:
        return(locs)
    else:
        return(locs[:,0:2])
    

def plot_raw_image(v):
    if v.pmw.periodBound.get():
        pad = v.pmw.sphereSize[0] * 5
    else:
        pad = 0

    if v.pmw.inFileType.get() in ["image", "assemblies"]:
        xsize, ysize = v.imgshape[0], v.imgshape[1]
        v.image = np.pad(v.image, pad, "wrap")
        v.plt_image = v.ax.imshow(v.image, 
                                  extent=[-pad, xsize + pad, -pad, ysize + pad],
                                  zorder=0,
                                  cmap="gist_gray")
    elif v.pmw.inFileType.get() == "particles":
        radius = v.pmw.sphereSize[0] / 2
        locations = pad_locations(v.locations, pad, v.imgshape)
        patches = [matplotlib.patches.CirclePolygon(xy, radius,
                        facecolor=colordict["particles"]) 
                                        for xy in locations]
        background = matplotlib.patches.Rectangle((-pad,-pad),
                        width=v.imgshape[0] + 2 * pad, 
                        height=v.imgshape[1] + 2 * pad,
                        facecolor=colordict["background"])
        patches.insert(0, background)
        coll = matplotlib.collections.PatchCollection(patches, 
                                zorder=0,match_original=True)
        v.plt_image = v.ax.add_collection(coll)

    border = matplotlib.patches.Rectangle((0,0),v.imgshape[0],v.imgshape[1],
                    linewidth=2,edgecolor=colordict["border"],facecolor='none',
                    zorder = 0.5)
    v.plt_border = v.ax.add_patch(border)

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
    
    
def plot_circles(v):
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

def plot_outer_vertices(v):
#For debugging purposes, these lines plot the vertices on the perimeter:
    v.outermost = v.ax.scatter(v.locations[v.tri.outer_vertices,0], 
                               v.locations[v.tri.outer_vertices,1], 
                         color='red', zorder=5)
#    print(v.locations[v.tri.outer_vertices])


def find_outer_vertices(v):
    #Finds outermost vertices of the triangulation, in no particular order.  
    #Note that there are many other vertices besides these that are CLOSE to 
    #to being outer vertices. 

    #tri.neighbors gives the neighbors of each triangle ("simplex").
    #Index of -1 means no neighboring simplex, so the other two vertices
    #in that triangle must be on the outside.
    where_outer_tri = np.where(v.tri.neighbors == -1)[0] #rows that contain -1
    outer_tri = v.tri.simplices[where_outer_tri,:].reshape(-1)
    outer_neighbors = v.tri.neighbors[where_outer_tri,:].reshape(-1)
    return(np.unique(outer_tri[np.where(outer_neighbors != -1)]))

def calculate_triangulation(v):
    #Performs Delaunay triangulation.
    #Also finds outermost vertices, which we'll want later.
    #Also gets orientation of each bond, which we'll use later for angle field.

    locations = v.locations.copy()
    #pad_width should be 0 for non-periodic boundaries, so location not changed.
    pad_width = int(v.pmw.periodBound.get()) * v.pmw.sphereSize[0] * 5
    locations = pad_locations(locations, pad_width, 
                              v.imgshape, with_indices=True)
    v.tri = spat.Delaunay(locations[:,0:2], qhull_options="QJ")
    v.tri.outer_vertices = find_outer_vertices(v)

    #See https://docs.scipy.org/doc/scipy/reference/generated
    #                 /scipy.spatial.Delaunay.vertex_neighbor_vertices.html
    v.tri.indptr = v.tri.vertex_neighbor_vertices[0]
    v.tri.indices = v.tri.vertex_neighbor_vertices[1]

    #figure out number of bonds.  For period boundaries, num_bonds is too big,
    #and the array v.tri.bond_ind will be truncated later.
    num_tri = len(v.tri.simplices)
    num_outer = len(v.tri.outer_vertices)
    num_bonds = int((num_tri * 3 + num_outer ) / 2)  
    v.tri.bond_ind = np.zeros((num_bonds,2), dtype=int)

    #Now step through and get vertices of the ends for each bond.
    #There's probably a cute way to do the next part using np.where, but since 
    #there's only order N of these to go through, I'll just iterate instead.
    bondi = 0
    #edge_bond_count = 0  #This was used for debugging purposes
    v.tri.cnum = np.zeros(len(v.locations)) #coordination number of each vertex
    for v1i in range(0,len(v.locations)): #with these limits, v1 is NOT in padding.
        v.tri.cnum[v1i] = v.tri.indptr[v1i+1] - v.tri.indptr[v1i]
        for v2i in v.tri.indices[v.tri.indptr[v1i]:v.tri.indptr[v1i + 1]] :
            if v2i > v1i:
                v.tri.bond_ind[bondi]=[v1i,v2i]
                bondi += 1
                #if v2i >= len(v.locations):
                #    edge_bond_count +=1
    v.tri.bond_ind = v.tri.bond_ind[0:bondi]
    #Note that the bonds included here may include some bonds that extend from
    #inside the box ("real" vertices) to virtual points outside the box,
    #due to periodic boundary conditions.  As a result, some physical bonds
    #are included twice, once on each side of the image.
                
    #Now Calculate information for all of the bonds
    xy0 = v.tri.points[v.tri.bond_ind[:,0]]
    xy1 = v.tri.points[v.tri.bond_ind[:,1]]
    v.tri.segs = np.stack([xy0,xy1],axis=1)  #Coords for plotting each line segment
    v.tri.bondsxy = (xy0 + xy1) / 2 #xy center of each bond.  
    diffs = xy1 - xy0
    v.tri.bondslength = np.linalg.norm(diffs, axis = 1)
    v.tri.bondsangle = np.arctan2(diffs[:,1],diffs[:,0]) + np.pi #from 0 to 2pi
    v.tri.bondsangle %= (np.pi / 3) #Orientation, from 0 to pi/3

    #Remove the few duplicate bonds on the edges (if periodic) 
    #and calculate the median bond length
    real_bondsxy = v.tri.bondsxy % v.imgshape
    u, u_idx = np.unique(real_bondsxy, axis=0, return_index=True)
    real_bondslength = v.tri.bondslength[u_idx]
    v.median_bondlength = np.median(real_bondslength)
#I should do the previous lines differently, simply selecting the ones where the
#xy center is actually within the original image.  Also need to save xy and angle to
#the viewer structure for later calculation of orientational correlation function.

    #These lines test the code for removing unique bonds
    #dup_bonds = np.setdiff1d(np.arange(0,bondi,dtype=int), u_idx)
    #print("number of bonds removed: ",len(dup_bonds))
    #print("duplicate these bonds:")
    #print(v.tri.bondsxy[dup_bonds])
    #print("edge bond count: ", edge_bond_count)
    #v.tri.segs=v.tri.segs[dup_bonds]

    #Finally, change the indices and indptr arrays for periodic boundary 
    #conditions.  For NON-periodic boundaries, these lines should have
    #no effect.
    v.tri.indices_notwrapped = v.tri.indices.copy()
    v.tri.indices = locations[v.tri.indices,2].astype(int)
    v.tri.indptr = v.tri.indptr[0:len(v.locations)+1]


def plot_triangulation(v):
    line_coll = matplotlib.collections.LineCollection(v.tri.segs, 
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

    
def plot_disclinations(v):
    #Note that for some reason, the np.where returns a tuple with one item,
    #so I need to access it with [0].
    disc = np.where(v.tri.cnum != 6)[0]
    radius = int(v.pmw.sphereSize[0]*0.45)

    #The line below should be appreciated, as it is particularly pythonic.
    patches = [matplotlib.patches.Circle(v.locations[i], radius, 
                                         facecolor=disc_color(v.tri.cnum[i])) 
                                    for i in disc]

    coll = matplotlib.collections.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_disc = v.ax.add_collection(coll)
    v.imgCanvas.draw()


def do_label_points(v):
    for i in range(0,len(v.tri.points)):
        label = str(i)
        v.ax.annotate(label,v.tri.points[i],zorder=10,color="#FFFFFF",size=20)
    
def plot_dislocations(v):
    #Now, step through the lists as we did in do_disclinations to compile
    #list of line segments to plot:
    num_dislocs = int(np.sum(v.tri.is_dislocation))
    segs = np.zeros((num_dislocs,2,2))
    disloci = 0
    for v1 in range(0,len(v.locations)):
        for index in range(v.tri.indptr[v1], v.tri.indptr[v1 + 1]):
            v2 = v.tri.indices_notwrapped[index]
            if v.tri.is_dislocation[index] > 0 and v2 > v1:
                x1 = v.tri.points[v1,0]
                x2 = v.tri.points[v2,0]
                y1 = v.tri.points[v1,1]
                y2 = v.tri.points[v2,1]
                #print(disloci,v1,v2)
                segs[disloci] = np.array([[x1,y1],[x2,y2]])
                disloci += 1
    segs = segs[0:disloci]
    line_coll = matplotlib.collections.LineCollection(segs, 
                                  color=colordict["dislocations"], zorder=6)
    v.plt_disloc = v.ax.add_collection(line_coll)
    v.imgCanvas.draw()
    
    
def plot_unbound_discs(v):
    unbound_discs = np.where(v.tri.unboundness != 0)[0]
    radius = int(v.pmw.sphereSize[0]*0.7)

    #The line below should be appreciated, as it is particularly pythonic.
    patches = [matplotlib.patches.Circle(v.locations[i], radius, 
                         edgecolor=disc_color(v.tri.unboundness[i]+6),
                         facecolor='None') 
                                    for i in unbound_discs]

    coll = matplotlib.collections.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_unbound = v.ax.add_collection(coll)
    v.imgCanvas.draw()
  



def calculate_angle_field(v):
    #This function ultimately produces a color image representing the local
    #orientation of the crystal.

#For debugging purposes, these lines plot the centers of each bond:
#    v.bondsplot = v.ax.scatter(v.tri.bondsxy[:,0], 
#                               v.tri.bondsxy[:,1], 
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
    cosarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                          cosangle, (grid_y, grid_x),method='linear')
    sinarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                          sinangle, (grid_y, grid_x),method='linear')
    v.anglearr = (np.arctan2(sinarr,cosarr) + np.pi ) / (2 * np.pi)

    v.angle_histogram = np.histogram(v.anglearr[np.where(np.isfinite(v.anglearr))], 
                                              bins=30,range=(0,1))[0]

def plot_angle_field(v):
    xsize, ysize = v.imgshape[0], v.imgshape[1]
    hsvimg = np.ones((ysize, xsize, 3))
    hsvimg[:,:,0] = v.anglearr
    w = np.where(np.isnan(hsvimg[:,:,0]))
    hsvimg[w] = np.array([0,0,1] )
    v.rgbimg = matplotlib.colors.hsv_to_rgb(hsvimg)
    v.rgbimg = np.flip(v.rgbimg, axis=0)

    v.plt_angleimg = v.ax.imshow(v.rgbimg, 
                              extent=[0, xsize, 0, ysize],
                              zorder=1, alpha = 0.3)
    v.imgCanvas.draw()


def write_logfile_entry(v):
    #print("writing log file")
    v.pmw.logfile.write(str(v.timestep) + ', ')
    v.pmw.logfile.write(np.array2string(v.angle_histogram, 
                    max_line_width=1000, separator=',')[1:-2] + "\n")
    v.pmw.logfile.flush()
    
def do_stats(v):
     v.pmw.stats.viewer=v
     st =  v.pmw.statsText
     st.delete(1.0, "end")
     formstr = "{:28}{:10}{:10}\n"
     st.insert("end", formstr.format(" ","     Total","     Inbounds"))
     st.insert("end", formstr.format(
             "Number of spheres found:",
               len(v.locations), np.sum(v.tri.inbounds,dtype=int) ))
     st.insert("end", formstr.format(
             "    with <= 4 neighbors:",
               len(np.where(v.tri.cnum <= 4)[0]), 
               len(np.where((v.tri.cnum <= 4) * v.tri.inbounds)[0]) ))
     st.insert("end", formstr.format(
             "    with 5 neighbors:",
               len(np.where(v.tri.cnum == 5)[0]), 
               len(np.where(v.tri.cnum * v.tri.inbounds == 5)[0]) ))
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
             "    Total <> 6 neighbors:",
               len(np.where(v.tri.cnum != 6)[0]), 
               len(np.where((v.tri.cnum != 6) * v.tri.inbounds)[0]) ))
     
     st.insert("end","\n")
     st.insert("end","Number of dislocations: "  + 
                   str(int(len(np.where(v.tri.is_dislocation != 0)[0]) / 2)) + "\n")
     st.insert("end","Number of unbound disclinations (inbounds only): "  + 
                   str(len(np.where((v.tri.unboundness != 0) & v.tri.inbounds)[0])) + "\n" )
     st.insert("end","Median bond length: {:.2f}\n".format(v.median_bondlength) )


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
