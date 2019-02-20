#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:35:41 2018

@author: mtrawick
"""

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import matplotlib
import matplotlib.collections as m_coll

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
#    "merge": '#00FFFF', #cyan
#    "split": '#FF00FF', #magenta
    "merge": np.array((1.0, 0.0, 1.0)), #magenta
    "split": np.array((0.0, 1.0, 1.0)), #cyan
    "spontaneous": np.array((1.0, 0.66, 0.0)), #orange
}

def pad_locations(locs, width, size, with_indices=False):
    #This function is analogous to the numpy.pad function in "wrap" mode.
    #The input locs is an array of (x,y) or (x,y,z) coordinates.
    #This function makes additional coordinates of the points just outside
    #the range (only in x and y), mimicking periodic boundaries.
    #If with_indices==True, then the function also appends
    #an additional column, which is a reference back to the index of
    #the original vertex.
    
    #First, create the additional column
    indices = np.arange(0,len(locs))
    locs = np.concatenate((locs, indices[np.newaxis].T), axis=1)

    #Now make copies of vertices within range of the edges
    #use them to make new points for padding
    offset = np.zeros(locs.shape[1])
    offset[0] = size[0]
    x1 = locs[np.where(locs[:,0] < width)] + offset
    x2 = locs[np.where(locs[:,0] > size[0] - width)] - offset
    locs = np.concatenate((locs, x1, x2))

    offset = np.zeros(locs.shape[1])
    offset[1] = size[1]
    y1 = locs[np.where(locs[:,1] < width)] + offset
    y2 = locs[np.where(locs[:,1] > size[1] - width)] - offset
    locs = np.concatenate((locs, y1, y2))
    if with_indices:
        return(locs)
    else:
        return(locs[:,0:-1])
    

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
        coll = m_coll.PatchCollection(patches, 
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
    
    try:  #v.sphere_properties may be undefined if the code was run using 2d sphere finding
        ellipse_array = np.hstack((v.locations, 
                                   v.ellipse_axes[:,[0,2]],
                                   v.ellipse_axis_rot.reshape(1,-1).T) )
        ellipse_tuple = tuple(map(tuple,ellipse_array))
        patches = [matplotlib.patches.Ellipse((x,y), maj_ax, min_ax, rot) 
           for x,y, min_ax, maj_ax, rot in ellipse_tuple]
        #print("used ellipses")
    except:
        #print("using circles")
        radius = int(v.pmw.sphereSize[0]*0.7)
        patches = [matplotlib.patches.CirclePolygon(xy, radius) 
                        for xy in v.locations]
        
    coll = m_coll.PatchCollection(patches, 
                            edgecolor=(colordict["circles"]), facecolor='None', zorder=2)
    v.plt_circles = v.ax.add_collection(coll)
    v.imgCanvas.draw()

def plot_outermost_vertices(v):
#For debugging purposes, these lines plot the vertices on the perimeter:
    v.outermost = v.ax.scatter(v.locations[v.tri.outer_vertices,0], 
                               v.locations[v.tri.outer_vertices,1], 
                         color='red', zorder=5)
#    print(v.locations[v.tri.outer_vertices])

def plot_outofbounds_vertices(v):
#For debugging purposes, these lines plot the vertices that are not inbounds:
    outofbounds = np.where(v.tri.inbounds == False)[0]
    v.outofbounds = v.ax.scatter(v.locations[outofbounds,0], 
                               v.locations[outofbounds,1], 
                         color='black', zorder=4.9)
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
    import scipy.spatial
    #Performs Delaunay triangulation.
    #Also finds outermost vertices, which we'll want later.
    #Also gets orientation of each bond, which we'll use later for angle field.

    locations = v.locations.copy()
    #pad_width should be 0 for non-periodic boundaries, so location not changed.
    pad_width = int(v.pmw.periodBound.get()) * v.pmw.sphereSize[0] * 5
    locations = pad_locations(locations, pad_width, 
                              v.imgshape, with_indices=True)
    v.tri = scipy.spatial.Delaunay(locations[:,0:2], qhull_options="QJ")
    v.tri.outer_vertices = find_outer_vertices(v)

    #See https://docs.scipy.org/doc/scipy/reference/generated
    #                 /scipy.spatial.Delaunay.vertex_neighbor_vertices.html
    v.tri.indptr = v.tri.vertex_neighbor_vertices[0]
    v.tri.indices = v.tri.vertex_neighbor_vertices[1]

    v.tri.cnum = v.tri.indptr[1:len(v.locations)+1] - v.tri.indptr[0:len(v.locations)]
    right_indices = v.tri.indices[0:v.tri.indptr[len(v.locations)]]
    is_first_neighbor = np.zeros(len(right_indices), dtype=int)
    is_first_neighbor[v.tri.indptr[0:len(v.locations)]] = 1
    v.tri.left_indices = np.cumsum(is_first_neighbor) - 1
    w = np.where(right_indices > v.tri.left_indices)
    bond_ind = np.vstack((v.tri.left_indices[w], right_indices[w])).T.copy()
    #Note that the bonds included here may include some bonds that extend from
    #inside the box ("real" vertices) to virtual points outside the box,
    #due to periodic boundary conditions.  As a result, some physical bonds
    #are included twice, once on each side of the image.
                
    #Now Calculate information for all of the bonds
    xy0 = v.tri.points[bond_ind[:,0]]
    xy1 = v.tri.points[bond_ind[:,1]]
    v.tri.segs = np.stack([xy0,xy1],axis=1)  #Coords for plotting each line segment
    v.tri.bondsxy = (xy0 + xy1) / 2 #xy center of each bond.  
    diffs = xy1 - xy0
    v.tri.bondslength = np.linalg.norm(diffs, axis = 1)
    v.tri.bondsangle = np.arctan2(diffs[:,1],diffs[:,0]) + np.pi #from 0 to 2pi
    v.tri.bondsangle %= (np.pi / 3) #Orientation, from 0 to pi/3

    #Remove the few duplicate bonds on the edges (if periodic) 
    #and calculate the median bond length
    #There are two possible ways to remove duplicate bonds.  Below, we simply
    #find all of the fonds that are unique
    real_bondsxy_with_dupes = v.tri.bondsxy % v.imgshape
    u, w = np.unique(real_bondsxy_with_dupes, axis=0, return_index=True)
    #Below is an alternate method, which selects only bonds which have their
    #centers within the real image.  However, it fails in the (rare) case of
    #two vertices just outside the corner of an image, but the bond is within
    #the image.  This leads to the wrong number of bonds.
    #w = np.where((0 <= v.tri.bondsxy[:,0]) &
    #             (v.tri.bondsxy[:,0] < v.imgshape[0]) &
    #             (0 <= v.tri.bondsxy[:,1]) &
    #             (v.tri.bondsxy[:,1] < v.imgshape[1]))
    #The real bondsxy and bondsangle will be used for the correlation function.
    v.tri.real_bondsxy = v.tri.bondsxy[w]
    v.tri.real_bondsangle = v.tri.bondsangle[w]
    v.median_bondlength = np.median(v.tri.bondslength[w])

    #These lines test the code for removing unique bonds
    #dup_bonds = np.setdiff1d(np.arange(0,bondi,dtype=int), w)
    #print(v.imgshape)
    #print("number of bonds removed: ",len(dup_bonds))
    #print("these are the bonds that were duplicates:")
    #print(v.tri.bondsxy[dup_bonds])
    #print("edge bond count: ", edge_bond_count)
    #print("number of real bonds is: ",len(v.tri.real_bondsxy))
    #print("number of real vertices is: ",len(v.locations))
    #v.tri.segs=v.tri.segs[dup_bonds]

    #Finally, change the indices and indptr arrays for periodic boundary 
    #conditions.  For NON-periodic boundaries, these lines should have
    #no effect.
    v.tri.indices_notwrapped = v.tri.indices.copy()
    v.tri.indices = locations[v.tri.indices,2].astype(int)
    v.tri.indptr = v.tri.indptr[0:len(v.locations)+1]


def plot_triangulation(v):
    line_coll = m_coll.LineCollection(v.tri.segs, 
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

    coll = m_coll.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_disc = v.ax.add_collection(coll)
    v.imgCanvas.draw()


def do_label_points(v):
    for i in range(0,len(v.tri.points)):
        label = str(i)
        v.ax.annotate(label,v.tri.points[i],zorder=10,color="#FFFFFF",size=20)
    
def plot_dislocations(v):
    right_indices = v.tri.indices_notwrapped[0:v.tri.indptr[len(v.locations)]]
    w = np.where((right_indices > v.tri.left_indices) &
                 (v.tri.is_dislocation[0:len(right_indices)] > 0))
    disloc_ind = np.vstack((v.tri.left_indices[w], right_indices[w])).T.copy()

    xy0 = v.tri.points[disloc_ind[:,0]]
    xy1 = v.tri.points[disloc_ind[:,1]]
    segs = np.stack([xy0,xy1],axis=1)  #Coords for plotting each line segment

    line_coll = m_coll.LineCollection(segs, 
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

    coll = m_coll.PatchCollection(patches, match_original=True, zorder=5)
    v.plt_unbound = v.ax.add_collection(coll)
    v.imgCanvas.draw()
  

def calculate_angle_field(v):
    import scipy.interpolate

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
    try:
        #The issue here is that when doing linear interpolation, the griddata routine calls
        #qhull, which throws an exception in the case of all input points being colinear.
        #This problem would be easily fixable using the "QJ" (joggle) option, but I don't
        #see a way to include that here.
        cosarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                              cosangle, (grid_y, grid_x),method='linear')
        sinarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                              sinangle, (grid_y, grid_x),method='linear')
    except:
        cosarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                              cosangle, (grid_y, grid_x),method='nearest')
        sinarr = scipy.interpolate.griddata((v.tri.bondsxy[:,0], v.tri.bondsxy[:,1]), 
                                              sinangle, (grid_y, grid_x),method='nearest')
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
                              zorder=1, alpha = 0.5)
    v.imgCanvas.draw()


def calculate_defect_stats(v):
    #This function builds a 2d array that holds the counts of various defects.
    #[:,0] holds the number that are "inbounds" (not on the edge).
    #[:,1] holds the total number.
    defects = np.zeros((9,2),dtype=int)

    #total number of spheres
    defects[0,0] = np.sum(v.tri.inbounds,dtype=int)
    defects[0,1] = len(v.locations)

    # <=4 nearest neighbors
    defects[1,0] = len(np.where((v.tri.cnum <= 4) * v.tri.inbounds)[0])
    defects[1,1] = len(np.where(v.tri.cnum <= 4)[0])

    # =5 nearest neighbors
    defects[2,0] = len(np.where((v.tri.cnum == 5) * v.tri.inbounds)[0])
    defects[2,1] = len(np.where(v.tri.cnum == 5)[0])

    # =6 nearest neighbors
    defects[3,0] = len(np.where((v.tri.cnum == 6) * v.tri.inbounds)[0])
    defects[3,1] = len(np.where(v.tri.cnum == 6)[0])

    # =7 nearest neighbors
    defects[4,0] = len(np.where((v.tri.cnum == 7) * v.tri.inbounds)[0])
    defects[4,1] = len(np.where(v.tri.cnum == 7)[0])

    # >=8 nearest neighbors
    defects[5,0] = len(np.where((v.tri.cnum >= 8) * v.tri.inbounds)[0])
    defects[5,1] = len(np.where(v.tri.cnum >= 8)[0])

    # total <> 6 nearest neighbors
    defects[6,0] = len(np.where((v.tri.cnum != 6) * v.tri.inbounds)[0])
    defects[6,1] = len(np.where(v.tri.cnum != 6)[0])

    # unbound disclinations
    defects[7,0] = len(np.where((v.tri.unboundness != 0) & v.tri.inbounds)[0])
    defects[7,1] = len(np.where(v.tri.unboundness != 0)[0])

    # dislocations (ALL of them have at least one vertex inbounds)
    defects[8,0] = int(len(np.where(v.tri.is_dislocation != 0)[0]) / 2)
    defects[8,1] = defects[8,0]


    # which dislocations have at least ONE end inbounds:
#    right_index_inbounds = v.tri.inbounds[v.tri.indices]
#    is_first_neighbor = np.zeros(len(v.tri.indices), dtype=int)
#    is_first_neighbor[v.tri.indptr[0:-2]] = 1
#    left_index = np.cumsum(is_first_neighbor)
#    left_index_inbounds = v.tri.inbounds[left_index]
#    inbounds_dislocation = ((v.tri.is_dislocation != 0) & 
#                            (left_index_inbounds | right_index_inbounds))
    
    v.defect_stats = defects
    
def display_defect_stats(v):
    v.pmw.stats.viewer=v
    st =  v.pmw.statsText
    st.delete(1.0, "end")

    defects = v.defect_stats
    formstr = "{:32}{:8}{:10}\n"
    st.insert("end", formstr.format(" ","  Inbounds","   Total"))
    st.insert("end", formstr.format(
            "Number of spheres found:", defects[0,0], defects[0,1]))
    st.insert("end", formstr.format(
            "    with <= 4 neighbors:", defects[1,0], defects[1,1]))
    st.insert("end", formstr.format(
            "    with 5 neighbors:", defects[2,0], defects[2,1]))
    st.insert("end", formstr.format(
            "    with 6 neighbors:", defects[3,0], defects[3,1]))
    st.insert("end", formstr.format(
            "    with 7 neighbors:", defects[4,0], defects[4,1]))
    st.insert("end", formstr.format(
            "    with >= 8 neighbors:", defects[5,0], defects[5,1]))
    st.insert("end", formstr.format(
            "    Total <> 6 neighbors:", defects[6,0], defects[6,1]))
    st.insert("end","\n")

    st.insert("end", formstr.format(
            "Number of unbound disclinations:", defects[7,0], defects[7,1]))
    st.insert("end", formstr.format(
            "Number of dislocations:", defects[8,0], defects[8,1]))
    st.insert("end","\n")

    st.insert("end","Median bond length: {:.2f}\n".format(v.median_bondlength) )


def plot_trajectories(v):
    dummy_segments = ()
    dummy_trajectories = m_coll.LineCollection(dummy_segments)
    v.plt_trajectories = v.ax.add_collection(dummy_trajectories)


def find_merge_events(pmw, tracks):
    #This function returns arrays that mark the beginning and ending
    #of each particle's trajectory.  The beginnings and ends are coded
    #according to whether it appears that a new particle was created as
    #part of a MERGE with another particle, or two partles were created
    #as part of a SPLIT of one particle into two.
    #Incidentally, this is a pretty nasty, low-level routine.  I tried to
    #do it with dataframes instead of numpy arrays, but I just couldn't
    #figure out how to get past the problems of "chain assignment".

    tr = tracks.values #convert dataframe to ndarray.
        #column 0: x
        #column 1: y
        #column 2: frame
        #column 3: particle number
        #column 4: (added later) color code for merge/split
        
    num_parts = int(np.max(tr[:,3])) + 1
    num_frames = int(np.max(tr[:,2])) + 1
    
    #compile lists of beginnings and endings of each trajectory.
    beglist = []
    endlist = []
    for p in range(0,num_parts):
        wp = np.where(tr[:,3]==p)
        thistrack = tr[wp]
        beginning = thistrack[0].copy()
        ending = thistrack[-1].copy()
        ending[2] +=1
        beglist.append(beginning)
        endlist.append(ending)

    #Convert lists to arrays, remove entries of trajectories that start on 
    #first frame or end of last frame, and add additional "column" for
    #code to denote merge/split/spontaneous.
    begarr = np.array(beglist)
    endarr = np.array(endlist)
    begarr = begarr[np.where(begarr[:,2] != 0)]
    endarr = endarr[np.where(endarr[:,2] != num_frames)]
    begarr = np.hstack((begarr, np.zeros((len(begarr),1))))
    endarr = np.hstack((endarr, np.zeros((len(endarr),1))))
    
    cutoff_dist = pmw.sphereSize[0] * 1.5 
    #min distance for events to be considered part of merge or split.
    #should really be related to sphere spacing, not size, but spacing
    #is technically different for every frame.  Note that this doesn't
    #account well for merges between a sphere and an already very long
    #ellipse.  If I really wanted to do it well, I'd use data on the
    #apparent aspect ratio and angle of each microdomain.

    #for each frame, look for new trajectories near two ending ones,
    #or old trajectories near two beginning ones.
    for f in range(0,num_frames):
        wfb = np.where(begarr[:,2]==f)
        begs = begarr[wfb]
        wfe = np.where(endarr[:,2]==f)
        ends = endarr[wfe]
        for b in range (0,len(begs)):
            dist = np.linalg.norm(ends[:,0:2] - begs[b,0:2], axis=1)
            w = np.where(dist <= cutoff_dist)
            if len(w[0]) == 2:
                endarr[wfe[0][w],4] = 1 #code for merge
                begarr[wfb[0][b],4] = 1 #code for merge
        for e in range (0,len(ends)):
            dist = np.linalg.norm(begs[:,0:2] - ends[e,0:2], axis=1)
            w = np.where(dist <= cutoff_dist)
            if len(w[0]) == 2:
                begarr[wfb[0][w],4] = 2 #code for split
                endarr[wfe[0][e],4] = 2 #code for split

    return(begarr, endarr)
        

def redo_trajectories(pmw):
    import copy
    import trackpy as tp
    
    search_range = pmw.sphereSize[0] / 2
    tracks = tp.link_df(pmw.feature_df, search_range)

    seglist = [np.array(tracks[["x","y"]].loc[tracks["particle"]==p])
               for p in range(tracks["particle"].max() + 1)]
    #I don't understand the line about traj_colors at all. Copied from example.
    traj_colors = [matplotlib.colors.to_rgba(c)
          for c in matplotlib.rcParams['axes.prop_cycle'].by_key()['color']]
    trajectories = m_coll.LineCollection(seglist, colors=traj_colors)

    beginnings, endings = find_merge_events(pmw, tracks)
    beg_colors = np.zeros((len(beginnings),3),dtype = float)
    beg_colors[np.where(beginnings[:,4]==1)] = colordict["merge"]
    beg_colors[np.where(beginnings[:,4]==2)] = colordict["split"]
    beg_colors[np.where(beginnings[:,4]==0)] = colordict["spontaneous"]

    end_colors = np.zeros((len(endings),3),dtype = float)
    end_colors[np.where(endings[:,4]==1)] = colordict["merge"]
    end_colors[np.where(endings[:,4]==2)] = colordict["split"]
    end_colors[np.where(endings[:,4]==0)] = colordict["spontaneous"]

    for v in pmw.viewers:
        cp_traj = copy.copy(trajectories)
        v.plt_trajectories.remove()
        v.plt_trajectories = v.ax.add_collection(cp_traj)
        
        w = np.where(beginnings[:,2]==v.idx)
        v.ax.scatter(beginnings[w,0], beginnings[w,1], marker="P", s=150,
                  c = beg_colors[w])

        w = np.where(endings[:,2]==v.idx)
        v.ax.scatter(endings[w,0], endings[w,1], marker="$*$", s=150,
                  c = end_colors[w])

        v.imgCanvas.draw()
    

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
        v.plt_circles.set_visible(False)
        v.plt_triang.set_visible(True)
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


def write_orientHist_entry(v):
    #print("writing log file")
    v.pmw.orientHistfile.write(str(v.timestep) + ', ')
    v.pmw.orientHistfile.write(np.array2string(v.angle_histogram, 
                    max_line_width=1000, separator=',')[1:-1] + "\n")
    v.pmw.orientHistfile.flush()

    
def do_zProfile(v, z_coords, box):
    bin_width = 0.1
    distribution, bin_edges = np.histogram(z_coords, 
                                bins = int(np.ceil(box[2] / bin_width)),
                                range = (-box[2]/2, box[2]/2))

    volume = bin_width * box[0] * box[1]
    distribution = distribution.astype(float) / volume
#    v.pmw.zProfilefile.write(str(v.timestep) + ', ')
#    v.pmw.zProfilefile.write(np.array2string(distribution, 
#                    max_line_width=1000, separator=',',
#                    formatter={'float_kind':lambda x: " %.3f" % x})[1:-1] + "\n")
    zdf = pd.DataFrame({"timestep": v.timestep,
                        "z-depth": bin_edges[0:-1],
                        "density": distribution})
    zdf.to_csv(v.pmw.zProfilefile, sep='\t', float_format='%.3f',
               index_label="#index")
    v.pmw.zProfilefile.flush()


def do_Sphere_Stats(v, m, ellipse_axes):
    #This function writes distribution information about both the mass and the
    #aspect ratio of the spheres
    
    #v.pmw.sphereStatsfile.write(str(v.timestep) + ', ')

    #first, the aspect ratio
    aspects = ellipse_axes[:,2] / ellipse_axes[:,0]
    aspect_distrib, aspect_bins = np.histogram(aspects, 
                                    bins = (0.5, 1.5, 2.5, 3.5, 4.5, 1000))
#    v.pmw.sphereStatsfile.write(np.array2string(aspect_distrib, 
#                    max_line_width=1000, separator=', ')[1:-1]+ ", ")
    asp_df = pd.DataFrame({"timestep": v.timestep,
                        "aspect_ratio": aspect_bins[0:-1],
                        "count": aspect_distrib})
    asp_df.to_csv(v.pmw.sphereStatsfile, sep='\t', float_format='%.3f',
                  index_label="#index")


    #now the masses
    bin_width = 10
    m_max = np.max(m)
    num_bins = int(np.ceil(m_max / bin_width))
    mass_distrib, mass_bins = np.histogram(m, bins = num_bins, 
                                range = (0, num_bins * bin_width))
#    v.pmw.sphereStatsfile.write(np.array2string(distribution, 
#                    max_line_width=1000, separator=', ')[1:-1] + "\n")
    asp_df = pd.DataFrame({"timestep": v.timestep,
                        "mass": mass_bins[0:-1],
                        "count": mass_distrib})
    asp_df.to_csv(v.pmw.sphereStatsfile, sep='\t', float_format='%.3f')

    v.pmw.sphereStatsfile.flush()
    

def write_defect_stats(v):
    output_array = v.defect_stats.flatten()
    v.pmw.defectStatsFile.write(str(v.timestep) + ', ')
    v.pmw.defectStatsFile.write(np.array2string(output_array, 
                    max_line_width=1000, separator=',')[1:-1])
    v.pmw.defectStatsFile.write(', %.2f\n' % v.median_bondlength)
    v.pmw.defectStatsFile.flush()

    
    
if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()
