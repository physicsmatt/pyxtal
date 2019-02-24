#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 08:58:49 2019

@author: mtrawick
"""
import numpy as np
import pyxtal_image_processing as pimg
import pandas as pd


def get_contiguous_regions_from_image(img):
    #This function is similar to trackpy in that it returns locatios of
    #features from an image.  But the "features" here are specifically
    #contiguous areas of non-zero density.
    #This routine is designed to be called in the same way as trackpy.locate,
    #in the sense that it doesn't read or write data directly from the
    #pyxtal viewer object.
    #Below, there is a similar and generally better routine called
    #get_locations_from_3d, which finds contiguously connected regions in
    #three dimensions.  It is to be preferred, in part because it also
    #records masses and major/minor axes of elliptical particles.
    #(I started to get data like mass from the spheres here, but I ended up
    #not doing anything with it, because doing so would require rewriting
    #this routine as to make calling it very different from calling
    #trackpy.locate.  I think there's currently no need to make this routine
    #any fancier; Think of it as an alternative to trackpy.locate.
    import skimage.measure as skim
    
    nonzero = (img > 0)
    labels = skim.label(nonzero, connectivity=2)
    props = skim.regionprops(labels, intensity_image = img, cache=False)
    locations = np.array([p.weighted_centroid for p in props])

    #mass data below is currently not used or stored except as a cutoff.  
    #See note above.
    masses = np.array([p.weighted_moments[0,0] for p in props])

    #Set the mass cutoff
    w = np.where(masses > 6)
    #print(np.sort(masses))
    return(locations[w])
    

def get_locations_from_image(v, method="trackpy"):
    #The main work in this function is done by trackpy.locate.
    #This function is a wrapper that fixes some quirks in the
    #coordinate system and deals with the case of periodic
    #boundary conditions
    import trackpy

    img_for_locs = v.image
    
    #Add periodic wraparound, if required:
    if v.pmw.periodBound.get():
        wrap_width = v.pmw.sphereSize[0] * 5
    else:
        wrap_width = 0
    img_for_locs = np.pad(img_for_locs, wrap_width, "wrap")
    
    # Apparently the locate function reverses the y coordinate, so
    # the [::-1] notation above verses the array top-to-bottom.
    img_for_locs = img_for_locs[::-1]

    if method == "contiguous":
        locations = get_contiguous_regions_from_image(img_for_locs)
    else: #presumably using "trackpy"
        full_locations = trackpy.locate(img_for_locs, v.pmw.sphereSize[0])
        #Trackpy gives dataframe with 8 columns. First two columns are y, x 
        locations = np.array(full_locations)[:,0:2]

    # Trackpy also puts y before x, so flip them back:
    locations = np.flip(locations, axis = 1)

    #Finally, remove the effects of the wraparound, and kill any out of bounds.
    locations -= np.array([wrap_width, wrap_width])
    w = np.where((0 <= locations[:,0]) &
                 (locations[:,0] < v.imgshape[0]) &
                 (0 <= locations[:,1]) &
                 (locations[:,1] < v.imgshape[1]) )

    return(locations[w])


def inertia_tensor(vert):
    #What this actually calculates is the moment of inertia tensor for a set
    #of points about its COM, normalized by the number of points.
    #(This seems to be what skimage.measure.regionprops does, and I'm
    #aiming for this output to behave similarly.)
    com = np.average(vert,axis=0)
    vertcom = vert - com
    i00 = np.sum(vertcom[:,1]**2 + vertcom[:,2]**2)
    i11 = np.sum(vertcom[:,0]**2 + vertcom[:,2]**2)
    i22 = np.sum(vertcom[:,0]**2 + vertcom[:,1]**2)
    i01 = -np.sum(vertcom[:,0]*vertcom[:,1])
    i12 = -np.sum(vertcom[:,1]*vertcom[:,2])
    i02 = -np.sum(vertcom[:,0]*vertcom[:,2])
    it = np.array([[i00,i01,i02],[i01,i11,i12],[i02,i12,i22]])
    #normalize by mass:
    it /= len(vert)
    return(it)


def get_locations_from_dbscan(v, part_locs3d, boxsize3d):
    import sklearn
    import sklearn.cluster

    #Add periodic wraparound, if required:
    dat = pimg.pad_locations(part_locs3d, v)

    labs = sklearn.cluster.dbscan(dat,eps = 2.0, min_samples = 7)[1]
    
    if np.max(labs) == -1: #No clusters found!
        v.masses = np.array([])
        v.ellipse_axes = np.array([]).reshape(0,3)
        v.ellipse_axis_rot = np.array([])
        return(np.array([]).reshape(0,2))
        
    clusters = np.arange(np.max(labs) + 1)
    locations = np.array([np.average(dat[np.where(labs==c)], axis=0) 
                                for c in clusters])

    #kill any out of bounds.
    w = np.where((0 <= locations[:,0]) &
                 (locations[:,0] < v.imgshape[0]) &
                 (0 <= locations[:,1]) &
                 (locations[:,1] < v.imgshape[1]) )
    locations = locations[w]
    clusters = clusters[w]

    masses = np.array([len(labs[np.where(labs==c)]) for c in clusters])
    it = np.array([inertia_tensor(dat[np.where(labs==c)]) for c in clusters])

    #find eigenvalues and vectors of inertia tensor
    iw, iv = np.linalg.eigh(it)
    v.ellipse_axes = 4 * np.sqrt(iw)
    
    #y and x components of zeroeth eigenvector.
    #The first eigenvector is the smallest moment of inertia, from eigh()
    v.ellipse_axis_rot = np.degrees(np.arctan2(iv[:,1,0],iv[:,0,0]))
    
    if v.pmw.doSphereStats.get():
        pimg.do_Sphere_Stats(v, masses, v.ellipse_axes)

    return(locations[:,0:2])


def get_locations_from_3d(v, part_locs3d, boxsize3d):
    import skimage.measure as skim

    image3d = np.histogramdd(part_locs3d, bins=boxsize3d,
                    range = np.stack((np.array([0.,0.,0.]),boxsize3d)).T )[0]

    #Add periodic wraparound, if required:
    if v.pmw.periodBound.get():
        wrap_width = v.pmw.sphereSize[0] * 5
    else:
        wrap_width = 0
    pad_spec = ((wrap_width,wrap_width),(wrap_width,wrap_width),(0,0))
    image3d = np.pad(image3d, pad_spec, "wrap")

    nonzero = (image3d > 0)
    labels = skim.label(nonzero, connectivity=3)
    props = skim.regionprops(labels, intensity_image = image3d, cache=False)
    locations = np.array([p.weighted_centroid for p in props])

    #remove the effects of the wraparound, and kill any out of bounds.
    locations -= np.array([wrap_width, wrap_width, 0])
    locations = locations[:,0:2]
    w = np.where(pimg.in_orig_box(locations, v))
    locations = locations[w]

    #use props[w] to get masses, moments
    masses = np.array([p.weighted_moments[0,0,0] for p in props])[w]

    #Set a mass cutoff
    w_mc = np.where(masses > 6)
    v.sphere_masses = masses[w_mc]
    locations = locations[w_mc]

    #find eigenvalues and vectors of inertia tensor
    it = np.array([p.inertia_tensor for p in props])[w][w_mc]
    iw, iv = np.linalg.eigh(it)
    v.ellipse_axes = 4 * np.sqrt(iw)
    
    #y and x components of zeroeth eigenvector.
    #The first eigenvector is the smallest moment of inertia, from eigh()
    v.ellipse_axis_rot = np.degrees(np.arctan2(iv[:,1,0],iv[:,0,0]))
    
    if v.pmw.doSphereStats.get():
        pimg.do_Sphere_Stats(v, masses, v.ellipse_axes)

    return(locations)


def pad_with_bullshit_points(v):
    #The Delauney triangulation requires at least 4 points to function, 
    #and other parts of the code presumably also have a minimum.  This
    #function ensures that there are always at least 4 points, creating
    #up to four completely bullshit points out of thin air as necessary.
    
    if len(v.locations) >=4: return()

    print("Warning: Number of points found is less than 4.")
    print("Filling in with extra completely bogus points.")
    x1 = v.imgshape[0] * 0.1
    x2 = v.imgshape[0] * 0.9
    y1 = v.imgshape[1] * 0.1
    y2 = v.imgshape[1] * 0.9
    bs_locations = np.array([[x1,y1],[x2,y1],[x1,y2],[x2,y2]])
#    bs_masses = [1.0, 1, 1, 1]
    bs_ellipse_axes = [[1.0,1, 1], [1,1,1], [1,1,1], [1,1,1]]
    bs_ellipse_axis_rot = [0,0,0,0]

    last_bs_index = 4 - len(v.locations)
    v.locations = np.append(v.locations, bs_locations[0:last_bs_index], axis=0)
#    np.append(v.sphere_masses, bs_masses[0:last_bs_index])
    v.ellipse_axes = np.append(v.ellipse_axes, bs_ellipse_axes[0:last_bs_index], axis=0)
    v.ellipse_axis_rot = np.append(v.ellipse_axis_rot, bs_ellipse_axis_rot[0:last_bs_index])
    

def load_images_and_locations(viewer):
    #Based on the input file type, this function reads the file.
    #If File is an image, it adds the location data.
    #If File is location data, it adds a fake "image" of spheres.
    #If File is assemblies, it calcultes both an image and location data.

    import gsd.hoomd
    import os

    full_filename = os.path.join(viewer.pmw.path, viewer.filename)
    if viewer.pmw.inFileType.get() == "image":
        #use code from colloid group.
        viewer.image = plt.imread(full_filename)
#below I can clip the image to something smaller, just for debugging purposes
#        viewer.image = viewer.image[0:600,0:900]
        viewer.imgshape = np.flip(np.shape(viewer.image))  #Note that order now [x, y]
        viewer.locations = get_locations_from_image(viewer)        

        
    elif viewer.pmw.inFileType.get() == "particles":
        #read gsd file
        s = gsd.hoomd.open(name=full_filename, mode='rb')
        fn = viewer.framenum  #this value was passed in from pmw
        if fn > len(s):
            print("ERROR: frame number out of range")
        viewer.locations = s[fn].particles.position[:,0:2].copy() #do I need a copy here?

        #The hoomd box seems to be always centered on 0,0;
        #I shift the particle locations so 0,0 is the lower left corner.
        boxsize = s[fn].configuration.box[0:2]  #z-component not used. assuming x,y.
        viewer.locations += boxsize / 2
        viewer.timestep = s[fn].configuration.step

        #I also scale all locations by an arbitrary scale factor.
        #Each new unit corresponds to a pixel in the displayed orientation image.
        scale_fact = 10
        viewer.locations *= scale_fact
        boxsize = np.ceil(boxsize) * scale_fact

        #set the global sphere size from diameter of first particle
        viewer.pmw.sphereSize[0] = s[fn].particles.diameter[0] * scale_fact
        viewer.imgshape = np.array([int(boxsize[0]),int(boxsize[1])])

    else: #must be a gsd assembly
        #read gsd file
        s = gsd.hoomd.open(name=full_filename, mode='rb')
        fn = viewer.framenum  #this value was passed in from pmw
        if fn >= len(s):
            print("ERROR: frame number out of range")
        boxsize3d = s[fn].configuration.box[0:3]  #z-component not used. assuming x,y.
        boxsize3d = np.ceil(boxsize3d)
        viewer.imgshape = np.array([int(boxsize3d[0]),int(boxsize3d[1])])
        viewer.timestep = s[fn].configuration.step
        
        #Get locations of only 'B' particles, or whatever type is specified.
        typestring = viewer.pmw.partTypeStr.get()
        typenum = s[fn].particles.types.index(typestring)
        #Note that for some reason, the np.where returns a tuple with one item,
        #so I need to access it with [0].
        w = np.where(s[fn].particles.typeid == typenum)[0]
        part_locs3d = s[fn].particles.position[w,0:3]
        
        if viewer.pmw.doZProfile.get():
            #for all particles:
            pimg.do_zProfile(viewer, s[fn].particles.position[:,2], boxsize3d)
            #for particles of only selected type ("B" or whatever):
            #pimg.do_zProfile(viewer, part_locs3d[:,2], boxsize3d)
        
        #now create an "image" based on the densities of particles at x,y locations
        part_locs3d += boxsize3d / 2
        part_locs = part_locs3d[:,0:2]
        image = np.histogram2d(part_locs[:,1], part_locs[:,0],
                               bins = np.flip(viewer.imgshape),
                               range = ((0,boxsize3d[1]),(0,boxsize3d[0])))[0]
        image = np.flip(image,axis=0)
        viewer.image = image.copy()
        
        #Below are four possible methods to extract the location of spheres:
        #viewer.locations = get_locations_from_image(viewer, method="trackpy")        
        #viewer.locations = get_locations_from_image(viewer, method="contiguous")        
        #viewer.locations = get_locations_from_3d(viewer, part_locs3d, boxsize3d)
        viewer.locations = get_locations_from_dbscan(viewer, part_locs3d, boxsize3d)

        pad_with_bullshit_points(viewer)

    wrapped_locations = pimg.pad_locations(viewer.locations, viewer)
    feature_df = pd.DataFrame({"x": wrapped_locations[:,0],
                               "y": wrapped_locations[:,1],
                               "frame": viewer.idx})
    
    viewer.pmw.feature_df = viewer.pmw.feature_df.append(feature_df, ignore_index=True)


if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()

