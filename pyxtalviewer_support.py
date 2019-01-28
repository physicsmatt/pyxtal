#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# Support module generated by PAGE version 4.17
# In conjunction with Tcl version 8.6
#    Oct 15, 2018 09:31:29 PM CEST  platform: Linux

import numpy as np
import sys
import pyxtal_image_processing as pimg
import pyxtal_dislocations as dislocs
import pyxtal_support
import matplotlib.pyplot as plt
#import matplotlib.backends.tkagg as tkagg
#import matplotlib.backends.backend_agg 
#import FigureCanvasAgg
#https://matplotlib.org/gallery/user_interfaces/embedding_in_tk_canvas_sgskip.html
#from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

def set_Tk_var():
    pass
    #This function is a place holder, originally created by PAGE
    #It defined a zillion global Tk variables that were used for all
    #the widgets.  That functionality is now within the creation function.

    
def write_zoom_to_globals(viewer):
    #print("writing global zoom")
    viewer.pmw.global_corners_set = True
    viewer.pmw.global_corners = viewer.corners.copy()
    viewer.pmw.global_zoom = viewer.zoom
    
def read_zoom_from_globals(viewer):
    if viewer.pmw.global_corners_set:
        #print("actually reading global zoom")
        viewer.corners = viewer.pmw.global_corners.copy()
        viewer.zoom = viewer.pmw.global_zoom
        set_limits_to_corners(viewer)
        zoom_linewidths(viewer)

    
def write_views_to_globals(viewer):
    viewer.pmw.whichImage = viewer.whichImage.get()
    viewer.pmw.invertImage = viewer.invertImage.get()
    viewer.pmw.showCircles = viewer.showCircles.get()
    viewer.pmw.showTriang = viewer.showTriang.get()
    viewer.pmw.showDefects = viewer.showDefects.get()
    viewer.pmw.showOrientation = viewer.showOrientation.get()
    viewer.pmw.showStats = viewer.showStats.get()
    
def read_views_from_globals(viewer):
    viewer.whichImage.set(viewer.pmw.whichImage)
    viewer.invertImage.set(viewer.pmw.invertImage)
    viewer.showCircles.set(viewer.pmw.showCircles)
    viewer.showTriang.set(viewer.pmw.showTriang)
    viewer.showDefects.set(viewer.pmw.showDefects)
    viewer.showOrientation.set(viewer.pmw.showOrientation)
    viewer.showStats.set(viewer.pmw.showStats)

def changeVisibleImage(viewer):
    #First, decide whether to show ANYTHING:
    imagetype = viewer.whichImage.get()        
    show_image = imagetype in ["raw", "filtered"]
    viewer.plt_image.set_visible(show_image)

    #filtered and inverted options only matter for images and assemblies
    if viewer.pmw.inFileType.get() in ["image", "assemblies"]:
        invert = viewer.invertImage.get()
        cmaps = ("gist_gray","gist_yarg")
        viewer.plt_image.set_cmap(cmaps[invert])
        
        if imagetype == "raw":
            viewer.plt_image.set_data(viewer.image)
        if imagetype == "filtered":
            viewer.plt_image.set_data(viewer.filtered_image)
    
    viewer.imgCanvas.draw()

def changeVisibleAnnotations(viewer):
    viewer.plt_circles.set_visible(viewer.showCircles.get())
    viewer.plt_triang.set_visible(viewer.showTriang.get())
    viewer.plt_angleimg.set_visible(viewer.showOrientation.get())
    viewer.plt_disc.set_visible(viewer.showDefects.get())
    viewer.plt_disloc.set_visible(viewer.showDefects.get())
    viewer.plt_unbound.set_visible(viewer.showDefects.get())
    viewer.imgCanvas.draw()
#The variables below still need to be implemented and eventually included
#in the list above:
#        self.showTraject = BooleanVar()
    
#def changeStatsVisibility
    if viewer.showStats.get():
        viewer.pmw.stats.top.deiconify()
    else:
        viewer.pmw.stats.top.withdraw()
    viewer.top.focus_force()
    

def showStatsWin():
    print('pyxtalviewer_support.showStatsWin')
    sys.stdout.flush()

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

    #Add periodic wraparound, if required:
    if v.pmw.periodBound.get():
        dat = pimg.pad_locations(part_locs3d, v.pmw.sphereSize[0] * 5, boxsize3d)
    else:
        dat = part_locs3d

    labs = sklearn.cluster.dbscan(dat,eps = 2.0, min_samples = 7)[1]
    
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
    w = np.where((0 <= locations[:,0]) &
                 (locations[:,0] < v.imgshape[0]) &
                 (0 <= locations[:,1]) &
                 (locations[:,1] < v.imgshape[1]) )

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
            #pimg.do_zProfile(viewer, s[fn].particles.position[:,2], boxsize3d)
            pimg.do_zProfile(viewer, part_locs3d[:,2], boxsize3d)
        
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


def dev_to_data(xy, viewer):
    # This routine translates "device" coordinates (in pixels)
    # to "data coordinates", which are whatever is on the x and y axes.
    # Input xy can be either a tuple, list, or ndarray.
    # Output type matches input type.
    # It assumes value of pixel is same on both axes
    xpix = xy[0]
    ypix = xy[1]
    canv_xmin = viewer.canvWidget.winfo_x() 
    canv_ymin = viewer.canvWidget.winfo_y()
    canv_w = viewer.canvWidget.winfo_width()
    canv_h = viewer.canvWidget.winfo_height()

    xlims = viewer.ax.get_xlim()
    ylims = viewer.ax.get_ylim()

    ratio = (xlims[1] - xlims[0]) / canv_w
    
    xcoor = xpix * ratio + xlims[0]
    ycoor = ylims[1] - ypix * ratio
    if isinstance(xy, tuple):
        return((xcoor, ycoor))
    if isinstance(xy, np.ndarray):
        return(np.array([xcoor, ycoor]))
    if isinstance(xy, list):
        return([xcoor, ycoor])


def zoom(event, viewer):
    if str(event.type) == "ButtonPress": #Linux mouse wheel. Is there a better way?
        if event.num == 4:
            zoom_by = 1.25
        elif event.num == 5:
            zoom_by = 0.8
        else:
            print("error: zoom button not 4 or 5")
    else:
        print("error: probably a windows machine. Need to code mousewheel")
        return()
    # To zoom, first change the axes:
    viewer.zoom *= zoom_by
    mouse_xy = dev_to_data(np.array([event.x, event.y]), viewer)
    viewer.corners[0] = mouse_xy - (mouse_xy - viewer.corners[0]) / zoom_by
    viewer.corners[1] = mouse_xy + (viewer.corners[1] - mouse_xy) / zoom_by
    set_limits_to_corners(viewer)
    #Also, scale various line thicknesses accordingly
    zoom_linewidths(viewer)
    viewer.imgCanvas.draw()

    
def zoom_linewidths(v): 
    #I'll suppose for now I want the thickness to be one image pixel
    #For now, I'm just going to make a starting guess and scale it.
    #On second thought that fails for small images with big spheres.
    lw = 1 * int(np.ceil(v.zoom))

    #Instead, let's make it always be sphereSize/10.  But how do I get
    #a sphere size in points?
    lw = convert_data_to_points(v.pmw.sphereSize[0]/10, v)
    v.plt_circles.set_linewidth(lw)
    v.plt_triang.set_linewidth(lw)
    v.plt_unbound.set_linewidth(lw)
    v.plt_disloc.set_linewidth(2*lw)
    #v.imgCanvas.draw()
   
    
def convert_data_to_points(data, viewer):
    #converts data coordinates to printers' "points", used for linewidths.
    length = viewer.fig.bbox_inches.width * viewer.ax.get_position().width
    value_range = np.diff(viewer.ax.get_xlim())
    # Convert length of axis to points
    length *= 72
    # Return data_coords converted to points
    return (data * (length / value_range))



def translate(event,v):
    #translates (moves) image with mouse, when button held down.
    #Also handles double clicking.
    import datetime
    
    time_now = datetime.datetime.now()
    xy_now = np.array([event.x, event.y])
    if str(event.type) in ("ButtonPress", "Motion") and v.mousebuttondown == False:
        #Button has apparently been clicked.
        
        #First, detect if it was a double click.  If so, recenter image.
        timenow = datetime.datetime.now()
        if v.prev_button_time != None:
            ms_elapsed = (timenow - v.prev_button_time).total_seconds()
            if ms_elapsed < 0.250:
                #Yes, this was a double click.
                v.zoom = 1.0
                v.corners = np.array([[0,0],v.imgshape])
                set_limits_to_corners(v)
                zoom_linewidths(v)
                v.imgCanvas.draw()
        v.prev_button_time = timenow

        #Now set "home position" relative to which we move the image on motion.
        v.mousebuttondown = True
        v.corners_home = v.corners.copy()
        v.xy_home = xy_now.copy()
        return()
    elif str(event.type) == "ButtonRelease":
        v.mousebuttondown = False
        return()
    else: #must be a motion event with button down
        delta_data = (dev_to_data(xy_now, v) - dev_to_data(v.xy_home, v) )
        v.corners = v.corners_home - delta_data
        set_limits_to_corners(v)
        v.imgCanvas.draw()
 

def key_event(event,viewer):
    # This function parses and handles keyboard events.
    # So far, the only keyboard keys enabled are the left and right arrows,
    # used for flipping through the viewer windows.
    if event.keysym in ("Left", "Right"): #right or left arrows
        thisidx = viewer.idx
        if event.keysym == ("Left"):
            targetidx = thisidx - 1
        else:
            targetidx = thisidx + 1
            if targetidx >= len(viewer.pmw.viewers):
                targetidx = 0
        viewer.pmw.viewers[targetidx].top.lift()
        viewer.pmw.viewers[targetidx].top.focus_force()


def focus_in_event(event,viewer):
    #print("Viewer number ", viewer.idx," Gained focus.")
    pimg.display_defect_stats(viewer)
    if viewer.pmw.lockViews.get():
        read_views_from_globals(viewer)
        changeVisibleAnnotations(viewer)
        changeVisibleImage(viewer)
    if viewer.pmw.lockZoom.get():
        read_zoom_from_globals(viewer)
    changeVisibleAnnotations(viewer)
    viewer.imgCanvas.draw()
    

def focus_out_event(event,viewer):
    #print("Viewer number ", viewer.idx," Lost focus.")
    if viewer.pmw.lockViews.get():
        write_views_to_globals(viewer)
    if viewer.pmw.lockZoom.get():
        write_zoom_to_globals(viewer)


def set_limits_to_corners(viewer):
    viewer.ax.set_xlim(viewer.corners[0,0], viewer.corners[1,0])
    viewer.ax.set_ylim(viewer.corners[0,1], viewer.corners[1,1])
    #viewer.imgCanvas.draw()
    

def setup_canvas_and_axes(viewer):
    viewer.fig, viewer.ax = plt.subplots()
    viewer.ax.axis('off')
    viewer.fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)

    #Replace current placeholder canvas with new canvas object that will
    #hold all of the pyplot figures.
    viewer.top.update() #possibly required for the winfo calls to work below.
    #Get current position of existing canvs:
    canv_relx = viewer.imgCanvas.winfo_x() / viewer.top.winfo_width()
    canv_rely = viewer.imgCanvas.winfo_y() / viewer.top.winfo_height()
    canv_relw = viewer.imgCanvas.winfo_width() / viewer.top.winfo_width()
    canv_relh = viewer.imgCanvas.winfo_height() / viewer.top.winfo_height()
    #Should I delete the old canvas? I think yes, below.
    viewer.imgCanvas.destroy()
    #Here's where we actually put the plot on the tk canvas:
    viewer.imgCanvas = FigureCanvasTkAgg(viewer.fig, master=viewer.top)
    viewer.canvWidget = viewer.imgCanvas.get_tk_widget()
    viewer.canvWidget.place(relx=canv_relx, 
                            rely=canv_rely,
                            relwidth=canv_relw,
                            relheight=canv_relh)

    #Now bind the canvas to mouse and keyboard events
    viewer.canvWidget.bind("<Button-4>", lambda e:zoom(e, viewer))
    viewer.canvWidget.bind("<Button-5>", lambda e:zoom(e, viewer))
    viewer.canvWidget.bind('<B1-Motion>',lambda e:translate(e, viewer))
    viewer.canvWidget.bind('<Button-1>', lambda e:translate(e, viewer))
    viewer.canvWidget.bind('<ButtonRelease-1>', lambda e:translate(e, viewer))
    viewer.canvWidget.bind('<Double-Button-1>', lambda e:translate(e, viewer))
    viewer.top.bind("<Key>", lambda e:key_event(e, viewer))

    viewer.top.bind("<FocusIn>", lambda e:focus_in_event(e, viewer))
    viewer.top.bind("<FocusOut>", lambda e:focus_out_event(e, viewer))

    #Add some other housekeeping parts to the viewer, to keep track of
    #zooming and translation
    viewer.corners = np.array([ [0,0], viewer.imgshape ])
    set_limits_to_corners(viewer)
    viewer.prev_button_time = None
    viewer.zoom = 1.00


def init(top, viewer, *args, **kwargs):
    viewer.top = top
    viewer.pmw = args[0]
    viewer.filename = args[1]
    viewer.idx = args[2]
    viewer.framenum = args[3]
    read_views_from_globals(viewer)
    viewer.top.title("Pyxtal Viewer "
                     + "[" + str(viewer.idx) + "]: "
                     + "        " + viewer.filename + "          "
                     + "(frame: " + str(viewer.framenum) + ")"  )
    viewer.top.protocol("WM_DELETE_WINDOW", lambda: destroy_viewer(viewer))
    if viewer.pmw.batchmode.get():
        viewer.top.withdraw()

    #as this is a work in progress, I'm disabling controls that are
    #not implemented yet:
    pyxtal_support.set_widget_state('disabled', viewer.trajectCheck)
#    pyxtal_support.set_widget_state('disabled', viewer.statsCheck)

    viewer.top.update()
    viewer.mousebuttondown = False
    if viewer.pmw.inFileType.get() == "particles":
        viewer.filteredButton.configure(state='disabled')
        viewer.invertCheck.configure(state='disabled')
    load_images_and_locations(viewer)
    setup_canvas_and_axes(viewer)
    pimg.plot_raw_image(viewer)
    pimg.do_filtered_image(viewer)
#    if len(viewer.locations) == 0:
#        print("No particles found.  Killing viewer.")
#        #destroy_viewer(viewer)
    pimg.calculate_triangulation(viewer)
    pimg.calculate_angle_field(viewer)
    dislocs.calculate_dislocations(viewer)
    #When running in batchmode, and with no image output, 
    #Don't create all the plots, as they take time.
    if (not viewer.pmw.batchmode.get() or
            (viewer.pmw.outCircles.get() or viewer.pmw.outTriang.get()
                                        or viewer.pmw.outAll.get())):
        pimg.plot_circles(viewer)
        pimg.plot_triangulation(viewer)
        pimg.plot_disclinations(viewer)
        pimg.plot_dislocations(viewer)
        pimg.plot_unbound_discs(viewer)
        pimg.plot_angle_field(viewer)
        zoom_linewidths(viewer)
        changeVisibleAnnotations(viewer)
    #pimg.plot_outermost_vertices(viewer)
    #pimg.plot_outofbounds_vertices(viewer)
    pimg.calculate_defect_stats(viewer)
    pimg.display_defect_stats(viewer)
    if viewer.pmw.doOrientHist.get():
        pimg.write_orientHist_entry(viewer)
    if viewer.pmw.doDefectStats.get():
        pimg.write_defect_stats(viewer)
    pimg.do_output_files(viewer)
#    pimg.do_label_points(viewer)
    #print(dir(viewer))
    #print(dir(viewer.tri))

def destroy_viewer(viewer):
    # Function which closes the individual viewer.
    viewer.pmw.viewers.remove(viewer) #remove from main list of viewers
    plt.close(viewer.fig) #keeps the plot from reappearing in the console.
    top = viewer.top
    top.destroy()
    viewer=None

if __name__ == '__main__':
    #import pyxtalviewer
    #pyxtalviewer.vp_start_gui()
    #print("This file is not runnable as main.  Run pyxtal.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()

