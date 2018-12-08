#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 12:41:02 2018

@author: mtrawick
"""


import numpy as np

def create_inbounds_list(v):
    #creates a list of NEAR edge vertices, which probably have weird
    #numbers of neighbors.
    if v.pmw.periodBound.get():
        v.tri.inbounds = np.ones(len(v.locations), dtype=bool)
        return()
    
    bvertsidx = v.tri.outer_vertices
    bverts = v.tri.points[bvertsidx]
    centerx = (np.max(bverts[:,0]) + np.min(bverts[:,0])) / 2
    centery = (np.max(bverts[:,1]) + np.min(bverts[:,1])) / 2
    disp = bverts - np.array([centerx,centery])
    theta = np.arctan2(disp[:,1],disp[:,0])
    #sort these.
    order = np.argsort(theta)
    bverts = bverts[order]
    
    #There's probably a cuter way to do this:
    numpoints = len(v.tri.points)
    #numverts = len(bverts)
    v.tri.inbounds = np.zeros(numpoints, dtype=bool)
    threshhold = v.median_bondlength * 0.7
    p1 = bverts
    p2 = np.roll(bverts, 1, axis = 0)
    denom = 1/ np.linalg.norm(p2-p1, axis=1)
    for i in range(0,numpoints):
        p3 = v.tri.points[i]
        ds = np.abs(np.cross(p2-p1,p3-p1, axis = 1) * denom)
        min_d = np.min(ds)
        v.tri.inbounds[i] = (min_d > threshhold)

    #For debugging purposes, these lines plot the vertices on the perimeter:
#    w = np.where(v.tri.inbounds)
#    v.plt_inbounds = v.ax.scatter(v.locations[w,0], 
#                               v.locations[w,1], 
#                         color='red', zorder=5)

def neighbor_inds(p,v):
    return(np.array(range(v.tri.indptr[p], v.tri.indptr[p+1])))

    
def neighbors(p,v):
    return(v.tri.indices[neighbor_inds(p,v)])


def make_bond_between(p1, p2, v):
    #print("making bond between ",p1,p2)
    sign_p1 = np.sign(v.tri.unboundness[p1])
    sign_p2 = np.sign(v.tri.unboundness[p2])
    v.tri.unboundness[p1] += sign_p2
    v.tri.unboundness[p2] += sign_p1
    
    wp2 = np.where(neighbors(p1,v) == p2)
    v.tri.is_dislocation[neighbor_inds(p1,v)[wp2]] +=1

    wp1 = np.where(neighbors(p2,v) == p1)
    v.tri.is_dislocation[neighbor_inds(p2,v)[wp1]] +=1


def break_bond_between(p1, p2, v):
    v.tri.unboundness[p1] += np.sign(v.tri.cnum[p1] - 6)
    v.tri.unboundness[p2] += np.sign(v.tri.cnum[p2] - 6)

    wp2 = np.where(neighbors(p1,v) == p2)
    v.tri.is_dislocation[neighbor_inds(p1,v)[wp2]] -=1

    wp1 = np.where(neighbors(p2,v) == p1)
    v.tri.is_dislocation[neighbor_inds(p2,v)[wp1]] -=1


def can_retract_from(p1, p2, v, edges_ok):
    #print("Recursion level: ", v.tri.recursion_level)
    if v.tri.recursion_level > 500:
        return(False)
    v.tri.inprocess[p1]=True
    break_bond_between(p1, p2, v)

    if not v.tri.inbounds[p2]:
        v.tri.inprocess[p1]=False
        return(True)
        
    v.tri.recursion_level += 1
    
    if find_a_mate_nicely(p2, v, edges_ok):
        v.tri.inprocess[p1]=False
        v.tri.recursion_level -= 1
        return(True)

    if find_a_mate_rudely(p2, v, edges_ok):
        v.tri.inprocess[p1]=False
        v.tri.recursion_level -= 1
        return(True)
        
    v.tri.recursion_level -= 1
    v.tri.inprocess[p1]=False
    make_bond_between(p1, p2, v)
    return(False)
    

def can_butt_in_on(p1, p2, v, edges_ok):
    if (not v.tri.inbounds[p2]) and (not edges_ok): return(False)
    if not v.tri.buttInOnAble[p2, int(edges_ok)]: return(False)
    v.tri.inprocess[p1]=True
    for i in neighbor_inds(p2,v):
        p3 = v.tri.indices[i]
        if v.tri.is_dislocation[i] and (v.tri.inbounds[p3] or edges_ok):
            if can_retract_from(p2, p3, v, edges_ok):
                v.tri.inprocess[p1]=False
                return(True)
    v.tri.inprocess[p1]=False
    v.tri.buttInOnAble[p2,int(edges_ok)] = False
    return(False)
    #Note: I really wonder about setting p1 in process.  Isn't it already in process?
    #And don't I screw stuff up if I make it NOT inprocess at the end?
                

def helpful_to_butt_in_on(p1, p2, v, edges_ok):
    helpful = ((np.sign(v.tri.cnum[p1] - 6) == -np.sign(v.tri.cnum[p2] - 6))
                and not v.tri.inprocess[p2]
                and v.tri.unboundness[p1] != 0)  #added to algorith, 2018, probably redundant.
    #Also, what if a 4 is already bonded to a 7.  Is it still helpful to butt in on itself?
    #if helpful: print("helpful:", p1, p2)
    return(helpful)
               

def find_a_mate_rudely(p1, v, edges_ok):
   # print("finding mate RUDELY for", p1)
    for p2 in neighbors(p1,v) :
        if helpful_to_butt_in_on(p1, p2, v, edges_ok):
            if can_butt_in_on(p1, p2, v, edges_ok):
                make_bond_between(p1, p2, v)
                return(True)
    return(False)


def helpful_to_bond_nicely(p1,p2,v,edges_ok):
    helpful = ((np.sign(v.tri.unboundness[p1]) == -np.sign(v.tri.unboundness[p2]))
                    and not v.tri.inprocess[p2]
                    and (v.tri.inbounds[p2] or edges_ok))
    return(helpful)


def find_a_mate_nicely(p1, v, edges_ok):
    #print("finding mate nicely for",p1)
    for p2 in neighbors(p1,v) :
        if helpful_to_bond_nicely(p1, p2, v, edges_ok):
            make_bond_between(p1, p2, v)
            return(True)
    return(False)


def find_mates_for(p,v):
    if v.tri.inbounds[p]:
        v.tri.inprocess[p] = True
        #I suspect that these inprocess statements are unnecessary in this funciton.
        while v.tri.unboundness[p] != 0:
            if not find_a_mate_nicely(p, v, edges_ok=False): break

        while v.tri.unboundness[p] != 0:
            if not find_a_mate_rudely(p, v, edges_ok=False): break

        while v.tri.unboundness[p] != 0:
            if not find_a_mate_nicely(p, v, edges_ok=True): break

        while v.tri.unboundness[p] != 0:
            if not find_a_mate_rudely(p, v, edges_ok=True): break

        v.tri.inprocess[p] = False


def calculate_dislocations(v):
    create_inbounds_list(v)
    v.tri.is_dislocation = np.zeros(len(v.tri.indices), dtype=np.int16)
    #the bond between a 4 and an 8 would be a "double bond", with a 2 here.
    v.tri.inprocess=np.full(len(v.locations), False)
    v.tri.buttInOnAble=np.full((len(v.locations), 2), True)
    v.tri.unboundness = v.tri.cnum.copy() - 6
    v.tri.recursion_level = 0
    for i in range(0, len(v.locations)):
        find_mates_for(i, v)

    
if __name__ == '__main__':
    #print("This file is not runnable as main.  Run Pyxtalmain.py instead.")
    import pyxtal
    pyxtal.vp_start_gui()
