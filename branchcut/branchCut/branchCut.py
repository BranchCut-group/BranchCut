#!/usr/bin/env python

"""branchCut.py: Library of functions to compute branch cuts from phase residues.
Author: Linus Ravn Gudmundsson
Co-Authors: Emil Haaber Tellefsen, Niels Schøtt Hvidberg

Date: 19/11/2024

"""

import numpy as np

def line(i1,j1,i2,j2):
    """
    Function for drawing a line between two points in a 2D array
    
    Parameters
    ----------
    i1 : int
        Row coordinate of point 1
    j1 : int
        Column coordinate of point 1
    i2 : int
        Row coordinate of point 2
    j2 : int
        Column coordinate of point 2

    Returns
    -------
    i : array_like
       1D array of row indices of pixels on line
    j : array_like
       1D array of column indices of pixels on line

    """
    if abs(i2-i1) == abs(j2-j1):
        if i1 > i2:
            i1,j1,i2,j2 = i2,j2,i1,j1
        i = np.arange(i1,i2+1)
        if j1 > j2:
            j = (j2-j1)/(i2-i1)*(i-i1)+j1+1
        else:
            j = (j2-j1)/(i2-i1)*(i-i1)+j1
                   
    elif abs(i2-i1) >= abs(j2-j1):
        if i1 > i2:
            i1,j1,i2,j2 = i2,j2,i1,j1
        i = np.arange(i1,i2+1)
        j = (j2-j1)/(i2-i1)*(i-i1)+j1
    else:
        if j1 > j2:
            i1,j1,i2,j2 = i2,j2,i1,j1
        j = np.arange(j1,j2+1)
        i = (i2-i1)/(j2-j1)*(j-j1)+i1
    
    return np.ceil(i).astype(int),np.ceil(j).astype(int)

def branch_cut(residue,mask=None,max_box_size=None):
    """
    Implementation of Goldsteins algorithm for placing branch cuts on a set of phase residues

    Parameters
    ----------
    residue : array_like
        2D array of residues
    mask : array_like
        2D boolean array of edge mask 
    max_box_size : int
        Maximum size of seach box

    Returns
    -------
    branch_cuts : array_like
        Boolean array of branch cut pixels
    """
    branch_cuts = np.zeros(residue.shape).astype(bool)
    balanced = np.zeros(residue.shape).astype(bool)
    #found_res = []

    if mask is None:
        mask = np.zeros(residue.shape).astype(bool)

    if max_box_size is None:
        max_box_size = np.min(residue.shape)

    for i in range(residue.shape[0]):
        for j in range(residue.shape[1]):

            if residue[i,j] and not balanced[i,j]:

                #found_res.append((i,j))

                # Reset active residues  and mark found residue as active
                active = np.zeros(residue.shape).astype(bool)
                active_list = []
                active[i,j] = True
                active_list.append((i,j))
                charge = residue[i,j]

                # Loop over box sizes
                for box_size in range(3,max_box_size+2,2):
                    
                    r = int((box_size-1)/2) # Radius of search (Loop over den her i stedet?)

                    # Loop throuch list of active pixels
                    N = 0
                    while N < len(active_list):

                        ia,ja = active_list[N] # Current active pixel

                        # Loop over box pixels
                        for ib in range(ia-r,ia+r+1):
                            for jb in range(ja-r,ja+r+1):

                                # Detect edges
                                if ib < 0:
                                    branch_cuts[0:ia+1,ja] = True #Top edge
                                    charge = 0

                                elif ib >= residue.shape[0]:
                                    branch_cuts[ia+1:,ja] = True #Bottom edge
                                    charge = 0

                                elif jb < 0:
                                    branch_cuts[ia,0:ja+1] = True #Left edge
                                    charge = 0

                                elif jb >= residue.shape[1]:
                                    branch_cuts[ia,ja+1:] = True #Right edge
                                    charge = 0

                                elif mask[ib,jb]:
                                    # Place branch cut to edge pixel
                                    i_bc,j_bc = line(ia,ja,ib,jb)
                                    branch_cuts[i_bc,j_bc] = True
                                    charge = 0

                                elif residue[ib,jb] and not active[ib,jb]:

                                    # Mark found residue as active
                                    active[ib,jb] = True
                                    active_list.append((ib,jb))

                                    # Mark as balanced if not already 
                                    if not balanced[ib,jb]:
                                        charge += residue[ib,jb]
                                        balanced[ib,jb] = True

                                    # Place branch cut
                                    i_bc,j_bc = line(ia,ja,ib,jb)
                                    branch_cuts[i_bc,j_bc] = True
                                if not charge:
                                    break
                            if not charge:
                                break                                          
                        N += 1
                        if not charge:
                            break

                    if not charge:
                        break
                
                # If charge still nonzero place branch cut to border (max box size exceeded)
                if charge:
                    closest_edge = np.argmin([i,residue.shape[0]-i,j,residue.shape[1]-j])
                    if closest_edge == 0:
                        branch_cuts[0:ia+1,ja] = True #Top edge
                    elif closest_edge == 1:
                        branch_cuts[ia+1:,ja] = True #Bottom edge
                    elif closest_edge == 2:
                        branch_cuts[ia,0:ja+1] = True #Left edge
                    elif closest_edge == 3:
                        branch_cuts[ia,ja+1:] = True #Right edge

                # Mark residue as balanced
                balanced[i,j] = True
      
    return branch_cuts#, found_res