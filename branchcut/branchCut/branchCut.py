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

def box_search(A,i,j,box_size):
    """
    Function for searching an n x n box around a point in a 2D array

    Parameters
    ----------
    A : 2D array
        Array to be searched
    i,j : int
        Coordinates of center of search box
    box size : int
        Size of search box 
        
    Returns
    -------
    edge : int
        Integer indicating an edge, 0 = top edge, 1 = left edge, 2 = bottom edge, 3 = right edge

    or

    inds : array_like
        2D array of coordinates of nonzero elements of A within search box

    """
    r = int((box_size-1)/2)

    # Edge detection
    edges = np.nonzero([i-r < 0, j-r < 0, i+r >= A.shape[0], j+r >= A.shape[1]])[0]
    if edges.size:
            return edges[0]
    
    # If no edge was found, return indices of all nonzero points in search box 
    else:
        inds = np.array(np.nonzero(A[i-r:i+r+1,j-r:j+r+1]))
        inds += np.array([[i-r],[j-r]])
        return inds
        
def branch_cut(residue,mask=None,max_box_size=None):
    """
    Implementation of Goldstein algorithm for placing branch cuts on a set of phase residues

    Parameters
    ----------
    residue : array_like
        2D array of residues
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

                    # Loop throuch list of active pixels
                    N = 0
                    while N < len(active_list):

                        m,n = active_list[N] # Current active pixel

                        # Search box around active pixel for non-active residues
                        res_inds = box_search(residue.astype(bool) & np.logical_not(active),m,n,box_size)

                        if type(res_inds) == np.ndarray: # No edge found
                            
                            # Loop over all found residues in search box
                            for k in range(res_inds.shape[1]):

                                p = tuple(res_inds[:,k])

                                # Mark found residue as balanced if not already
                                if not balanced[p]:
                                    charge += residue[p]
                                    balanced[p] = True

                                # Mark residue as active and add to list of active pixels
                                active[p] = True
                                active_list.append(p)

                                # Place branch cut
                                i_bc,j_bc = line(m,n,p[0],p[1])
                                branch_cuts[i_bc,j_bc] = True
                                #branch_cuts[i_bc,j_bc] += 1
                                #found_res.append(p)

                                if not charge:
                                    break
                        else: # Edge found
                            if res_inds == 0:
                                branch_cuts[0:m+1,n] = True #Top edge
                            elif res_inds == 1:
                                branch_cuts[m,0:n+1] = True #Left edge
                            elif res_inds == 2:
                                branch_cuts[m+1:,n] = True #Bottom edge
                            elif res_inds == 3:
                                branch_cuts[m,n+1:] = True #Right edge
                            charge = 0
                        
                        N += 1
                        if not charge:
                            break

                    if not charge:
                        break

                balanced[i,j] = True

                # TODO If charge still nonzero place branch cut to border (max box size exceeded)
    return branch_cuts#, found_res