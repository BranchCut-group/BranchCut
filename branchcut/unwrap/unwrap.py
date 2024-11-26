#!/usr/bin/env python

"""unwrap.py: library of functions for unwrapping SAR interferiometry images given branchcuts and phase.

Author: Emil Haaber Tellefsen
Co-Authors: Linus Ravn Gudmundsson, Niels Schøtt Hvidberg

Date: 12/11/2024
"""

# -- Third party --
import numpy as np

def unwrap(phase: np.ndarray, seed: tuple, branchCuts: np.ndarray = None, mode: str = 'bfs', unwrapBranchPixels: bool = True):
    """
    Function for unwrapping the phase of a 2D image, given a reference seed as location and branch cuts defining walls
    the unwrapping is not allowed to cross.

    Parameters
    ----------
    phase : array_like
        2-d float array with wrapped phase values
    seed : tuple
        (row,column) starting point from which phase is unwrapped from
    branchCuts : array_like, optional
        2-d boolean array where True represents the location of the branch cuts which will be avoided
    mode: string, optional
        dfs or bfs; unwrapping method for which priority order unwrapping should occur in. 
            dfs is depth-first search and unwraps most resent entry in list (stack).
            bfs is breath-first search and unwraps oldest entry in list (queue).
    unwrapBranchPixels: bool, optional
        Notes whether branchcut pixels should be unwrapped or not. Default is true.

    Returns
    -------
    array_like
        2-d float array with the unwrapped phase. Branch cuts and unwrapped zones blocked by cuts have been
        removed. Seed location is set to phase=0.
    """

    # Makes a copy of phase which will be used as unwrapped output
    unwrapped_phase = np.zeros(phase.shape)*np.nan
    
    # Initializing branchcut array as empty if none is specified
    if branchCuts is None:
        branchCuts = np.zeros(phase.shape, dtype=bool)

    # prealocating adjoin list and stack/queue structure for book keeping
    adjoin = np.zeros(phase.shape, dtype=bool)
    structure = []

    # adding seed to adjoin list and structure
    # tuple contains ((coord_r,coord_c), parentPhase, parentValue)
    adjoin[seed[0],seed[1]]=True
    structure.append((seed, phase[seed[0],seed[1]], 0))

    while len(structure)>0:
        if mode == 'dfs':
            # Depth-First Search - treating structure as stack
            cObject = structure[-1]
            structure.pop()
        
        elif mode == 'bfs':
            # Breadth-First Search - treationg structure as queue
            cObject = structure[0]
            structure.pop(0)

        # Extract tuple info
        r = cObject[0][0]
        c = cObject[0][1]
        parPhase = cObject[1]
        parVal = cObject[2]

        # unwrapping via Itoh's method (Ghiglia and Pritt (1998) p. 21 )
        phase_diff =  phase[r,c] - parPhase
        wrapped_phase_diff = np.arctan2(np.sin(phase_diff),np.cos(phase_diff))
        unwrapped_phase[r,c] = parVal + wrapped_phase_diff
        parVal = unwrapped_phase[r,c]      
        
        # Checking if neigbouring pixels are valid, does not intersect branchcut and is not already registered
        for dv in [(0,-1),(0,1),(-1,0),(1,0)]:
            r_new = r+dv[0]
            c_new = c+dv[1]
            if 0 <= r_new < phase.shape[0] and 0 <= c_new < phase.shape[1]:
                if not adjoin[r_new, c_new]:
                    if not branchCuts[r, c]: # if new and current pixel is not branch cut
                        if not branchCuts[r_new, c_new]:
                            # adding to stack/queue and registering
                            structure.append(((r_new,c_new), phase[r,c], parVal))
                            adjoin[r_new,c_new] = True
                    
                        # if new pixel is branchcut but comming from right direction and current pixel is not branchcut
                        elif unwrapBranchPixels and dv in [(0,1),(1,0)]:
                            # adding to stack/queue and registering
                            structure.append(((r_new,c_new), phase[r,c], parVal))
                            adjoin[r_new,c_new] = True    

    # Removing pixels not unpacked
    unwrapped_phase[~adjoin] = np.nan
    return unwrapped_phase



class unWrapDisplayer:
    def __init__(self, phase: np.ndarray, seed: tuple, branchCuts: np.ndarray = None, mode: str = 'dfs', unwrapBranchPixels: bool = True):
    # Makes a copy of phase which will be used as unwrapped output
        self.f_phase = np.copy(phase)
        self.seed = seed
        self.mode = mode
        self.unwrapBranchPixels = unwrapBranchPixels

        # Initializing branchcut array as empty if none is specified
        if branchCuts is None:
            self.branchCuts = np.zeros(phase.shape, dtype=bool)
        else:
            self.branchCuts = branchCuts

        # prealocating adjoin and stack/queue structure for book keeping
        self.adjoin = np.zeros(phase.shape, dtype=bool)
        self.unwrapped = np.zeros(phase.shape, dtype=bool)
        self.structure = []

        # adding seed to adjoin and structure
        # tuple contains ((coord_r,coord_c), numPeriods, parentValue)
        self.adjoin[seed[0],seed[1]]=True
        self.structure.append((seed, self.f_phase[seed[0],seed[1]], 0))
    
    def update(self):
        if len(self.structure) > 0:
            if self.mode == 'dfs':
                # Depth-First Search - treating structure as stack
                cObject = self.structure[-1]
                self.structure.pop()
            
            elif self.mode == 'bfs':
                # Breadth-First Search - treationg structure as queue
                cObject = self.structure[0]
                self.structure.pop(0)

            # Extract tuple info
            r = cObject[0][0]
            c = cObject[0][1]
            parPhase = cObject[1]
            parVal = cObject[2]
            self.unwrapped[r,c] = True

            # checking if phase difference is larger than modulus/2 and adding modulus in that case
            phase_diff =  self.f_phase[r,c] - parPhase
            wrapped_phase_diff = np.arctan2(np.sin(phase_diff),np.cos(phase_diff))
            self.f_phase[r,c] = parVal + wrapped_phase_diff
            parVal = self.f_phase[r,c]

            # Checking if neigbouring pixels are valid, does not intersect branchcut and is not already registered
            for dv in [(0,-1),(0,1),(-1,0),(1,0)]:
                r_new = r+dv[0]
                c_new = c+dv[1]
                if 0 <= r_new < self.f_phase.shape[0] and 0 <= c_new < self.f_phase.shape[1]:
                    if not self.adjoin[r_new, c_new]:
                    # if new and current pixel is not branch cut
                        if not self.branchCuts[r_new, c_new] and not self.branchCuts[r, c]:
                            # adding to stack/queue and registering
                            self.structure.append(((r_new,c_new), self.f_phase[r,c], parVal))
                            self.adjoin[r_new,c_new] = True
                        # if new pixel is branchcut but comming from right direction and current pixel is not branchcut
                        elif self.unwrapBranchPixels and dv in [(0,1),(1,0)] and not self.branchCuts[r, c]:
                            # adding to stack/queue and registering
                            self.structure.append(((r_new,c_new), self.f_phase[r,c], parVal))
                            self.adjoin[r_new,c_new] = True    
            return False
        else:
            return True 