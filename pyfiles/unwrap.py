#!/usr/bin/env python

"""unwrap.py: library of functions for unwrapping SAR interferiometry images given branchcuts and phase.

Author: Emil Haaber Tellefsen
Co-Authors: Linus Ravn Gudmundsson, Niels Schøtt Hvidberg

Date: 09/11/2024
"""

# -- Third party --
import numpy as np

def unwrap(phase: np.ndarray, seed: tuple, branchCuts: np.ndarray = None, mode: str = 'dfs', modulus: float = 2*np.pi):
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
    modulus: float, optional
        modulus of wrapped phase. Default is 2*pi

    Returns
    -------
    array_like
        2-d float array with the unwrapped phase. Branch cuts and unwrapped zones blocked by cuts have been
        removed. Seed location is set to phase=0.
    """

    # Makes a copy of phase which will be used as unwrapped output
    f_phase = np.copy(phase)
    
    # Initializing branchcut array as empty if none is specified
    if branchCuts is None:
        branchCuts = np.zeros(phase.shape, dtype=bool)

    # prealocating registry and stack/queue structure for book keeping
    registry = np.zeros(phase.shape, dtype=bool)
    structure = []

    # adding seed to registry and structure
    # tuple contains ((coord_r,coord_c), numPeriods, parentValue)
    registry[seed[0],seed[1]]=True
    structure.append((seed, 0, f_phase[seed[0],seed[1]]))

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
        numPeriods = cObject[1]
        parVal = cObject[2]

        # adding already unwrapped periods
        f_phase[r,c] += numPeriods*modulus

        # checking if phase difference is larger than modulus/2 and adding modulus in that case
        phase_diff = parVal - f_phase[r,c]
        if abs(phase_diff) > modulus/2:
            f_phase[r,c] +=np.sign(phase_diff)*modulus
            numPeriods+=np.sign(phase_diff)
        
        # Checking if neigbouring pixels are valid, does not intersect branchcut and is not already registered
        for dv in [(0,-1),(0,1),(-1,0),(1,0)]:
            r_new = r+dv[0]
            c_new = c+dv[1]
            if 0 <= r_new < f_phase.shape[0] and 0 <= c_new < f_phase.shape[1]:
                if not branchCuts[r_new, c_new]:
                    if not registry[r_new, c_new]:
                        # adding to stack/queue and registering
                        structure.append(((r_new,c_new), numPeriods, f_phase[r,c]))
                        registry[r_new,c_new] = True
    
    # Removing branchcuts and pixels not unpacked
    f_phase[branchCuts] = np.nan   
    f_phase[~registry] = np.nan
    f_phase -= f_phase[seed[0],seed[1]]
    return f_phase



class unWrapDisplayer:
    def __init__(self, phase: np.ndarray, seed: tuple, branchCuts: np.ndarray = None, mode: str = 'dfs', modulus: float = 2*np.pi):
    # Makes a copy of phase which will be used as unwrapped output
        self.f_phase = np.copy(phase)
        self.seed = seed
        self.mode = mode
        self.modulus = modulus

        # Initializing branchcut array as empty if none is specified
        if branchCuts is None:
            self.branchCuts = np.zeros(phase.shape, dtype=bool)
        else:
            self.branchCuts = branchCuts

        # prealocating registry and stack/queue structure for book keeping
        self.registry = np.zeros(phase.shape, dtype=bool)
        self.unwrapped = np.zeros(phase.shape, dtype=bool)
        self.structure = []

        # adding seed to registry and structure
        # tuple contains ((coord_r,coord_c), numPeriods, parentValue)
        self.registry[seed[0],seed[1]]=True
        self.structure.append((seed, 0, self.f_phase[seed[0],seed[1]]))
    
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
            numPeriods = cObject[1]
            parVal = cObject[2]
            self.unwrapped[r,c] = True

            # adding already unwrapped periods
            self.f_phase[r,c] += numPeriods*self.modulus

            # checking if phase difference is larger than modulus/2 and adding modulus in that case
            phase_diff = parVal - self.f_phase[r,c]
            if abs(phase_diff) > self.modulus/2:
                self.f_phase[r,c] +=np.sign(phase_diff)*self.modulus
                numPeriods+=np.sign(phase_diff)
            
            # Checking if neigbouring pixels are valid, does not intersect branchcut and is not already registered
            for dv in [(0,-1),(0,1),(-1,0),(1,0)]:
                r_new = r+dv[0]
                c_new = c+dv[1]
                if 0 <= r_new < self.f_phase.shape[0] and 0 <= c_new < self.f_phase.shape[1]:
                    if not self.branchCuts[r_new, c_new]:
                        if not self.registry[r_new, c_new]:
                            # adding to stack/queue and registering
                            self.structure.append(((r_new,c_new), numPeriods, self.f_phase[r,c]))
                            self.registry[r_new,c_new] = True   
            return False
        else:
            return True 
    