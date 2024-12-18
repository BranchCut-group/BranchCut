# This module contains the different functions for the BranchCut algorithm
import numpy as np

def makephase(N, M=None):
    """
    This Function creates a linear increasing phase diagram of size (N x M), or (N x N) if M is not given
    """
    if M == None:
        M = N
    
    data = np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            data[i,j] = j*np.pi*0.4+i*np.pi*0.4

    return data

def wrap(data):
    """
    This function can wrap a phase matrix.
    == input ==
    data:       ndarray with phases

    == output == 
    data_out:   ndarray of the wrapped or unwrapped values. Same
                size as data.
    """

    # Wrap the phase
    return (data+np.pi)%(2*np.pi) - np.pi

def addresidues(data, location, percentage=False):
    """
    This Function add residues to a phase diagram. If line is set to True add a horizontal line with center in the given locations.
    == Input ==
    data:       (N x M) ndarray of phase diagram.

    location:   (N x 2) ndarray where N is number of residues. Contains (x, y) location of residues or (x, L) location and length 
                of residue lines.

    line:       Boolean (default False) determines wether the residue is a line or point. If the residue is a line second row of 
                location is used as length.

    percentage: Boolean (default False) determines wehter all values should be interpreted as percentage of the input dimensions.

    == Output ==
    data_out:   (N x M) ndarray of phase diagram with residues.
    """
    data_out = np.copy(data)

    for row in location:
        if row[2] == 1:
            data_out[row[0], row[1]] += np.pi

        else:
            data_out[row[0], row[1]] -= np.pi
    
    return data_out

def find_residues(data):
    # Can be parallelized
    up      = data[:-1,1:]-data[:-1,:-1]
    right   = data[1:,1:]-data[:-1,1:]
    down    = data[1:,:-1]-data[1:,1:]
    left    = data[:-1,:-1]-data[1:,:-1]

    residue_map = up + right + down + left

    print(residue_map)

if __name__ == "__main__":
    print("This is a Goldstein Branch Cut algorithm module.")