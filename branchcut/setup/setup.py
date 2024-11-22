# This module contains the different functions for the BranchCut algorithm
import numpy as _np


def makephase(N, M=None):
    """
    This Function creates a linear increasing phase diagram of size (N x M), or (N x N) if M is not given
    """
    if M == None:
        M = N
    
    data = _np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            data[i,j] = j*_np.pi*0.4+i*_np.pi*0.4

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
    return (data+_np.pi)%(2*_np.pi) - _np.pi


def add_residues(data, location, percent=False):
    data_out = _np.copy(data)

    if percent:
        loc = [int(loc[0] * size[0]), int(loc[1] * size[1])]

    for row in location:
        if row[2] == 1:
            data_out[row[0], row[1]] += 1.2*_np.pi
    
        else:
            data_out[row[0], row[1]] -= 1.2*_np.pi
    
    return data_out


def create_mask(size: list or tuple, loc: list = [0.5,0.5], shape: int = 'line', percent: bool = True):
    """
    Create a mask in the shape of a line, square, triangle, or circle.
    """
    if percent:
        loc = [int(loc[0] * size[0]), int(loc[1] * size[1])]
    
    mask = _np.zeros(size, dtype=_np.uint8)
    
    match shape:
        case 'line':
            pos = int(size[1]/2)
            hl  = int(loc[1]/2)
            mask[pos, loc[0]-hl:loc[0]+hl] = 1

        case 'square':
            half_size = min(size) // 7
            mask[loc[0]-half_size:loc[0]+half_size, loc[1]-half_size:loc[1]+half_size] = 1

        case 'triangle':
            pos = int(size[1]/2)
            width = loc[1] + (loc[1]%2 - 1)
            if width >= size[1]:
                width = size[1] - (size[1]%2 + 1)
            
            height = (width + 1) // 2
            if loc[0] + _np.ceil(height/2) >= size[0]:
                height = loc[0] - size[0]

            for i in range(height):
                mask[loc[0]-height//2 + i, pos-i:pos+i+1] = 1

        case 'circle':
            rr, cc = _np.ogrid[:size[0], :size[1]]
            circle = (rr - loc[0]) ** 2 + (cc - loc[1]) ** 2 <= (min(size) // 6) ** 2
            mask[circle] = 1

        case _:
            raise ValueError("Invalid shape. Choose from 'line', 'square', 'triangle', or 'circle'.")

    return mask

def find_residues(data):
    """
    Find residues in a wrapped phase diagram
    """
    # Can be parallelized
    up      = data[:-1,1:]-data[:-1,:-1]
    right   = data[1:,1:]-data[:-1,1:]
    down    = data[1:,:-1]-data[1:,1:]
    left    = data[:-1,:-1]-data[1:,:-1]
    out = _np.zeros(data.shape,dtype=float)
    out[:-1,:-1] = (wrap(up) + wrap(right) + wrap(down) + wrap(left))/(2*_np.pi)
    return _np.round(out).astype(int)

if __name__ == "__main__":
    print("This is a Goldstein Branch Cut algorithm module.")