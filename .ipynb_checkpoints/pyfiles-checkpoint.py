# This module contains the different functions for the BranchCut algorithm
import numpy as np

def wrap(data, unwrap=False):
    """
    This function can wrap a phase matrix.
    input:
        data - ndarray with phases

    output: 
        data_out - ndarray of the wrapped or unwrapped values. Same
                   size as data.
    """
    # Wrap the phase
    data_out = np.copy(data)
    data_out = (data_out+np.pi)%(2*np.pi) - np.pi

    return data_out


if __name__ == "__main__":
    print("This is a Goldstein Branch Cut algorithm module.")