#!/usr/bin/env python

"""exampleCreator.py: library of functions for generating unwrapped and wrapped examples and baselines
Author: Emil Haaber Tellefsen
Co-Authors: Linus Ravn Gudmundsson, Niels Schøtt Hvidberg

Date: 15/11/2024
"""

# -- Third party --
import numpy as np


def createExample(unwrapped_baseline:  np.ndarray, add_image:  np.ndarray = None , mult_image:  np.ndarray = None, noise_sigma: float = 0.0, add_first: bool = True):
    """
    Creating unwrapped and wrapped case example for demonstrating Goldsteins algorithm, through adding and multiplying images on top of an unwrapped baseline.

    Parameters
    ----------
    unwrapped_baseline : array_like
        (n x m) Float array of the unwrapped baseline image which is to be modified
    add_image: array_like, optional
        (n x m) Float array of values to be added to the baseline. Set to 0 for all per default.
    mult_image: array_like, optional
        (n x m) Float array of values to be multiplied to the baseline. Set to 1 for all per default.
    noise_sigma: float, optional
        Standard deviation of applied gaussian noise on image. Set to 0 per default meaning no noise is applied.
    add_first: bool, optional
        Denotes whether add_image is added before mult_image is multiplied to baseline, or not. Default is True.

    Returns
    -------
    unwrapped_phase: array_like
        (n x m) Float array with altered unwrapped phase.

    wrapped_phase: array_like
        (n x m) Float array with adapted wrapped phase.

    """
    N, M = unwrapped_baseline.shape

    # Creating default images if nothing is specified
    if add_image is None:
        add_image = np.zeros((N,M))
    if mult_image is None:
        mult_image = np.ones((N,M))
    
    # creating noise
    noise = np.random.randn(N,M)*noise_sigma

    # applying arrays
    if add_first:
        phase_unwrapped = (unwrapped_baseline+add_image)*mult_image + noise
    else:
        phase_unwrapped = unwrapped_baseline*mult_image + add_image + noise
    
    # wrapping
    phase_wrapped = np.mod(phase_unwrapped, 2*np.pi)

    return phase_unwrapped, phase_wrapped



def createUnwrappedBaseline(shape: tuple, UnwrappedPhaseLimists: tuple, format: str = 'diagonal'):
    """
    Making baseline unwrapped phase given image shape, value limits of unwrapped phase, and template format. 
    Five formats are specified:
    - 'diagonal': phase increases from upper left corner to lower right.
    - 'horizontal': phase increases from left  to right.
    - 'vertical': phase increases from up  to down.
    - 'parabola_peak': 2D quadratic polynomial with max value in center and min value in corner.
    - 'parabola_valley': 2D quadratic polynomial with min value in center and max value in corner.

    Parameters
    ----------
    shape : tuple
        (nrow,ncolumn) int values for output image dimensions
    UnwrappedPhaseLimists: tuple
        (minphase,maxphase) float values for minimum and maximum phase
    format:
        format for output as specified above.

    Returns
    -------
    array_like
        (n x m) Float array of the unwrapped baseline image
    """   
    if format in ['diagonal', 'horizontal', 'vertical']:
        # making linear meshgrid
        x = np.linspace(UnwrappedPhaseLimists[0],UnwrappedPhaseLimists[1],shape[0])
        y = np.linspace(UnwrappedPhaseLimists[0],UnwrappedPhaseLimists[1],shape[1])
        xv, yv = np.meshgrid(x,y)

        # returning chosen format
        if format == 'diagonal':
            return (xv+yv)/2
        
        elif format == 'horizontal':
            return xv
        
        elif format == 'vertical':
            return yv
    
    elif format in ['parabola_peak', 'parabola_valley']:
        center_row, center_col = shape[0] // 2, shape[1] // 2

        # Create a grid of indices
        y, x = np.ogrid[:shape[0], :shape[1]]

        # Calculate squared distances from the center
        dist_squared = (y - center_row)**2 + (x - center_col)**2

        # Normalize the distances so corners are at max distance
        max_dist_squared = center_row**2 + center_col**2
        normalized_dist = dist_squared / max_dist_squared

        if format == 'parabola_peak':
            return UnwrappedPhaseLimists[1] - (UnwrappedPhaseLimists[1] - UnwrappedPhaseLimists[0]) * normalized_dist
        elif format == 'parabola_valley':
            return UnwrappedPhaseLimists[0] - (UnwrappedPhaseLimists[0] - UnwrappedPhaseLimists[1]) * normalized_dist
    
    else:
        raise ValueError('Unknown format')


def TrueBranchCuts(phase_unwrapped: np.ndarray):
    """
    Creating the theoretical "True Branch Cuts" given the known unwrapped phase

    Parameters
    ----------
    phase_unwrapped : array_like
        (n x m) Float array with altered unwrapped phase.

    Returns
    -------
    array_like
        (n x m) boolean array where theoretical branchcuts are True
    """      
    dhorz = phase_unwrapped[:,1:] - phase_unwrapped[:,:-1]
    dvert = phase_unwrapped[1:,:] - phase_unwrapped[:-1,:]
    BC = np.zeros(phase_unwrapped.shape, dtype=bool)
    BC[:,:-1][abs(dhorz)>np.pi] = True
    BC[:-1,:][abs(dvert)>np.pi] = True
    return BC

