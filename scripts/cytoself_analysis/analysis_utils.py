
import numpy as np

def autoscale(im, percentile=None, dtype='uint8', gamma=None):
    '''
    Set the blackpoint and whitepoint of an image,
    and scale the values to span the dynamic range of `dtype`

    im : numpy array (any dimension)
    percentile : percentile threshold used for choosing the black and white points
    dtype : dtype of the scaled image (float, uint8, or uint16)
    gamma : optional gamma to apply to the scaled image
    '''

    max_values = {'float': 1.0, 'uint8': 255, 'uint16': 65535}

    if percentile is None:
        percentile = 0

    im = im.copy().astype(float)
    minn, maxx = np.percentile(im.flatten(), (percentile, 100 - percentile))
    if minn == maxx:
        return (im * 0).astype(dtype)

    im = im - minn
    im[im < 0] = 0
    im = im/(maxx - minn)
    im[im > 1] = 1
    if gamma is not None:
        im = im**gamma

    im = (im * max_values[dtype]).astype(dtype)
    return im
