import math
import numpy as np
import matplotlib.pyplot as plt


def gauss(x, a, b, c):
    """
    Creates values of gaussian function:
    f(x) = ae**(((x-b)**2)/(2c**2))
    Takes as x list of values
    Returns list of y
    Requires coeffs a, b, c
    Also takes h if this value exists
    """
    return a*np.exp(-( ( (x-b) ** 2)/(2 * (c**2))) )
