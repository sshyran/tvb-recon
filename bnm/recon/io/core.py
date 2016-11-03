
"""
Core I/O.

"""

import numpy as np

try:
    from cStringIO import StringIO
except ImportError: # Py 3
    from io import BytesIO as StringIO


def np_save_strio(arr, fmt):
    "Save a NumPy array to text in a StringIO object."
    sio = StringIO()
    np.savetxt(sio, arr, fmt)
    return sio
