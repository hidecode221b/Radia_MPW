# ----------------------------------------------------------
# PyRadiaUndulators library
# radia_util.py
# Useful methods
#
# Gael Le Bec, ESRF, 2019
# ----------------------------------------------------------

import radia as rad


# -------------------------------------------
def save(obj, filename, ext='.rad', path=''):
    """
    Save the object to filename.rad (another extensions can be specified)
    Use the load() method for loading the object.

    The file can be  opened with the Mathematica interface using the RadUtiLoad[] function

    :param filename: file name without extension
    :param ext: file extension (default = '.rad')
    :param path: absolute path if specified (default = '', i.e. relative path)
    :return: True
    """
    # --- Dump the undulator
    dmp = rad.UtiDmp(obj, 'bin')
    # --- Write to a file
    f = open(path + filename + ext, 'wb')
    f.write(dmp)
    f.close()

    return True

# -------------------------------------------
def load(filename, ext='.rad', path=''):
    """
    Load a Radia object from a file with .rad extension (or other extension if specified)
    :param filename: file name without extension
    :param ext: file extension (default = '.rad')
    :param path: absolute path if specified (default = '', i.e. relative path)
    :return: a Radia object
    """
    # --- Read the file
    f = open(path + filename + ext, 'rb')
    und_bin = f.read()
    f.close()
    # --- Build a Radia objec
    obj = rad.UtiDmpPrs(und_bin)
    # --- Return
    return obj
