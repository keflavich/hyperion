from __future__ import print_function, division

import string
import random
import os
import shutil
import sys
import tempfile

import h5py
import numpy as np

from .constants import h, c, k

TMPDIR = tempfile.mkdtemp()


def link_or_copy(group, path, link, copy, absolute_paths=False):
    '''
    Link or copy a dataaset or group

    Parameters
    ----------
    group: h5py.Group
        The group to create the link, dataset, or group in
    path: str
        The path of the link, dataset, or group in the new file
    link: h5py.ExternalLink
        A link to the group or dataset to include
    copy: bool
        Whether to copy or link to the dataset
    absolute_paths: bool
        If copy=False, then if absolute_paths is True, absolute filenames
        are used in the link, otherwise the path relative to the file is used.
    '''
    if copy:
        group.copy(h5py.File(link.filename, 'r')[link.path], path)
    else:
        if absolute_paths:
            group[path] = h5py.ExternalLink(os.path.abspath(link.filename), link.path)
        else:
            group[path] = h5py.ExternalLink(os.path.relpath(link.filename, os.path.dirname(group.file.filename)), link.path)


class FreezableClass(object):

    _frozen = False
    _final = False
    _attributes = []

    def _freeze(self):
        object.__setattr__(self, '_frozen', True)

    def _finalize(self):
        object.__setattr__(self, '_final', True)

    def isfrozen(self):
        return self._frozen

    def isfinal(self):
        return self._final

    def __setattr__(self, key, value):
        if self._final:
            raise Exception("Attribute %s can no longer be changed" % key)
        if self._frozen and not key in self._attributes:
            raise AttributeError("Attribute %s does not exist" % key)
        self._attributes.append(key)
        object.__setattr__(self, key, value)


def planck_nu_range(tmin, tmax=None):

    # Set constant for Wien displacement law
    alpha = 2.821439

    # Find peak frequencies for lower and upper temperature
    nu_peak_min = alpha / h * k * tmin
    if tmax is None:
        nu_peak_max = alpha / h * k * tmin
    else:
        nu_peak_max = alpha / h * k * tmax

    # Find frequency range if we want to go below 1 thousandth of the Planck
    # function peak

    nu_min = np.log10(nu_peak_min / 100.)
    nu_max = np.log10(nu_peak_max * 10.)

    # If we use at least 100 frequencies per order of magnitude of frequency
    # then we will achieve differences around 2% flux-to-flux

    n_nu = int((nu_max - nu_min) * 100.)

    return np.logspace(nu_min, nu_max, n_nu)


def nu_common(nu1, nu2):

    # Combine frequency lists/arrays and sort
    nu_common = np.hstack([nu1, nu2])
    nu_common.sort()

    # Remove unique elements (can't just use np.unique because also want
    # to remove very close values)
    keep = (nu_common[1:] - nu_common[:-1]) / nu_common[:-1] > 1.e-10
    keep = np.hstack([keep, True])
    nu_common = nu_common[keep]

    # Return common frequency range
    return nu_common


class extrap1d_log10(object):

    def __init__(self, x, y):

        self.x = np.log10(x)
        self.y = np.log10(y)

    def __call__(self, x):

        xval = np.log10(x)

        if self.x[-1] > self.x[0]:
            top = xval > self.x[-1]
            bot = xval < self.x[0]
        else:
            top = xval < self.x[-1]
            bot = xval > self.x[0]

        if type(xval) == np.ndarray:

            yval = np.zeros(xval.shape)

            yval[top] = 10. ** (self.y[-1] + (xval[top] - self.x[-1]) * (self.y[-1] - self.y[-2]) / (self.x[-1] - self.x[-2]))
            yval[bot] = 10. ** (self.y[0] + (xval[bot] - self.x[0]) * (self.y[0] - self.y[1]) / (self.x[0] - self.x[1]))

        else:

            if top:
                yval = 10. ** (self.y[-1] + (xval - self.x[-1]) * (self.y[-1] - self.y[-2]) / (self.x[-1] - self.x[-2]))
            elif bot:
                yval = 10. ** (self.y[0] + (xval - self.x[0]) * (self.y[0] - self.y[1]) / (self.x[0] - self.x[1]))
            else:
                raise Exception("xval should lie outside x array")

        return yval


def B_nu(nu, T):
    return 2. * h * nu ** 3. / c ** 2. / (np.exp(h * nu / k / T) - 1.)


def filename2fits(filename):

    ext = os.path.splitext(filename)[1]
    if ext in ['.par', '.dat', '.txt', '.ascii']:
        return filename.replace(ext, '.fits')
    else:
        raise Exception("Unknown extension: %s" % ext)


def filename2hdf5(filename):

    ext = os.path.splitext(filename)[1]
    if ext in ['.par', '.dat', '.txt', '.ascii']:
        return filename.replace(ext, '.hdf5')
    else:
        raise Exception("Unknown extension: %s" % ext)


def random_id(length=32):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def random_filename():
    return os.path.join(TMPDIR, random_id())


def create_dir(dir_name):
    delete_dir(dir_name)
    os.mkdir(dir_name)


def delete_dir(dir_name):

    if os.path.exists(dir_name):
        reply = raw_input("Delete directory %s? [y/[n]] " % dir_name)
        if reply == 'y':
            shutil.rmtree(dir_name)
        else:
            print("Aborting...")
            sys.exit()


def delete_file(file_name):

    if os.path.exists(file_name):
        reply = raw_input("Delete file %s? [y/[n]] " % file_name)
        if reply == 'y':
            os.remove(file_name)
        else:
            print("Aborting...")
            sys.exit()

def is_numpy_array(variable):
    return type(variable) in [np.ndarray,
                              np.core.records.recarray,
                              np.ma.core.MaskedArray]


def monotonically_increasing(array):
    for i in range(len(array)-1):
        if not array[i+1] > array[i]:
            return False
    return True
