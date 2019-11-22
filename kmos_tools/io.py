# -*- coding: utf-8 -*-
"""I/O functions for KMOS data


"""

__author__ = ["Charlotte Mason (UCLA)"]

import glob


def insensitive_glob(pattern):
    """Case insensitive find file names

    Args:
        pattern (str): pattern to search for

    Returns:
        file_list (list): sorted list of filenames which match pattern
    """
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c

    file_list = sorted(glob.glob(''.join(map(either, pattern))))

    return file_list


def find_exposures(dir='*', prefix='*'):
    """Find all the SCI_RECONSTRUCTED fits files in a directory

    Args:
        dir (str): Directory to look for exposures in
        prefix (str): prefix for filenames

    """
    dir_prefix = '%s/%s' % (dir, prefix)

    sci_exposures = set(glob.glob('%s*SCI_RECONSTRUCTED*.fits' % dir_prefix)) \
                - set(glob.glob('%s*star*.fits' % dir_prefix)) \
                - set(glob.glob('%s*SKYFIX*.fits' % dir_prefix)) \
                - set(glob.glob('%s*SKYSPEC*.fits' % dir_prefix)) \
                - set(glob.glob('%s*MASK*.fits' % dir_prefix)) \
                - set(glob.glob('%s*COMBINE*.fits' % dir_prefix)) \
                - set(glob.glob('%s*COLL*.fits' % dir_prefix))  \
                - set(glob.glob('%s*BADCALIB*.fits' % dir_prefix))  \
                - set(glob.glob('%s*FLUXFIX*.fits' % dir_prefix))

    print(len(sci_exposures), 'science exposures found')

    return sci_exposures
