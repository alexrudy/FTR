# -*- coding: utf-8 -*-
"""
Setup extension helpers.
"""
from __future__ import absolute_import
import sys
import os
import glob
import copy

__all__ = ['get_external_libraries', 'get_libFTR_extensions']

library_name = "libFTR"

def get_external_libraries():
    return ['ftr']
    
def get_libFTR_include_directory():
    """Get the libFTR include directory, relative to the standard location of setup.py
    
    The path is roughly ./cextern/libFTR/src
    
    """
    root_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    include_directory = os.path.join(root_directory, 'cextern', library_name, 'src')
    return include_directory

def get_libFTR_extensions(filename, modulename, pattern="*.pyx", **kwargs):
    """Make distutils extensions for libFTR files in the FTR tree.
    
    :param filename: Filename to be used as a root for finding extensions. Often ``__file__``.
    :param modulename: Module name to be used as an extension root. Often ``__name__``.
    :param pattern: A glob pattern to be used to find extension source files.
    :kwargs: Remaining keyword arguments provide the default values for the distutils extension.
    :returns: Return a list of distutils extensions.
    
    """
    from distutils.core import Extension
    from astropy_helpers import setup_helpers
    
    this_directory = os.path.dirname(filename)
    this_name = modulename.split(".")[:-1]
    
    include_directory = get_libFTR_include_directory()
    libraries = ['fftw3', 'fftw3_threads', 'pthread']
    
    extension_args = {
        'include_dirs' : ['numpy', include_directory ],
        'libraries' : libraries,
        'sources' : []
    }
    extension_args.update(kwargs)
    
    extensions = []
    
    for component in glob.iglob(os.path.join(this_directory, pattern)):
        # Component name and full module name.
        this_extension_args = copy.deepcopy(extension_args)
        cname = os.path.splitext(os.path.basename(component))[0]
        if cname.startswith("_"):
            cname = cname[1:]
            name = ".".join(this_name + ["_{0:s}".format(cname)])
        else:
            name = ".".join(this_name + [cname])
        this_extension_args['sources'].append(component)
        
        # Library checks.
        if setup_helpers.use_system_library('ftr'):
            this_extension_args['libraries'].append('ftr')
        else:
            this_extension_args['sources'].extend(glob.glob(os.path.join(include_directory, "*.c")))
            this_extension_args['sources'] = list(set(this_extension_args['sources']))
        # Extension object.
        extension = Extension(name, **this_extension_args)
        extensions.append(extension)
    
    return extensions