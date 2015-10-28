# -*- coding: utf-8 -*-
"""
Setup extension helpers.
"""
from __future__ import absolute_import
import sys
import os
import glob
from distutils.core import Extension
from astropy_helpers import setup_helpers

__all__ = ['get_external_libraries', 'get_libFTR_extensions']

library_name = "libFTR"

def get_external_libraries():
    return ['ftr']
    
def get_libFTR_include_directory():
    """Get libFTR include directory."""
    root_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    include_directory = os.path.join(root_directory, 'cextern', library_name, 'src')
    return include_directory

def get_libFTR_extensions(filename, modulename, pattern="*.pyx", **kwargs):
    """docstring for get_libFTR_extensions"""
    
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
        cname = os.path.splitext(os.path.basename(component))[0]
        if cname.startswith("_"):
            cname = cname[1:]
            name = ".".join(this_name + ["_{0:s}".format(cname)])
        else:
            name = ".".join(this_name + [cname])
        extension_args['sources'] = [component]
        
        # Library checks.
        if setup_helpers.use_system_library('ftr'):
            libraries.append('ftr')
        else:
            extension_args['sources'].extend(glob.glob(os.path.join(include_directory, "*" + cname + "*.c")))
        # Extension object.
        extension = Extension(name, **extension_args)
        extensions.append(extension)
    
    return extensions