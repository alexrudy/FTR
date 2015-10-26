# -*- coding: utf-8 -*-
from __future__ import absolute_import
import sys
import os
import glob
from distutils.core import Extension
from astropy_helpers import setup_helpers

def get_extensions():
    """Return appropriate cython extensions"""
    # Some basic constants
    library_name = "libFTR"
    
    this_directory = os.path.dirname(__file__)
    this_name = __name__.split(".")[:-1]
    
    root_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    include_directory = os.path.join(root_directory, 'cextern', library_name, 'src')
    
    libraries = ['fftw3', 'fftw3_threads', 'pthread']
    
    extension_args = {
        'include_dirs' : ['numpy', include_directory ],
        'libraries' : libraries,
        'sources' : []
    }
    
    extensions = []
    
    for component in glob.iglob(os.path.join(this_directory, '_*.pyx')):
        # Component name and full module name.
        cname = os.path.basename(component)[1:-len(".pyx")]
        name = ".".join(this_name + ["_{0:s}".format(cname)])
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