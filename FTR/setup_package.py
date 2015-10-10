# -*- coding: utf-8 -*-

import os

def get_extensions():
    """Return appropriate cython extensions"""
    from distutils.core import Extension
    
    # Get include directory relative to setup.py
    ftr_include = os.path.abspath(os.path.join('cextern','FTR','src'))
    ftr_library = os.path.abspath(os.path.join('cextern','FTR'))
    local_dir = os.path.dirname(__file__)
    local_name = ".".join([os.path.basename(os.path.dirname(__file__)), "_ftr"])
    
    # Set up sources.
    sources = [os.path.join(local_dir,"_ftr.pyx")]
    sources.extend([ os.path.join(ftr_include,filename) for filename in os.listdir(ftr_include) if filename.endswith(".c") ])
    
    # Set up libraries
    libraries = ['fftw3']
    
    # Cython _ftr extension
    ext_modules = [
        Extension(
            local_name, sources, include_dirs=['numpy', ftr_include], libraries=libraries
        )
    ]
    return ext_modules