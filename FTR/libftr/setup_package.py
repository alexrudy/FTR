# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
from distutils.core import Extension


def get_extensions():
    """Return appropriate cython extensions"""
    
    ftr_library_name = "libFTR"
    
    # Get include directory relative to setup.py
    ftr_include = os.path.abspath(os.path.join('cextern',ftr_library_name,'src'))
    ftr_library = os.path.abspath(os.path.join('cextern',ftr_library_name))
    local_dir = os.path.dirname(__file__)
    local_name = ".".join(["FTR", os.path.basename(os.path.dirname(__file__)), "_ftr"])
    
    # Set up sources.
    sources = [os.path.join(local_dir,"_ftr.pyx")]
    sources.extend([ os.path.join(ftr_include, filename) for filename in os.listdir(ftr_include) if filename.endswith(".c") ])
    
    # Set up libraries
    libraries = ['fftw3', 'fftw3_threads', 'pthread']
    
    _ftr = Extension(
        local_name, sources, include_dirs=['numpy', ftr_include], libraries=libraries
    )
    
    local_name = ".".join(["FTR", os.path.basename(os.path.dirname(__file__)), "_slopemanage"])
    
    # Set up sources.
    sources = [os.path.join(local_dir,"_slopemanage.pyx")]
    sources.extend([ os.path.join(ftr_include, filename) for filename in os.listdir(ftr_include) if filename.endswith(".c") ])
    
    _slopemanage = Extension(
        local_name, sources, include_dirs=['numpy', ftr_include], libraries=libraries
    )
    
    # Cython _ftr extension
    ext_modules = [
        _ftr, _slopemanage
    ]
    return ext_modules