Fourier Transform Reconstructor
===============================

This package implements the Fourier Transform Reconstructor and a few additional reconstruction tools. The Fourier Transform Reconstructor is a method for reconstructing a 2-D function from the first derivatives of that function, using the fast fourier transform.

The Fourier Transform Reconstructor is useful in adaptive optics, where a wavefront sensor measures the slope of the incoming wavefront, but the useful quantity to the adaptive optics system is the phase.

This module implements the Fourier Transform Reconstructor, as well as a few additional reconstruction tools and filters which are helpful for adaptive optics systems. There are two implementations, a C implementation, found in cextern/libFTR, and a python implementation in FTR. The python implementation also includes interfaces into the C implementation, via cython.

For real-time impelementations, using the C implemetnation is appropriate. For development and diagnostics, the python implementation is more appropriate.

Documentation is available at http://alexrudy.github.io/FTR/.

Available Reconstructors
------------------------

The python implementation contains implementations of the Fourier Transform Reconstructor and a simple vector-matrix-multiply reconstructor. The C library only contains an implementation of the Fourier Transform Reconstructor, but an example implementation of a vector-matrix-multiply reconstructor can be found in the examples in libFTR.

Additional Features
-------------------

Additionally, this library implements the slope-management techniques described by Lisa Poyneer in her dissertation. Slope management helps the Fourier Transform Reconstructor account for the finite aperture discontinuity between opposite edges.

Installation
------------

With ``pip``, you can install this package::
    
    $ pip install git+https://github.com/alexrudy/ftr.git
    

It requires ``numpy`` to run. ``astropy`` and ``h5py`` are required for optional I/O support for saving and loading filters.

Status reports for developers
-----------------------------

.. image:: https://travis-ci.org/alexrudy/FTR.svg?branch=master
    :target: https://travis-ci.org/alexrudy/FTR
    :alt: Test Status
