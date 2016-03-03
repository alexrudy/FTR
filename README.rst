Fourier Transform Reconstructor
===============================

This package implements the Fourier Transform Reconstructor and a few additional reconstruction tools. The Fourier Transform Reconstructor is a method for reconstructing a 2-D function from the first derivatives of that function, using the fast Fourier transform.

The Fourier Transform Reconstructor is useful in adaptive optics, where a wavefront sensor measures the slope of the incoming wavefront, but the useful quantity to the adaptive optics system is the phase.

This module implements the Fourier Transform Reconstructor, as well as a few additional reconstruction tools and filters which are helpful for adaptive optics systems. There are two implementations, a C implementation, found in cextern/libFTR, and a python implementation in FTR. The python implementation also includes interfaces into the C implementation, via Cython.

For real-time implementations, using the C implementation is appropriate. For development and diagnostics, the python implementation is more appropriate.

Documentation is available at http://alexrudy.github.io/FTR/.

Python implementation
---------------------

Available Reconstructors
************************

The python implementation contains implementations of the Fourier Transform Reconstructor and a simple vector-matrix-multiply reconstructor.

Additional Features
*******************

Additionally, this library implements the slope-management techniques described by Lisa Poyneer in her dissertation. Slope management helps the Fourier Transform Reconstructor account for the finite aperture discontinuity between opposite edges.

C implementation
----------------

The C implementation, ``libFTR``, contains the Fourier Transform Reconstructor, an implementation of slope management, and a simple aperture management utility.

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
