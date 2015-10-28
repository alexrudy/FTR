.. highlight:: c

.. _libftr:

*******************************
libFTR: C implementation of FTR
*******************************

This is a library implementation of FTR using pure C, suitable for inclusion in additional code.

The C library is meant to be more light-weight than the python implementation,
and does not provide functions to generate filters out of the box. It is designed to do the reconstruction work in-place, so that arrays can be allocated only once for a single reconstructor.

The design is meant to mimic the design of FFTW plans. This was chosen because FTRlib uses FFTW to implement the fast Fourier transform, making the plan architecture a reasonable choice.

.. toctree::
    :maxdepth: 2
    
    ftr
    slopemanagement
    aperture