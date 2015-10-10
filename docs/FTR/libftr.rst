.. highlight:: c

.. _libftr:

*********************************
C implementation of FTR in libFTR
*********************************

This is a library implementation of FTR using pure C, suitable for inclusion in additional code.

The C library is meant to be more light-weight than the python implementation,
and does not provide functions to generate filters out of the box. It is designed to do the reconstruction work in-place, so that arrays can be allocated only once for a single reconstructor.

The design is meant to mimic the design of FFTW plans. This was chosen because FTRlib uses FFTW to implement the fast Fourier transform, making the plan architecture a reasonable choice.

Integrating libFTR
==================

If you have a pointer to x slopes, a pointer to y slopes, and you want to create a function which does the reconstruction once do::
    
    ftr_plan recon;
    int nx = 10, ny = 10;
    double * sx, * sy, * est;
    est = malloc(sizeof(double) * nx * ny);
    
    // Set up the reconstructor
    recon = ftr_plan_reconstructor(nx, ny, sx, sy, est);
    
    // Set the filter
    ftr_set_filter(recon, gx, gy);
    
    // Do the reconstruction as many times as necessary.
    ftr_reconstruct(recon);
    
At the end of this code snippet, the reconstructed phase is always stored in the ``est`` variable defined above, and the operation is in place. Subsequent reconstructions should replace the contents of ``sx`` and ``sy`` rather than allocating a new reconstructor via :c:func:`ftr_plan_reconstructor`.

To change the filter used by a particular :c:type:`ftr_plan`, it is safe to call :c:func:`ftr_set_filter`.

Refrence / API
==============

.. c:type:: ftr_plan

    This is the core persistance type for reconstruction data. It is a struct which can be treated as an opaque object to the user, which maintains pointers to the re-used variables in the reconstruction process.

.. c:function:: ftr_plan ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est)

    This function allocates a :c:type:`ftr_plan` struct with the correct members.

    :param int nx: The number of points in the x direction.
    :param int ny: The number of points in the y direction.
    :param double sx: A pointer to the x slope data.
    :param double sy: A pointer to the y slope data.
    :param double est: A pointer to the estimated phase output data.
    :returns: A :c:type:`ftr_plan` with allocated data arrays.

    This function also serves to initialize the FFTW plans which will be used to perform the reconstruction.

.. c:function:: void ftr_set_filter(ftr_plan recon, fftw_complex *gx, fftw_complex *gy)

    This function sets the filter in the :c:type:`ftr_plan` struct to point to the provided filter arrays. It also computes the filter denominator.

    :param ftr_plan recon: The :c:type:`ftr_plan` struct for this reconstructor.
    :param complex gx: The x spatial filter.
    :param complex gy: The y spatial filter.

    Changing the values in ``gx`` and ``gy`` after calling this function will leave the incorrect denominator stored in the :c:type:`reconstructor` struct.

.. c:function:: void ftr_reconstruct(ftr_plan recon)

    Perform the reconstruction. Reconstruction results are stored in the data assigned to ``est`` with :c:func:`ftr_plan_reconstructor`.

    :param ftr_plan recon: The :c:type:`ftr_plan` struct for this reconstructor.
    

