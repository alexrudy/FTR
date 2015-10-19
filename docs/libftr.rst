.. highlight:: c

.. _libftr:

*******************************
libFTR: C implementation of FTR
*******************************

This is a library implementation of FTR using pure C, suitable for inclusion in additional code.

The C library is meant to be more light-weight than the python implementation,
and does not provide functions to generate filters out of the box. It is designed to do the reconstruction work in-place, so that arrays can be allocated only once for a single reconstructor.

The design is meant to mimic the design of FFTW plans. This was chosen because FTRlib uses FFTW to implement the fast Fourier transform, making the plan architecture a reasonable choice.

Integrating libFTR
==================

If you have a pointer to x slopes, a pointer to y slopes, and you want to create a function which does the reconstruction once do::

    #include <ftr.h>
    ftr_plan recon;
    int nx = 10, ny = 10;
    double * sx, * sy, * est;
    est = malloc(sizeof(double) * nx * ny);
    sx = malloc(sizeof(double) * nx * ny);
    sy = malloc(sizeof(double) * nx * ny);

    // Set up the reconstructor
    recon = ftr_plan_reconstructor(nx, ny, sx, sy, est);

    // Set the filter
    ftr_set_filter(recon, gx, gy);

    // Do the reconstruction as many times as necessary.
    ftr_reconstruct(recon);

    // Free allocated memory at the end.
    ftr_plan_destroy(recon);

At the end of this code snippet, the reconstructed phase is always stored in the ``est`` variable defined above, and the operation is in place. Subsequent reconstructions should replace the contents of ``sx`` and ``sy`` rather than allocating a new reconstructor via :c:func:`ftr_plan_reconstructor`.

To change the filter used by a particular :c:type:`ftr_plan`, it is safe to call :c:func:`ftr_set_filter`.

Fourier Transform Reconstructor API
===================================

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

.. c:type:: ftr_estimate_callback

    This is the callback type for functions which can serve as callbacks in :c:func:`ftr_reconstruct_with_callback`. The callback signature must match ``void (*ftr_estimate_callback)(void * data, fftw_complex * est_ft)``. The inclusion of the ``void * data`` pointer allows for an arbitrary structure of user data to be passed in to the Fourier Transform Reconstructor.

.. c:function:: void ftr_reconstruct_with_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)

    :param ftr_plan recon: The :c:type:`ftr_plan` struct for this reconstructor.
    :param ftr_estimate_callback callback: A callback function to be applied to the fourier transform of the phase estimate.
    :param void* data: A pointer to data required by `callback`.

.. c:function:: void ftr_destroy

C Impelmentation of Slope Management for non-periodic domains
=============================================================

Slope management corrects finite aperture slope measurements for the Fourier transform reconstructor. Fourier transforms are implicitly carried out on a fully periodic domain. This assumption does not hold when looking at a typical wavefront sensor. For more details about slope management, see :ref:`slopemanagement`.

A minimial example of slope management::

    #include <slopemanage.h>
    sm_plan plan;
    int nx = 10, ny = 10;
    double * sx, * sy;
    int * ap, i, j;

    sx = malloc(sizeof(double) * nx * ny);
    sy = malloc(sizeof(double) * nx * ny);
    ap = malloc(sizeof(int) * nx * ny);

    // Set up an aperture with a border.
    // At least one border row/column is required for slope
    // management, so that there is enough room to put the
    // fixed slope values.
    for (i = 0; i < nx; ++i)
    {
        for (j = 0; j < ny; ++j)
        {
            if(i == 0 || j == 0 || i == nx - 1 || j == ny - 1)
            {
                ap[i + (j * nx)] = 0;
            }else{
                ap[i + (j * nx)] = 1;
            }
        }
    }


    // Set up the slope management plan
    plan = slope_management_plan(nx, ny, ap);

    // Do the slope managmenet as many times as necessary.
    slope_management_execute(plan, sy, sx);
    // Unlike FTR plans, slope management plans can operate
    // on different arrays each time.

    // Free allocated memory at the end.
    slope_management_destroy(plan);


Slope Management API
====================

.. c:type:: sm_plan

    The slope management plan, which contains the memory allocation for a single slope management scheme. The plan is generated by :c:func:`slope_management_plan`, and is an opaque structure containing the relevant pointers for performing slope management.

.. c:function:: sm_plan slope_management_plan(int ny, int nx, int *ap)

    Prepare a slope management scheme. This function creates a :c:type:`sm_plan` object which contains the memory allocation for the slope management scheme.

    :param int ny: Number of y positions (rows).
    :param int nx: Number of x positions (columns).
    :param int* ap: A pointer to the aperture (which should be `nx` by `ny` in size). Apertures are defined as 1 where light is transmissive.
    :returns: :c:type:`sm_plan`, the slope management plan.

.. c:function:: void slope_management_execute(sm_plan plan, double * sy, double * sx)

    Execute the slope managment plan, adjusting slopes in the `sx` and `sy` pointers. This method adjusts slopes in-place.

    :param sm_plan plan: The slope management plan to execute. A :c:type:`sm_plan` can be created using :c:func:`slope_management_plan`.
    :param double* sy: A pointer to the y slopes, as an array. Must conform to the dimensions set during the planning process.
    :param double* sx: A pointer to the x slopes, as an array. Must conform to the dimensions set during the planning process.
    :returns: No return value is provided, as the function acts on `sx` and `sy` in place.

.. c:function:: void slope_management_destroy(sm_plan plan)

    Deallocate the slope management plan. Memory allocated using :c:func:`slope_management_plan` will be freed.

    :param sm_plan plan: The slope management plan to execute. A :c:type:`sm_plan` can be created using :c:func:`slope_management_plan`.

.. c:function:: void slope_management(int ny, int nx, int *ap, double * sy, double * sx)

    Conduct the entire slope management process in a single call. This will allocate memory and determine aperture settings on the fly.

    Slope management happens in-place on the original arrays. No copy is performed.

    :param int ny: Number of y positions (rows).
    :param int nx: Number of x positions (columns).
    :param int* ap: A pointer to the aperture (which should be `nx` by `ny` in size). Apertures are defined as 1 where light is transmissive.
    :param double* sy: A pointer to the y slopes, as an array. Must conform to the dimensions set during the planning process.
    :param double* sx: A pointer to the x slopes, as an array. Must conform to the dimensions set during the planning process.
    :returns: No return value is provided, as the function acts on `sx` and `sy` in place.

