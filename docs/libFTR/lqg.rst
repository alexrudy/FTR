.. highlight:: c
.. default-domain:: c

Linear–quadratic–Gaussian control in libFTR
*******************************************

LQG Control in libFTR is implemented as an LQG controller in the Fourier domain.

Using the LQG filter
====================

Filters in libFTR use a plan architecture, in the same vein as :ref:`libFTR's reconstructors <libftr-ftr>` .The LQG Filter in libFTR is implemented as an in-place function which acts on a complex array (using the FFTW complex type :type:`fftw_complex`.)

The user is responsible for allocating and providing the input array, and the filtering occurs in-place on the input array.

Given a pointer to an array of complex values to filter, to use the modal LQG::

    #include <lqg.h>
    lqg_filter filter;
    // A 10x10, 2-layer filter
    int nx = 10, ny = 10, nl = 2;
    fftw_complex *est_ft;
    fftw_complex *gain, *alpha, *hp_coeffs;

    // Phase estimate array.
    est_ft = fftw_malloc(sizeof(fftw_complex) * nx * ny);

    // Generate the filter coefficients.
    gain = fftw_malloc(sizeof(fftw_complex) * nx * ny * nl);
    alpha = fftw_malloc(sizeof(fftw_complex) * nx * ny * nl);
    hp_coeffs = fftw_malloc(sizeof(fftw_complex) * nx * ny);
    // You need to provide values for gain, alpha and hp_coeffs.

    // Set up the filter
    filter = lqg_new_filter(nl, ny, nx, gains, alphas, hp_coeffs);

    // Reset the filter
    lqg_reset(filter);

    // Do the reconstruction as many times as necessary.
    lqg_apply_filter(filter, est_ft);

    // Free allocated memory at the end.
    lqg_destroy(filter);

At the end of this code snippet, the filtered phase is always stored in the ``est_ft`` variable defined above, and the operation is in place. Subsequent reconstructions should not allocate a new filter via :func:`lqg_new_filter`. When the integrators in the LQG filter must be reset, use :func:`lqg_reset`

To change the coefficients used by the LQG filter, you can change them directly in memory.

LQG Filter API
==============

The API relies heavily on the concept of a :type:`lqg_filter`. The filter is an opaque pointer to an internal structure which manages the geometry and memory involved in performing the filtering and retaining integrator state. Filters are created with :func:`lqg_new_filter`, and destroyed with :func:`lqg_destroy`.

.. type:: lqg_filter

    This is the core persistance type for LQG filters. It is a struct which can be treated as an opaque object to the user, which maintains pointers to the re-used variables in the LQG filter process.

Creating and destroying filters
-------------------------------

.. function:: lqg_filter lqg_new_filter(const int nl, const int ny, const int nx, fftw_complex * gains, fftw_complex * alphas, fftw_complex * hp_coefficients)

    This function allocates a :type:`lqg_filter` struct with the correct members.

    :param int nl: The number of LQG filter layers.
    :param int nx: The number of points in the x direction.
    :param int ny: The number of points in the y direction.
    :param fftw_complex gains: A pointer to the gain coefficients (nx by ny by nl).
    :param fftw_complex alphas: A pointer to the alphas (nx by ny by nl).
    :param fftw_complex hp_coefficients: A pointer to the high pass coefficients (nx x ny).
    :returns: A :type:`lqg_filter` with allocated data arrays.

.. function:: void lqg_destroy(lqg_filter filter)

    Destroy an LQG filter, deallocating memory as necessary. The user is still responsible for freeing the coefficient arrays (`gains`, `alphas` and `hp_coefficients`).

    :param lqg_filter filter: The :type:`lqg_filter` to deallocate.

Applying and managing filters
-----------------------------

.. function:: void lqg_apply_filter(lqg_filter filter, fftw_complex * est_ft)

    Perform the LQG filter. Filter results are computed in-place in ``est_ft``.

    :param lqg_filter filter: The :type:`lqg_filter` struct for this filter.
    :param fftw_complex est_ft: A pointer to the estimated phase data, which will be overwritten with the result.

.. function:: void lqg_reset(lqg_filter, filter)

    Reset the LQG filter integrator internal states, by zeroing the arrays.

    :param lqg_filter filter: The :type:`lqg_filter` to reset.

Integrating the LQG filter with the Fourier Transform Reconstructor
-------------------------------------------------------------------

:ref:`libftr-ftr` accepts a callback function to apply a filter or series of filters. The callback is fairly general, but it can be used to include an LQG filter in the middle of the Fourier Transform Reconstructor. The function :func:`lqg_filter_callback` satisfies the :type:`ftr_estimate_callback` required to be a callback in the middle of the Fourier Transform Reconstructor. To use this function, simply cast the :type:`lqg_filter` to a void pointer, ``filter = <void *>my_lqg_filter`` and use it as the `data` argument to :func:`ftr_reconstruct_with_callback`.

.. function:: void lqg_filter_callback(void * filter, const int ny, const int nx, fftw_complex * est_ft)
    
    :param void * filter: The LQG filter, should be a :type:`lqg_filter` item.
    :param int nx: The number of points in the x direction.
    :param int ny: The number of points in the y direction.
    :param fftw_complex * est_ft: A pointer to the estimated phase. Filtering takes place in place in the `est_ft` array.
    

