.. highlight:: c
.. default-domain:: c

Fourier Transform Reconstructor in libFTR
*****************************************

Using the Fourier Transform Reconstructor
=========================================

The Fourier Transform Reconstructor uses a plan architecture, similar to what fftw_ uses. The intent is that the user can allocate the necessary memory and describe the operations required by the reconstructor once in a plan object, and then repeately re-use the plan to pefrom the reconstruction many times.

The user is still responsible for allocating the input and output arrays, to give the maximum flexibility inside a reconstruction loop. However, the underliyng fourier transforms performed with fftw_ will be fastest if the arrays are byte-aligned and SIMD compatible. The easiest way to do this is to allocate the input and output arrays with ``fftw_malloc``.

If you have a pointer to x slopes, a pointer to y slopes, and you want to create a function which does the reconstruction once do::

    #include <ftr.h>
    ftr_plan recon;
    int nx = 10, ny = 10;
    double * sx, * sy, * est;
    est = fftw_malloc(sizeof(double) * nx * ny);
    sx = fftw_malloc(sizeof(double) * nx * ny);
    sy = fftw_malloc(sizeof(double) * nx * ny);

    // Set up the reconstructor
    recon = ftr_plan_reconstructor(nx, ny, sx, sy, est);

    // Set the filter
    ftr_set_filter(recon, gx, gy);

    // Do the reconstruction as many times as necessary.
    ftr_reconstruct(recon);

    // Free allocated memory at the end.
    ftr_plan_destroy(recon);

At the end of this code snippet, the reconstructed phase is always stored in the ``est`` variable defined above, and the operation is in place. Subsequent reconstructions should replace the contents of ``sx`` and ``sy`` rather than allocating a new reconstructor via intfunc:`ftr_plan_reconstructor`.

To change the filter used by a particular inttype:`ftr_plan`, it is safe to call intfunc:`ftr_set_filter`.

Fourier Transform Reconstructor API
===================================

The API relies heavily on the concept of a :type:`ftr_plan`. The plan is an opaque pointer to an internal structure which manages the geometry and memory involved in performing the filtering and Fourier transforms. Plans are created with :func:`ftr_plan_reconstructor`, and destroyed with :func:`ftr_destroy`.

.. type:: ftr_plan

    This is the core persistance type for reconstruction data. It is a struct which can be treated as an opaque object to the user, which maintains pointers to the re-used variables in the reconstruction process.

Creating and destroying plans
-----------------------------

.. function:: ftr_plan ftr_plan_reconstructor(int nx, int ny, double *sx, double *sy, double *est)

    This function allocates a inttype:`ftr_plan` struct with the correct members.

    :param int nx: The number of points in the x direction.
    :param int ny: The number of points in the y direction.
    :param double sx: A pointer to the x slope data.
    :param double sy: A pointer to the y slope data.
    :param double est: A pointer to the estimated phase output data.
    :returns: A inttype:`ftr_plan` with allocated data arrays.

    This function also serves to initialize the FFTW plans which will be used to perform the reconstruction.

.. function:: void ftr_set_filter(ftr_plan recon, fftw_complex *gx, fftw_complex *gy)

    This function sets the filter in the inttype:`ftr_plan` struct to point to the provided filter arrays. It also computes the filter denominator.

    :param ftr_plan recon: The inttype:`ftr_plan` struct for this reconstructor.
    :param complex gx: The x spatial filter.
    :param complex gy: The y spatial filter.

    Changing the values in ``gx`` and ``gy`` after calling this function will leave the incorrect denominator stored in the inttype:`reconstructor` struct.

.. function:: void ftr_destroy(ftr_plan recon)

    Destroy an FTR plan, deallocating memory as necessary.

    :param ftr_plan recon: The inttype:`ftr_plan` to deallocate.

Reconstruction, with and without callbacks
------------------------------------------

.. function:: void ftr_reconstruct(ftr_plan recon)

    Perform the reconstruction. Reconstruction results are stored in the data assigned to ``est`` with intfunc:`ftr_plan_reconstructor`.

    :param ftr_plan recon: The inttype:`ftr_plan` struct for this reconstructor.

.. type:: ftr_estimate_callback

    This is the callback type for functions which can serve as callbacks in intfunc:`ftr_reconstruct_with_callback`. The callback signature must match ``void (*ftr_estimate_callback)(void * data, fftw_complex * est_ft)``. The inclusion of the ``void * data`` pointer allows for an arbitrary structure of user data to be passed in to the Fourier Transform Reconstructor.

.. function:: void ftr_reconstruct_with_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)

    Perform the reconstruction, and use a callback to adjust the fourier transform of the estimate. See inttype:`ftr_estimate_callback` for a descritpion of the callback function required to apply additional filters to the fourier transform of the estimate.

    :param ftr_plan recon: The inttype:`ftr_plan` for this reconstructor.
    :param ftr_estimate_callback callback: A callback function to be applied to the fourier transform of the phase estimate.
    :param void* data: A pointer to data required by `callback`.

Individual reconstruction steps
-------------------------------

.. function:: void ftr_forward_transform(ftr_plan recon)

    Perform only the forward FFTs to transform the slopes into the Fourier domain.

    :param ftr_plan recon: The inttype:`ftr_plan` to use for the forward transform.

.. function:: void ftr_apply_filter(ftr_plan recon)

    Only apply the filter to the transformed slopes, to estimate the phase.

    :param ftr_plan recon: The inttype:`ftr_plan` to use to apply the filter.

.. function:: void ftr_backward_transform(ftr_plan recon)

    Perform the backward FFT to transform the Fourier mode estimate of the phase into a real phase.

    :param ftr_plan recon: The inttype:`ftr_plan` to use for the backward transform.

.. function:: void ftr_apply_callback(ftr_plan recon, ftr_estimate_callback callback, void * data)

    Apply a callback function to the estimated phase in the Fourier domain.

    :param ftr_plan recon: The inttype:`ftr_plan` for this reconstructor.
    :param ftr_estimate_callback callback: A callback function to be applied to the fourier transform of the phase estimate.
    :param void* data: A pointer to data required by `callback`.

FFTW Halfcomplex Format Utilties
================================

.. function:: void ftr_map_half_complex(int ny, int nx, int * map, int * imap)

    Compute the mapping between half-complex transforms and the fully-expanded transform, for flattened arrays.

    In FFTW, the complex side of a real-to-complex transform (or vice-versa) does not include all the data points, and rather elimintates some of the data points which can be reconstructed based on the Hermitian symmetry of the output. See the FFTW documentation on `The Halfcomplex-format DFT <http://www.fftw.org/doc/The-Halfcomplex_002dformat-DFT.html#The-Halfcomplex_002dformat-DFT>`_ for more details on the exact specifics of this format. This function simply provides the mapping between a Halfcomplex array and a full array, as a pair of integer pointers.

    Using the pointer ``map``, you can map from a Halfcomplex array to a full array in 2D::

        int i;
        int *map, *imap;
        double *full_array, *half_array;
        full_array[i] = half_array[map[i]];
        half_array[i] = full_array[imap[i]];


    The pointers to `map` and `imap` should be allocated before calling this function.

    :param int ny: Number of y grid points.
    :param int nx: Number of x grid points.
    :param int* map: Mapping of full indices to halfcomplex indicies.
    :param int* imap: Mapping of halfcomplex indicies to full indices.

.. _fftw: http://www.fftw.org
