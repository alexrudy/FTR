FFT Normalization
*****************

Different FFT implementations often have different normalization parameters. Most FFTs will be defined such that a forward transform follwed by an inverse transform will result in the same values. However, implementations tend to apply the normalization at different points. This document describes the normalizations applied by each FFT, and their implications for the Fourier transform reconstructor.

The implementation of the FTR is taken from Lisa Poyneer's IDL implementation (which is not) included in this distribution. However, we can use the IDL FFT definition to clarify our FFT normalization. This document shows how to modify other FFT normalizations so that they match the IDL normalization, the canonical one used in this algorithm. [#fIDL]_

.. [#fIDL] There is nothing *more* correct about the IDL implementation, it just happens to be the canonical implementation as it is the referernce implementaiton provided to match Lisa's thesis.

IDL Normalization
-----------------

IDL defines the FFT (`here <http://www.physics.nyu.edu/grierlab/idl_html_help/F4.html>`_) as

.. math::
    
    F(u) = \frac{1}{N}\sum_{x = 0}^{N - 1} f(x) \exp\left[-j \frac{2 \pi u x}{N}\right]

and the inverse FFT as

.. math::
    
    f(x) = \sum_{u = 0}^{N - 1} F(u) \exp\left[j \frac{2 \pi u x}{N}\right]
    
The important part of this definition for our purposes is the factor of :math:`\tfrac{1}{N}` which is applied to the forward transform.

When we apply the spatial filters used in the Fourier transform reconstructor, we apply two forward transforms (one for x, one for y), and one backwards transform. Collecting only the factors of :math:`\tfrac{1}{N}`, we get two factors due to the two forward transforms which are added together. No factor is included in the inverse transform.

The normalization factor for the full FTR filter is :math:`\tfrac{2}{N}`.

Numpy Normalization
-------------------

:mod:`numpy` defines the forward FFT (`see here <http://docs.scipy.org/doc/numpy/reference/routines.fft.html#implementation-details>`_) as

.. math::
   A_k =  \sum_{m=0}^{n-1} a_m \exp\left\{-2\pi i{mk \over n}\right\}
   \qquad k = 0,\ldots,n-1.
   
and the inverse FFT as

.. math::
   a_m = \frac{1}{n}\sum_{k=0}^{n-1}A_k\exp\left\{2\pi i{mk\over n}\right\}
   \qquad m = 0,\ldots,n-1.
   
Comparing this definition to the one used by IDL, we can see that the normalization is applied as a :math:`\tfrac{1}{n}` factor in the inverse transform.

Comparing to the IDL implementation, when we apply the two forward transforms, then add them, and apply a single reverse transform, we find that the total normalization is :math:`\tfrac{1}{n}`. Therefore, we should multiply the Fourier transform of the phase estimate by 2 to match the normalization used by IDL.


  