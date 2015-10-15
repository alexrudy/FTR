Simple Fourier Transform Reconstructor
**************************************

The simple Fourier Transform Reconstructor uses the Fourier transform to provide an estimate of the phase of an adaptive optics system from the measured slope values.

Using this module only requires you to provide an aperture to the reconstructor class, and then call the reconstructor class::
    
    >>> ap = np.ones((10,10))
    >>> reconstructor = FourierTransformReconstructor(ap, filter='mod_hud')
    >>> sx = np.random.randn(ap.shape)
    >>> sy = np.random.randn(ap.shape)
    >>> est = reconstructor(sx, sy)
    

The Fourier Transform reconstructor assumes that the estimated phase is occuring on a square, fully-periodic grid. In order to correct for this assumption, it is better to use a technique called "Slope Management". Reconstructors that incorporate slope management are provided in :mod:`FTR.slopemanage`.

There are two similar reconstructor classes which users might want to take advantage of in this module, the :class:`FourierTransformReconstructor` and :class:`FastFTReconstructor`. The first class provides a pure-python implementation of the Fourier Transform reconstructor, and is generally easier to debug, in that it ensures that data passed wont cause problems. The second class, :class:`FastFTReconstructor`, uses :ref:`libFTR` under the hood, and so achieves fast speeds leveraging the FFTW library and the speed of pure c-code. However, it does less work to ensure variable saftey, and so is not a good choice for debugging.

Reference / API
===============

.. automodapi:: FTR.ftr
