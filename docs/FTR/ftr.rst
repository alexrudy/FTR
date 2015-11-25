Simple Fourier Transform Reconstructor
**************************************

The simple Fourier Transform Reconstructor uses the Fourier transform to provide an estimate of the phase of an adaptive optics system from the measured slope values.

Using this module only requires you to provide an aperture to the reconstructor class, and then call the reconstructor class::

    >>> import numpy as np
    >>> from FTR import FourierTransformReconstructor

    >>> ap = np.ones((10,10))
    >>> reconstructor = FourierTransformReconstructor(ap, filter='mod_hud')
    >>> sx = np.random.randn(*ap.shape)
    >>> sy = np.random.randn(*ap.shape)
    >>> est = reconstructor(sx, sy)


The Fourier Transform reconstructor assumes that the estimated phase is occuring on a square, fully-periodic grid. In order to correct for this assumption, it is better to use a technique called "Slope Management". Reconstructors that incorporate slope management are provided in :mod:`FTR.slopemanage`.

There are two similar reconstructor classes which users might want to take advantage of in this module, the :class:`FourierTransformReconstructor` and :class:`FastFTReconstructor`. The first class provides a pure-python implementation of the Fourier Transform reconstructor, and is generally easier to debug, in that it ensures that data passed wont cause problems. The second class, :class:`FastFTReconstructor`, uses :ref:`libFTR` under the hood, and so achieves fast speeds leveraging the FFTW library and the speed of pure c-code. However, it does less work to ensure variable saftey, and so is not a good choice for debugging.

Reconstruct and Invert Phase
============================

Using the Fourier Transform Reconstructor only requires some description of the illuminated aperture of your system. The aperture is specified as a nonzero numpy array::

    >>> import numpy as np
    >>> from FTR.utils import circle_aperture
    >>> from FTR import FourierTransformReconstructor
    
    >>> ap = circle_aperture((10,10), r=4.0)
    >>> reconstructor = FourierTransformReconstructor(ap, filter='mod_hud')
    >>> sx = np.random.randn(*ap.shape)
    >>> sy = np.random.randn(*ap.shape)
    >>> est = reconstructor(sx, sy)
    
The illuminated aperture is used for tip and tilt mode filtering. Tip-tilt filtering is controlled by the ``suppress_tt`` and ``manage_tt`` parameters, which can be specified when constructing an instance of :class:`FourierTransformReconstructor` or set as attributes on the reconstructor object::
    
    >>> import numpy as np
    >>> from FTR.utils import circle_aperture
    >>> from FTR import FourierTransformReconstructor
    
    >>> ap = circle_aperture((10,10), r=4.0)
    >>> reconstructor = FourierTransformReconstructor(ap, filter='mod_hud', suppress_tt=True)
    >>> reconstructor
    <FourierTransformReconstructor (10x10) filter='mod_hud' tt='suppressed'>
    >>> sx = np.random.randn(*ap.shape)
    >>> sy = np.random.randn(*ap.shape)
    >>> est = reconstructor(sx, sy)
    
    >>> reconstructor.suppress_tt = False
    >>> reconstructor.manage_tt = True
    >>> reconstructor
    <FourierTransformReconstructor (10x10) filter='mod_hud' tt='managed'>
    

Reference / API
===============

.. automodapi:: FTR.ftr
