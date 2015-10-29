Slope Periodicity Management
****************************

The Fourier Transfrom Reconstructor assumes a fully periodic square grid for the phase. Most telescopes are not square, and atmospheric turbulence is not periodic, so we have to use a strategy to account for the non-periodicity in the slope grid. We do this by adding pseudo-slopes outside the telescope aperture, which serve to make the slope domain appear periodic, despite the non-periodic nature of the observed slopes.

Reconstruction on a periodic domain
===================================

Lets provide some example slopes on a fully periodic domain to work with.::

    >>> import numpy as np
    >>> sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    >>> sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2

These slopes are both periodic on the full fourier grid, and can be reconstructed with a pure Fourier Transform Reconstructor.

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt

    sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2
    f, (ax_x, ax_y) = plt.subplots(1, 2, figsize=(10,5))
    im = ax_x.imshow(sx, cmap='hot', vmin=1, vmax=3)
    ax_x.set_title("$S_x$")
    f.colorbar(im, ax=ax_x, fraction=0.046, pad=0.04)
    im = ax_y.imshow(sy, cmap='hot', vmin=1, vmax=3)
    ax_y.set_title("$S_y$")
    f.colorbar(im, ax=ax_y, fraction=0.046, pad=0.04)

Reconstructing these slopes with the Fourier Transform Reconstructor is quite easy::

    >>> from FTR import FourierTransformReconstructor as FTRecon
    >>> recon = FTRecon(np.ones((20,20)), filter='mod_hud', suppress_tt=True)
    >>> phi = recon(sx, sy)


.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from FTR import FourierTransformReconstructor as FTRecon
    recon = FTRecon(np.ones((20,20)), filter='mod_hud', suppress_tt=True)

    sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2

    phi = recon(sx, sy)

    plt.imshow(phi, cmap='hot')
    plt.title(r"$\hat{\phi}$")
    plt.colorbar()


Reconstruction without slope management
=======================================

First, let's set up a circular aperture for our telescope. We'll skip any central obscuration, though in general this shouldn't provide a problem, as the secondary obscuration is continuous across the subapertures.::

    >>> from FTR.utils import circle_aperture
    >>> ap = circle_aperture((20, 20), 8)


.. plot::

    import matplotlib.pyplot as plt
    from FTR.utils import circle_aperture

    ap = circle_aperture((20, 20), 8)
    plt.imshow(ap, cmap='binary_r')
    plt.title("Aperture Mask")


We can naievely use the pure Fourier Transform Reconstructor on apertured slopes::

    >>> recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
    >>> phi_ap = recon(sx * ap, sy * ap)


.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from FTR import FourierTransformReconstructor as FTRecon

    sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2

    recon = FTRecon(np.ones((20,20)), filter='mod_hud', suppress_tt=True)
    phi = recon(sx, sy)

    from FTR.utils import circle_aperture, remove_piston
    ap = circle_aperture((20, 20), 8)
    
    recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
    phi_ap = recon(sx * ap, sy * ap)

    f, (ax_p, ax_r) = plt.subplots(1, 2, figsize=(10,5))
    im = ax_p.imshow(phi_ap * ap, cmap='hot')
    ax_p.set_title(r"$\hat{\phi_{\textrm{ap}}}$")
    f.colorbar(im, ax=ax_p, fraction=0.046, pad=0.04)
    phi_ap_r = remove_piston(ap, phi_ap - phi)[0]*ap
    dmax = np.max(np.abs(phi_ap_r))
    im = ax_r.imshow(phi_ap_r, cmap='bwr', vmin=-dmax, vmax=dmax)
    ax_r.set_title(r"$\hat{\phi_{\textrm{ap}}} - \hat{\phi}$")
    f.colorbar(im, ax=ax_r, fraction=0.046, pad=0.04)

This leaves a characteristic residual pattern around the edge of the aperture.

Reconstruction with slope management
====================================

To improve the reconstruction with an aperture, we apply additional slopes outside the aperture which correct for the periodicity of the system.

Using a the slope-managed reconstructor, we get a much better reconstruction::
    
    >>> from FTR.slopemanage import SlopeManagedFTR
    >>> recon = SlopeManagedFTR(ap, filter='mod_hud', suppress_tt=True)
    >>> phi_sm = recon(sx, sy)
    

.. plot:: 
    
    import numpy as np
    import matplotlib.pyplot as plt
    from FTR import FourierTransformReconstructor as FTRecon
    from FTR.slopemanage import SlopeManagedFTR
    
    sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2
    
    from FTR.utils import circle_aperture, remove_piston, remove_tiptilt
    ap = circle_aperture((20, 20), 8)
    
    recon = SlopeManagedFTR(ap, filter='mod_hud', suppress_tt=True)
    phi_sm = recon(sx * ap, sy * ap)
    
    f, ax_p = plt.subplots(1, 1)
    im = ax_p.imshow(phi_sm * ap, cmap='hot')
    ax_p.set_title(r"$\hat{\phi_{\textrm{sm}}}$")
    f.colorbar(im, ax=ax_p, fraction=0.046, pad=0.04)
    

The residuals from this reconstruction are much improved.

.. plot::
    
    import numpy as np
    import matplotlib.pyplot as plt
    from FTR import FourierTransformReconstructor as FTRecon
    from FTR.slopemanage import SlopeManagedFTR
    
    sx = -np.sin(np.linspace(-10,20,20)[None,:]) * np.ones((20,1)) + 2
    sy = (np.linspace(-1,1,20)[:,None]) * np.ones((1,20)) + 2
    
    recon = FTRecon(np.ones((20,20)), filter='mod_hud', suppress_tt=True)
    phi = recon(sx, sy)
    
    from FTR.utils import circle_aperture, remove_piston, remove_tiptilt
    ap = circle_aperture((20, 20), 8)
    
    recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
    phi_ap = recon(sx * ap, sy * ap)
    
    recon = SlopeManagedFTR(ap, filter='mod_hud', suppress_tt=True)
    phi_sm = recon(sx * ap, sy * ap)
    
    f, (ax_r, ax_ar, ax_fr) = plt.subplots(1, 3)

    phi_sm_fa = remove_tiptilt(ap, remove_piston(ap, phi_sm - phi)[0])[0] * ap
    phi_sm_ap = remove_tiptilt(ap, remove_piston(ap, phi_sm - phi_ap)[0])[0] * ap
    phi_ap_fa = remove_tiptilt(ap, remove_piston(ap, phi_ap - phi)[0])[0] * ap
    dmax = np.max(np.abs([phi_sm_fa, phi_sm_ap, phi_ap_fa]))

    im = ax_r.imshow(phi_sm_fa, cmap='bwr', vmin=-dmax, vmax=dmax)
    ax_r.set_title(r"$\hat{\phi_{\textrm{sm}}} - \hat{\phi}$")

    im = ax_ar.imshow(phi_sm_ap, cmap='bwr', vmin=-dmax, vmax=dmax)
    ax_ar.set_title(r"$\hat{\phi_{\textrm{sm}}} - \hat{\phi_\textrm{ap}}$")

    im = ax_fr.imshow(phi_ap_fa, cmap='bwr', vmin=-dmax, vmax=dmax)
    ax_fr.set_title(r"$\hat{\phi_{\textrm{ap}}} - \hat{\phi}$")
    f.colorbar(im, ax=[ax_r, ax_ar, ax_fr], fraction=0.046/3.0, pad=0.04)