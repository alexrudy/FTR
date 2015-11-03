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
ax_p.set_title(r"$\hat{\phi_{sm}}$")
f.colorbar(im, ax=ax_p, fraction=0.046, pad=0.04)