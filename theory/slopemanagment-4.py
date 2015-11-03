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
ax_p.set_title(r"$\hat{\phi_{ap}}$")
f.colorbar(im, ax=ax_p, fraction=0.046, pad=0.04)
phi_ap_r = remove_piston(ap, phi_ap - phi)[0]*ap
dmax = np.max(np.abs(phi_ap_r))
im = ax_r.imshow(phi_ap_r, cmap='bwr', vmin=-dmax, vmax=dmax)
ax_r.set_title(r"$\hat{\phi_{ap}} - \hat{\phi}$")
f.colorbar(im, ax=ax_r, fraction=0.046, pad=0.04)