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
ax_r.set_title(r"$\hat{\phi_{sm}} - \hat{\phi}$")

im = ax_ar.imshow(phi_sm_ap, cmap='bwr', vmin=-dmax, vmax=dmax)
ax_ar.set_title(r"$\hat{\phi_{sm}} - \hat{\phi_{ap}}$")

im = ax_fr.imshow(phi_ap_fa, cmap='bwr', vmin=-dmax, vmax=dmax)
ax_fr.set_title(r"$\hat{\phi_{ap}} - \hat{\phi}$")
f.colorbar(im, ax=[ax_r, ax_ar, ax_fr], fraction=0.046/3.0, pad=0.04)