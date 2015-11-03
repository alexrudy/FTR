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