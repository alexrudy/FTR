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