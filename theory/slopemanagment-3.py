import matplotlib.pyplot as plt
from FTR.utils import circle_aperture

ap = circle_aperture((20, 20), 8)
plt.imshow(ap, cmap='binary_r')
plt.title("Aperture Mask")