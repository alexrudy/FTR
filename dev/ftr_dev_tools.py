# -*- coding: utf-8 -*-
#
#  ftr_dev_tools.py
#  FTR
#
#  Created by Alexander Rudy on 2015-08-12.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#
import numpy as np
import numpy.fft

import matplotlib.pyplot as plt
import matplotlib.gridspec as mgrid

from FTR.utils import ignoredivide

def fftmag(a, piston=False, waffle=False):
    """FFT magnitude."""
    a_ft = np.real(np.abs(np.fft.fftn(np.asarray(a))))
    if not piston:
        a_ft[0,0] = 0.0
    a_ft = np.fft.fftshift(a_ft)
    if not waffle:
        a_ft[0,0] = 0.0
    return a_ft

def filtermag(a):
    """Filter FFT magnitude."""
    return np.real(np.abs(np.fft.fftshift(np.asarray(a))))

def reconstruct_and_display(reconstructor, xs, ys, phi, rtol = 0.001, atol = 0.025, **kwargs):
    """Reconstruct and display."""
    
    
    phi_r = reconstructor(xs, ys)
    display_reconstructor(phi, phi_r, **kwargs)
    
    with ignoredivide():
        
        if np.allclose(phi_r, phi, rtol=rtol, atol=atol):
            print("Success! Filter and np.gradient match each other.")
        else:
            ap = reconstructor.ap
            m_ok = (np.abs(phi_r - phi) * ap) <= atol + np.abs(phi * ap) * rtol
            p_ok = (m_ok).sum() / np.prod(phi.shape)
        
            d = (np.abs(phi_r - phi) * ap)
            d_r = (np.abs(phi_r - phi) * ap) / np.abs(phi * ap)
            d_m = (np.abs(phi_r) / np.abs(phi))
        
            print("{0:.2%} of points are O.K.".format(p_ok))
            print("rtol={0:.2g}".format(d_r[ap.astype(np.bool)].mean()))
            print("atol={0:.2g}".format(d[ap.astype(np.bool)].mean()))
            print("mtol={0:.2g}".format(np.nanmean(d_m[ap.astype(np.bool)])))
    
    return phi_r
    
    
def display_reconstructor(phi, phi_r, figure=None, cmap='hot'):
    """Reconstructor display."""
    if figure is None:
        figure = plt.figure()
    
    figure.suptitle("Phase")
    gs = mgrid.GridSpec(2, 3)
    
    ax = figure.add_subplot(gs[0,0])
    im = ax.imshow(phi, cmap=cmap)
    ax.set_title(r"$\phi$")
    figure.colorbar(im)
    
    ax = figure.add_subplot(gs[1,0])
    im = ax.imshow(fftmag(phi), cmap=cmap)
    ax.set_title(r"$F(\phi)$")
    figure.colorbar(im)
    
    ax = figure.add_subplot(gs[0,1])
    im = ax.imshow(phi_r, cmap=cmap)
    ax.set_title(r"$\hat{\phi}$")
    figure.colorbar(im)
    
    ax = figure.add_subplot(gs[1,1])
    im = ax.imshow(fftmag(phi_r), cmap=cmap)
    ax.set_title(r"$F(\hat{\phi})$")
    figure.colorbar(im)
    
    ax = figure.add_subplot(gs[0,2])
    ax.set_title(r"$\phi - \hat{\phi}$")
    im = ax.imshow(np.abs(phi - phi_r), cmap='hot')
    figure.colorbar(im)
    
    ax = figure.add_subplot(gs[1,2])
    ax.set_title(r"$F(\phi) - F(\hat{\phi})$")
    im = ax.imshow((fftmag(phi) - fftmag(phi_r)), cmap='hot')
    figure.colorbar(im)
    
    
    return figure