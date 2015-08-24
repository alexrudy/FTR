#!/usr/bin/env python
# -*- coding: utf-8 -*-

from FTR import FourierTransformReconstructor
import FTR.ftr
from FTR.utils import circle_aperture, shapegrid, ignoredivide, complexmp, make_hermitian
import numpy as np
import numpy.fft
import argparse
import six

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from ftr_dev_tools import reconstruct_and_display
    
def make_aperture(n):
    """Make a square aperture with a small border."""
    ap = np.ones((n, n), dtype=np.int)
    # ap[:2,:] = 0.0
    # ap[:,:2] = 0.0
    # ap[-2:,:] = 0.0
    # ap[:,-2:] = 0.0
    return ap
    
def prepare_grids(n):
    """Prepare grids."""
    y, x = shapegrid((n, n))
    ap = make_aperture(n)
    
    # phi = np.random.randn(N, N)
    phi = 2.0 * y **2.0 - 1.5 * x ** 2.0 + (x * -y)
    phi -= np.mean(phi)
    phi *= ap
    
    ys, xs = np.gradient(phi)
    return ap, phi, ys, xs
    
def prepare_fmode_grids(n):
    """Prepare grids in fourier space."""
    ap = make_aperture(n)
    phi_ft = complexmp(10.0, 2 * np.pi * np.random.random((n, n)))
    phi_ft = (phi_ft + np.triu(phi_ft.T)).T
    phi_ft = make_hermitian(phi_ft)
    phi_ft[0,:] = 0.0
    phi_ft[:,0] = 0.0
    phi_ft[-1,:] = 0.0
    phi_ft[:,-1] = 0.0
    plt.figure()
    plt.subplot(121)
    plt.imshow(np.real(phi_ft), cmap='hot')
    phi_ft = np.fft.ifftshift(phi_ft)
    plt.subplot(122)
    plt.imshow(np.real(phi_ft), cmap='hot')
    # This should remove waffle!
    phi_ft[0,0] = 0.0
    phi = np.abs(np.fft.ifftn(phi_ft))
    phi[0,0] = 0.0
    ys, xs = np.gradient(phi)
    return ap, phi, ys, xs

def main():
    """Main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument('filter', default="fried", type=six.text_type,
        choices=FourierTransformReconstructor.filters(), nargs="?",
        help="filter name for FTR reconstructor, default='fried'")
    parser.add_argument("--tt", help="leave tip/tilt in slopes",
         action='store_false', dest='manage_tt')
    parser.add_argument("-n", type=int, default=20,
        help="number of subapertures on a side.")
    opt = parser.parse_args()
    ap, phi, ys, xs = prepare_fmode_grids(opt.n)
    recon = FourierTransformReconstructor(ap, opt.filter,
        manage_tt=opt.manage_tt, suppress_tt=False)
    print(recon)
    xs, ys = recon.invert(phi)
    phi_r = reconstruct_and_display(recon, xs, ys, phi)
        
    plt.show()
    
if __name__ == '__main__':
    main()