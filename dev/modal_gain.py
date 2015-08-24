#!/usr/bin/env python
# -*- coding: utf-8 -*-
from FTR import FourierTransformReconstructor
from FTR.slopemanage import SlopeManagedFTR
from FTR.utils import circle_aperture, shapegrid, ignoredivide, complexmp, make_hermitian
import numpy as np
import numpy.fft
import numpy.random
import argparse
import six

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def rms_with_ap(data, ap):
    """RMS with an aperture"""
    return np.sqrt(np.sum((data * ap) ** 2.0) / np.sum(ap))

def fourier_modal_gain(recon, fx=slice(None, None), fy=slice(None, None)):
    """Collect the Fourier modal gains for a reconstructor."""
    FY, FX = shapegrid(recon.shape, centered=False)
    
    gain = np.zeros(recon.shape)
    rms = np.zeros(recon.shape)
    
    for _fy, _fx in zip(FY[fy, fx], FX[fy, fx]):
        grid = np.zeros(recon.shape, dtype=np.complex)
        grid[_fy,_fx] = 1.0
        # complexmp(np.random.uniform(0.0,10.0), np.random.uniform(0.0, 2*np.pi))
        phi = np.real(np.fft.ifftn(grid))
        sx, sy = recon.invert(phi)
        phi_r = recon(sx * recon.ap, sy * recon.ap)
        phi_ft = np.fft.fftn(phi_r)
        gain[_fy, _fx] = np.abs(phi_ft[_fy, _fx]) / np.abs(grid[_fy, _fx])
        rms[_fy, _fx] = rms_with_ap((phi_r - phi), recon.ap)
    return (gain, rms)

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
    
    ap = circle_aperture((opt.n, opt.n), r=(opt.n//2) - 1)
    
    sm_recon = SlopeManagedFTR(ap, opt.filter, manage_tt=opt.manage_tt, suppress_tt=False)
    recon = FourierTransformReconstructor(ap, opt.filter,
        manage_tt=opt.manage_tt, suppress_tt=False)
    
    fy = 0
    fig, (ax_gain, ax_rms) = plt.subplots(2, 1)
    ax_gain.set_ylabel("Gain")
    ax_rms.set_ylabel("RMS")
    ax_rms.set_xlabel("$f_x$")
    for fy in range(0, 1):
        for r, color in zip([recon, sm_recon], "bg"):
            gain, rms = fourier_modal_gain(r, fy=fy)
            ax_gain.plot(gain[fy,:], color=color)
            ax_rms.plot(rms[fy,:], color=color)
        
    
    # fig, (ax_gain, ax_rms) = plt.subplots(1, 2)
    
    
    plt.show()
    
if __name__ == '__main__':
    main()