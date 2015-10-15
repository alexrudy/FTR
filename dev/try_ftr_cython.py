#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Try the Cython-based FTR.
"""

import numpy as np
from FTR.ftr import FastFTReconstructor, FourierTransformReconstructor

def index_view(index, recon, xs, ys):
    """Produce a printed view of a specific index."""
    x, y = np.fft.fftn(xs)[index], np.fft.fftn(ys)[index]
    gd = 1.0 /recon.denominator[index] * 2.0/np.prod(recon.shape)
    est = (x * np.conj(recon.gx[index]) + y * np.conj(recon.gy[index])) * gd
    _x,_y = index
    print("py: sx[{x:d},{y:d}]={0:.4f} sy[{x:d},{y:d}]={1:.4f} es[{x:d},{y:d}]={2:.4f}".format(x,y,est,x=_x,y=_y))
    print("py: gx[{x:d},{y:d}]={:.4f} gy[{x:d},{y:d}]={:.4f} gd[{x:d},{y:d}]={:.4f}".format(recon.gx[index], recon.gy[index], gd, x=_x,y=_y))

def main():
    """Main function."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("shape", type=int, nargs=2, help="Shape to pack/unpack.")
    parser.add_argument('-r','--random', action='store_true', help="Use random numbers.")
    parser.add_argument('-f','--filter', type=str, choices=FastFTReconstructor.filters(), default='hud', help="Filter name.")
    parser.add_argument('-i','--iter', type=int, default=1, help="Number of iterations.")
    opt = parser.parse_args()
    
    opt.shape = tuple(opt.shape)
    if opt.random:
        phase = np.random.randn(*opt.shape)
        p_ft = np.fft.fftn(phase)
        p_ft[0,0] = 0.0
        p_ft[5,5] = 0.0
        index = (3,3)
        phase = np.real(np.fft.ifftn(p_ft))
    else:
        index = (4,5)
        phase_ft = np.zeros(opt.shape, dtype=np.complex)
        phase_ft[index] = 5.0
        phase_ft[2,4] = 2.0
        phase_ft[1,6] = 2.0
        phase = np.real(np.fft.ifftn(phase_ft))
    
    ap = np.ones(opt.shape, dtype=np.int)
    
    print("Initializing CFTRBase...")
    recon = FastFTReconstructor(ap, filter=opt.filter)
    recon2 = FourierTransformReconstructor(ap, filter=opt.filter)
    print("Checking Attribute Access...")
    assert recon.name == opt.filter
    assert recon.shape == opt.shape
    assert recon.gx.shape == opt.shape
    
    assert recon.gy.shape == opt.shape
    assert recon.denominator.shape == opt.shape
    
    # Check against second reconstructor.
    assert recon2.name == recon.name
    assert recon2.shape == recon.shape
    assert recon2.gx.shape == recon.gx.shape
    assert recon2.gy.shape == recon.gy.shape
    assert recon2.denominator.shape == recon.denominator.shape
    np.testing.assert_allclose(recon2.gx, recon.gx)
    np.testing.assert_allclose(recon2.gy, recon.gy)
    np.testing.assert_allclose(recon2.denominator, recon.denominator)
    
    print("Calling reconstructor {} times...".format(opt.iter))
    success = True
    for _ in range(opt.iter):
        x2, y2 = recon2.invert(phase)
        xs, ys = recon.invert(phase)
        phi = recon(xs, ys)
        ph2 = recon2(xs, ys)
        success = success and np.allclose(phi, phase)
    if success:
        print("Success!")
    else:
        np.set_printoptions(precision=4, linewidth=120, suppress=True)
        index_view(index, recon, xs, ys)
        index_view((2,4), recon, xs, ys)
        index_view((1,6), recon, xs, ys)
        # index_view((9,4), recon, xs, ys)
        print("Finished, but round-tripping failed.")
        
        
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1,3)
        for ax, data, label in zip(axes, [phase, phi, ph2], ["Original", "Estimate", "Python"]):
            im = ax.imshow(np.imag(np.fft.fftn(data)), cmap='hot')
            ax.set_title(label)
            fig.colorbar(im, ax=ax)
        
        plt.pause(0.1)
        print(np.fft.fftn(phase))
        print(np.fft.fftn(phi))
        try:
            np.testing.assert_allclose(xs, x2)
            np.testing.assert_allclose(ys, y2)
            np.testing.assert_allclose(phase, ph2)
            np.testing.assert_allclose(phi, ph2)
        except AssertionError:
            plt.show()
            raise
        
        
if __name__ == '__main__':
    main()