#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Try the Cython-based FTR.
"""

import numpy as np
from FTR.ftr import FastFTReconstructor

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
        # p_ft[0,-1] = 0.0
        # p_ft[-1,0] = 0.0
        # p_ft[-1,-1] = 0.0
        phase = np.real(np.fft.ifftn(p_ft))
    else:
        index = (5,3)
        phase_ft = np.zeros(opt.shape, dtype=np.complex)
        phase_ft[index] = 5.0
        phase = np.real(np.fft.ifftn(phase_ft))
    
    ap = np.ones(opt.shape, dtype=np.int)
    
    print("Initializing CFTRBase...")
    recon = FastFTReconstructor(ap, filter=opt.filter)
    print("Checking Attribute Access...")
    assert recon.name == opt.filter
    assert recon.shape == opt.shape
    assert recon.gx.shape == opt.shape
    assert recon.gy.shape == opt.shape
    assert recon.denominator.shape == opt.shape
    
    
    
    print("Calling reconstructor {} times...".format(opt.iter))
    success = True
    for _ in range(opt.iter):
        xs, ys = recon.invert(phase)
        phi = recon(xs, ys)
        success = success and np.allclose(phi, phase)
    if success:
        print("Success!")
    else:
        np.set_printoptions(precision=4, linewidth=120, suppress=True)
        x, y = np.fft.fftn(xs)[index], np.fft.fftn(ys)[index]
        est = (x * np.conj(recon.gx[index]) + y * np.conj(recon.gy[index]))/recon.denominator[index]
        _x,_y = index
        print("py: sx[{x:d},{y:d}]={:.4f} sy[{x:d},{y:d}]={:.4f} es[{x:d},{y:d}]={:.4f}".format(x,y,est,x=_x,y=_y))
        print("py: gx[{x:d},{y:d}]={:.4f} gy[{x:d},{y:d}]={:.4f} gd[{x:d},{y:d}]={:.4f}".format(recon.gx[index], recon.gy[index], recon.denominator[index],x=_x,y=_y))
        print("Finished, but round-tripping failed.")
        
        
        import matplotlib.pyplot as plt
        fig, (ax_1, ax_2, ax_r) = plt.subplots(1,3)
        for ax, data, label in zip([ax_1, ax_2], [phase, phi], ["Original", "Estimate"]):
            im = ax.imshow(np.real(np.fft.fftn(data)), cmap='hot')
            ax.set_title(label)
            fig.colorbar(im, ax=ax)
        print(phase/phi)
        nf = opt.shape[-1] // 2 + 1
        print(nf,opt.shape[-1],nf/opt.shape[-1])
        # fig, ((a_x1, a_x2),(a_y1,a_y2)) = plt.subplots(2,2)
        # for ax, data, label in zip((a_x1, a_x2,a_y1,a_y2),(xs,sx,ys,sy),"xs,sx,ys,sy".split(",")):
        #     ax.imshow(np.real(np.fft.fftn(data)), cmap='hot')
        #     ax.set_title(label)
        plt.show()
        
        
if __name__ == '__main__':
    main()