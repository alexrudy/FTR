#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Try the Cython-based FTR.
"""
import time
import numpy as np
from FTR.ftr import FastFTReconstructor, FourierTransformReconstructor

def random_phase(shape):
    """Generate phase."""
    phase = np.random.randn(*shape)
    p_ft = np.fft.fftn(phase)
    p_ft[0,0] = 0.0
    p_ft[5,5] = 0.0
    return np.real(np.fft.ifftn(p_ft))

def main():
    """Main function."""
    global start
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("shape", type=int, nargs=2, help="Shape to pack/unpack.")
    parser.add_argument('-f','--filter', type=str, choices=FastFTReconstructor.filters(), default='hud', help="Filter name.")
    parser.add_argument('-i','--iter', type=lambda s : int(float(s)), default=1000, help="Number of iterations.")
    opt = parser.parse_args()
    
    opt.shape = tuple(opt.shape)
    ap = np.ones(opt.shape, dtype=np.int)
    phase = random_phase(opt.shape)
    
    start = time.time()
    reconstructors = {}
    results = {}
    
    reconstructors['libFTR'] = FastFTReconstructor(ap, filter=opt.filter)
    reconstructors['fftpack'] = FourierTransformReconstructor(ap, filter=opt.filter)
    xs, ys = reconstructors['fftpack'].invert(phase)
    
    for name, recon in reconstructors.items():
        rstart = time.time()
        for _ in range(opt.iter):
            phi = recon(xs, ys)
        rduration = (time.time() - rstart) / opt.iter
        print("Reconstructor {} took {:.0f} ns per iteration = {:.0f} kHz.".format(name, rduration * 1e9, 1e-3 / rduration))
        results[name] = rduration
    
    print("libFTR was {:.1f} times faster than fftpack.".format(results['fftpack'] / results['libFTR']))
    
if __name__ == '__main__':
    main()