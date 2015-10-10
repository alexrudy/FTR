#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Try the Cython-based FTR.
"""

import numpy as np
from FTR._ftr import CFTRBase

if __name__ == '__main__':
    shape = (20,20)
    
    # sx = np.random.randn(*shape)
    # sy = np.random.randn(*shape)
    sx = np.zeros(shape, dtype=np.float)
    sy = np.zeros(shape, dtype=np.float)
    
    print("Initializing CFTRBase...")
    recon = CFTRBase(shape)
    recon.gx = np.zeros(shape, dtype=np.complex)
    recon.gy = np.zeros(shape, dtype=np.complex)
    print("Calling reconstructor...")
    phi = recon(sx, sy)
    print("Success!!!")
    