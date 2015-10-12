# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor from FTRLib bound via cython.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .. import ftr
from ..ftr import mod_hud_filter
from ..utils import circle_aperture, remove_piston, complexmp


class TestFTRLibCython(ReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = ftr.FastFTReconstructor
    
    repr = "<CFTRBase {shape}>"
    
    @pytest.fixture
    def shape(self):
        """Return the shape."""
        return (10,10)
    
    @pytest.fixture
    def ap(self, shape):
        """Aperture."""
        return np.ones(shape, dtype=np.bool)
    
    @pytest.fixture
    def phase(self, shape):
        """Random phase."""
        ny, nx = shape
        # phi_ft = complexmp(10.0, 2 * np.pi * np.random.random(shape))
        phi_ft = complexmp(0.0, np.zeros(shape)) * 0.0
        phi_ft[3,3] = 2.0
        phi_ft[7,7] = 2.0
        
        # We zero everything beyond nyquist, so that we can make it hermitian.
        # phi_ft[ny//2:,:] = 0.0
        # phi_ft[:,nx//2:] = 0.0
        
        phi = np.real(np.fft.ifftn(phi_ft))
        # phi = np.ones(shape)
        # phi[::3,::3] = 2.0
        return phi
    
    @pytest.fixture
    def reconstructor(self, ap):
        """Construct the reconstructor."""
        recon = self.cls(ap, filter='mod_hud')
        gs = np.ones(ap.shape, dtype=np.complex)
        recon.gx = gs
        recon.gy = gs
        return recon
        
    def test_init_filter(self, ap):
        """Test initialization with different filters."""
        reconstructor = self.cls(ap)
        assert reconstructor.gx.shape == ap.shape
        assert reconstructor.gy.shape == ap.shape
        
    def test_filter_roundtrip(self, reconstructor, phase, shape):
        """Test filter roundtripping."""
        xs, ys = reconstructor.invert(phase)
        phase_rt = reconstructor.reconstruct(xs, ys)
        assert np.allclose(phase_rt, phase)