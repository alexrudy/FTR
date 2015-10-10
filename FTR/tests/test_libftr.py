# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor from FTRLib bound via cython.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .. import _ftr
from ..ftr import mod_hud_filter
from ..utils import circle_aperture, remove_piston, complexmp


class TestFTRLibCython(ReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = _ftr.CFTRBase
    
    repr = "<CFTRBase {shape}>"
        
    @pytest.fixture
    def phase(self, shape):
        """Random phase."""
        ny, nx = shape
        phi_ft = complexmp(10.0, 2 * np.pi * np.random.random(shape))
        phi_ft[0,0] = 0.0
        
        # We zero everything beyond nyquist, so that we can make it hermitian.
        phi_ft[ny//2:,:] = 0.0
        phi_ft[:,nx//2:] = 0.0
        
        phi = np.real(np.fft.ifftn(phi_ft))
        return phi
    
    @pytest.fixture
    def reconstructor(self, shape):
        """Construct the reconstructor."""
        recon = self.cls(shape)
        filter = mod_hud_filter(shape)
        recon.gx = filter.gx
        recon.gy = filter.gy
        return recon
        
    def test_init_filter(self, shape):
        """Test initialization with different filters."""
        reconstructor = self.cls(shape)
        assert reconstructor.gx.shape == shape
        assert reconstructor.gy.shape == shape
        
    # def test_filter_roundtrip(self, reconstructor, phase, shape):
    #     """Test filter roundtripping."""
    #     xs, ys = reconstructor.invert(phase)
    #     phase_rt = reconstructor.reconstruct(xs, ys)
    #     assert np.allclose(phase_rt, phase)