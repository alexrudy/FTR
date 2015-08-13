# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .. import ftr
from ..utils import circle_aperture, remove_piston, complexmp

@pytest.fixture(params=ftr.FourierTransformReconstructor.filters())
def filter_name(request):
    """An FTR filter name."""
    return request.param

class TestFourierTransformReconstructor(ReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = ftr.FourierTransformReconstructor
    
    repr = "<FourierTransformReconstructor {shape} filter='mod_hud'>"
    
    @pytest.fixture
    def ap(self, shape):
        """Aperture."""
        return np.ones(shape, dtype=np.bool)
        # r = np.min([np.floor((s / 2.0) - 2) for s in shape])
        # return circle_aperture(shape, r)
        
    @pytest.fixture
    def phase(self, shape):
        """Random phase."""
        ny, nx = shape
        phi_ft = complexmp(10.0, 2 * np.pi * np.random.random(shape))
        phi_ft[0,0] = 0.0
        phi_ft[ny//2:,:] = 0.0
        phi_ft[:,nx//2:] = 0.0
        
        phi = np.real(np.fft.ifftn(phi_ft))
        return phi
    
    @pytest.fixture
    def reconstructor(self, ap):
        """Construct the reconstructor."""
        return self.cls(ap, filter="mod_hud")
        
    def test_init_filter(self, shape, ap, filter_name):
        """Test initialization with different filters."""
        reconstructor = self.cls(ap, filter=filter_name)
        assert reconstructor.name == filter_name
        assert reconstructor.gx.shape == shape
        assert reconstructor.gy.shape == shape
        assert reconstructor.denominator.shape == shape
        
    def test_filter_roundtrip(self, ap, filter_name, phase, shape):
        """Test filter roundtripping."""
        reconstructor = self.cls(ap, filter=filter_name)
        xs, ys = reconstructor.invert(phase)
        phase_rt = reconstructor.reconstruct(xs, ys)
        assert np.allclose(phase_rt, phase)
        