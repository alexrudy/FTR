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

class FourierTransformReconstructorTestBase(ReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = ftr.FourierTransformReconstructor
    
    repr = "<FourierTransformReconstructor {shape} filter='mod_hud'>"
    
    @pytest.fixture
    def ap(self, shape):
        """Aperture."""
        return np.ones(shape, dtype=np.bool)
        
    @pytest.fixture
    def phase(self, shape):
        """Random phase."""
        ny, nx = shape
        phi_ft = np.zeros(shape, dtype=np.complex)
        phi_ft[0,0] = 0.0
        phi_ft[3,3] = 1.0
        phi = np.real(np.fft.ifftn(phi_ft))
        return phi
    
    @pytest.fixture(scope='function')
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
        np.testing.assert_allclose(phase_rt, phase, atol=1e-7)
        
class TestFourierTransformReconstructor(FourierTransformReconstructorTestBase):
    pass

