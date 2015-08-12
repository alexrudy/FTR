# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .. import ftr
from ..utils import circle_aperture

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
        r = np.min([np.floor((s / 2.0) - 2) for s in shape])
        return circle_aperture(shape, r)
    
    @pytest.fixture
    def reconstructor(self, shape, ap):
        """Construct the reconstructor."""
        return self.cls(shape, ap, filter="mod_hud")
        
    def test_init_filter(self, shape, ap, filter_name):
        """Test initialization with different filters."""
        reconstructor = self.cls(shape, ap, filter=filter_name)
        assert reconstructor.name == filter_name
        assert reconstructor.gx.shape == shape
        assert reconstructor.gy.shape == shape
        assert reconstructor.denominator.shape == shape
        