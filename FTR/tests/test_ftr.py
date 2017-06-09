# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .test_io import IOTestBase
from .. import ftr
from ..utils import circle_aperture, remove_piston, complexmp

@pytest.fixture(params=ftr.FourierTransformReconstructor.filters())
def filter_name(request):
    """An FTR filter name."""
    return request.param

class TestFTRFilter(IOTestBase):
    """Test FTR Filter."""
    
    cls = ftr.FTRFilter
    
    @pytest.fixture
    def obj(self, filter_name, shape):
        """FTR filter object."""
        return ftr.FourierTransformReconstructor.get(filter_name, shape)
        

class FourierTransformReconstructorTestBase(ReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = ftr.FourierTransformReconstructor
    name = "FTR"
    
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
        
    def test_reconstructor_tt_flags_init(self, ap, filter_name):
        """Test tip-tilt control flags as __init__ args."""
        reconstructor = self.cls(ap, filter=filter_name)
        assert reconstructor.manage_tt is False
        assert reconstructor.suppress_tt is False
        
        reconstructor = self.cls(ap, filter=filter_name, manage_tt=True)
        assert reconstructor.manage_tt is True
        assert reconstructor.suppress_tt is False
        
        reconstructor = self.cls(ap, filter=filter_name, suppress_tt=True)
        assert reconstructor.manage_tt is True
        assert reconstructor.suppress_tt is True
        
        with pytest.raises(ValueError):
            reconstructor = self.cls(ap, filter=filter_name, suppress_tt=True, manage_tt=False)
        
    def test_reconstructor_tt_flags(self, ap, filter_name):
        """Test tip-tilt control flags as attributes."""
        reconstructor = self.cls(ap, filter=filter_name)
        assert reconstructor.manage_tt is False
        assert reconstructor.suppress_tt is False
    
        reconstructor = self.cls(ap, filter=filter_name)
        reconstructor.manage_tt = True
        assert reconstructor.manage_tt is True
        assert reconstructor.suppress_tt is False
    
        reconstructor = self.cls(ap, filter=filter_name)
        reconstructor.suppress_tt = True
        assert reconstructor.manage_tt is True
        assert reconstructor.suppress_tt is True
        
        reconstructor = self.cls(ap, filter=filter_name)
        reconstructor.manage_tt = False
        with pytest.raises(ValueError):
                reconstructor.suppress_tt = True
                
    def test_get_filter_attribute(self, ap, filter_name):
        """Test filter attribute."""
        reconstructor = self.cls(ap, filter=filter_name)
        filter_data = reconstructor.filter
        assert isinstance(filter_data, ftr.FTRFilter)
        np.testing.assert_allclose(filter_data.gx, reconstructor.gx)
        np.testing.assert_allclose(filter_data.gy, reconstructor.gy)
        assert filter_data.name == filter_name
        gx, gy, name = filter_data
        np.testing.assert_allclose(gx, reconstructor.gx)
        np.testing.assert_allclose(gy, reconstructor.gy)
        assert name == filter_name
            
class TestFourierTransformReconstructor(FourierTransformReconstructorTestBase):
    pass

