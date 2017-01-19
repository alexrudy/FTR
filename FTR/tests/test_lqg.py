# -*- coding: utf-8 -*-
"""
Tests for the LQG module.
"""

from __future__ import absolute_import
import pytest
import numpy as np

from .test_filter_base import FilterTestBase
from ..lqg import LQGFilter, FastLQGFilter

class TestLQGFilter(FilterTestBase):
    """Tests for the LQG filter."""
    
    cls = LQGFilter
    
    repr = "<LQGFilter nlayers=... shape={shape}>"
    
    @pytest.fixture(params=[1,3,8])
    def nlayers(self, request):
        """Return the number of layers."""
        return request.param
        
    @pytest.fixture
    def full_filter_shape(self, nlayers, shape):
        """Shape of the filter, accounting for number of layers."""
        return (nlayers,) + shape
        
    @pytest.fixture
    def gains(self, full_filter_shape):
        """Gains"""
        return np.random.randn(*full_filter_shape) + 1j * np.random.randn(*full_filter_shape)
        
    @pytest.fixture
    def alphas(self, full_filter_shape):
        """Alphas"""
        return np.random.randn(*full_filter_shape) + 1j * np.random.randn(*full_filter_shape)
        
        
    @pytest.fixture
    def hp_coeffs(self, shape):
        """High-pass coefficients."""
        return np.random.randn(*shape) + 1j * np.random.randn(*shape)
    
    @pytest.fixture
    def filter(self, gains, alphas, hp_coeffs):
        """Generate a random filter for testing."""
        return self.cls(gains, alphas, hp_coeffs)
        
    

class TestCLQG(TestLQGFilter):
    """Test the LQG filter."""
    
    cls = FastLQGFilter
    
    repr = "<FastLQGFilter nlayers=... shape={shape}>"
    
    @pytest.fixture
    def filter_py(self, gains, alphas, hp_coeffs):
        """Pure-python LQG filter"""
        return LQGFilter(gains, alphas, hp_coeffs)
    
    def test_compare_to_python(self, filter, filter_py, phase_ft):
        """docstring for test_compare_to_python"""
        np.testing.assert_allclose(filter.gains, filter_py.gains)
        np.testing.assert_allclose(filter.alphas, filter_py.alphas)
        np.testing.assert_allclose(filter.highpass_coefficients, filter_py.highpass_coefficients)
        
        filter.reset()
        phase_c1 = filter(phase_ft)
        phase_c2 = filter(phase_ft)
        
        filter_py.reset()
        phase_py1 = filter_py(phase_ft)
        phase_py2 = filter_py(phase_ft)
        
        np.testing.assert_allclose(phase_c1, phase_py1, atol=1e-7)
        np.testing.assert_allclose(phase_c2, phase_py2, atol=1e-7)
    