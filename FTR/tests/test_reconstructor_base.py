# -*- coding: utf-8 -*-
"""
Tests for the Abstract implementation details of the reconstructor.
"""
import pytest
import random
import six
import numpy as np

from ..utils import shapestr, remove_piston

class ReconstructorTestBase(object):
    """A base class for reconstructor tests."""
    
    cls = None
    
    repr = None
    
    @pytest.fixture
    def szero(self, shape):
        """Zero slopes"""
        x = np.zeros(shape, dtype=np.float)
        y = x.copy()
        return (x, y)
        
    @pytest.fixture
    def srandom(self, shape):
        """Random slopes"""
        x = np.random.randn(*shape)
        y = np.random.randn(*shape)
        return (x, y)
    
    @pytest.fixture
    def phase(self, shape):
        """Random phase."""
        return np.random.randn(*shape)
    
    @pytest.fixture
    def reconstructor(self, shape):
        """Create the reconstructor"""
        return self.cls(shape)
    
    def test_reconstructor_shape(self, reconstructor, shape):
        """Test reconstructor shape"""
        assert reconstructor.shape == shape
    
    def test_reconstructor_repr(self, reconstructor, shape):
        """Test the repr method."""
        r = repr(reconstructor)
        assert isinstance(r, six.string_types)
        if self.repr is not None:
            assert r == self.repr.format(shape=shapestr(shape))
    
    def test_reconstructor_zeros(self, reconstructor, szero):
        """Reconstruction with zeros should return all zeros."""
        sx, sy = szero
        phi = reconstructor(sx, sy)
        assert np.allclose(phi, 0.0)
        
    @pytest.mark.skipif("sys.version_info < (2,7)")
    def test_recon_bench(self, benchmark, reconstructor, srandom):
        """Benchmark reconstructor results."""
        sx, sy = srandom
        benchmark.name = benchmark.name.replace("test_recon_bench", "recon") + "[{0}]".format(self.name)
        phi = benchmark(reconstructor, sx, sy)
        
    