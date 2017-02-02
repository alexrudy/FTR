# -*- coding: utf-8 -*-
"""
Tests for the Abstract implementation details of post-reconstruction filters.
"""
import pytest
import random
import six
import numpy as np

from ..utils import shapestr, remove_piston

def _compare_attribute_numpy(a, b, attr):
    """Compare a single attribute"""


def _compare_attributes_numpy(a, b):
    """Compare numpy array-like attributes."""
    dirb = set(dir(b))
    
    for attr in dir(a):
        attra = getattr(a, attr)
        if not isinstance(attra, np.ndarray):
            continue
        assert attr in dir(b), "b={0!r} missing attribute {1:s}".format(b, attr)
        attrb = getattr(b, attr)
        np.testing.assert_allclose(attra, attrb)
        dirb.remove(attr)
    
    for attr in dirb:
        attrb = getattr(b, attr)
        if not isinstance(attrb, np.ndarray):
            continue
        assert attr in dir(a), "a={0!r} missing attribute {1:s}".format(a, attr)
        attra = getattr(a, attr)
        np.testing.assert_allclose(attra, attrb)
    

class FilterTestBase(object):
    """A base class for FTR filter tests."""
    
    cls = None
    
    repr = None
        
    @pytest.fixture
    def pzero(self, shape):
        """Zero slopes"""
        est_ft = np.zeros(shape, dtype=np.complex)
        return est_ft
    
    @pytest.fixture
    def phase_ft(self, shape):
        """Random phase."""
        return np.random.randn(*shape) + 1j * np.random.randn(*shape)
    
    @pytest.fixture
    def filter(self):
        """Create the reconstructor"""
        return self.cls()
        
    def pytest_generate_tests(self, metafunc):
        if 'format' in metafunc.fixturenames:
            formats = [ str(fmt) for fmt, read_write in self.cls.formats().items() if read_write == (True, True) ]
            metafunc.parametrize("format", formats)
    
    def test_filter_repr(self, filter, shape):
        """Test the repr method."""
        r = repr(filter)
        assert isinstance(r, six.string_types)
        if self.repr is not None:
            expected = self.repr.format(shape=shapestr(shape))
            if "..." in expected:
                start, end = expected.split("...", 1)
                assert r.startswith(start)
                assert r.endswith(end)
            else:
                assert r == expected
    
    def test_filter_zeros(self, filter, pzero):
        """Reconstruction with zeros should return all zeros."""
        est_ft = filter(pzero)
        
    def test_filter_roundtrip(self, filter, tmpdir, format):
        """Test a read-write filter roundtrip."""
        filename = str(tmpdir.join('test.{0}'.format(format)))
        filter.write(filename)
        new_filter = self.cls.read(filename)
        _compare_attributes_numpy(filter, new_filter)
        
def test_filter_abc():
    """Test filter is abstract."""
    from ..base import Filter
    with pytest.raises(TypeError):
        Filter()
    
from ..filters import StaticFilter
    
class TestStaticFilters(FilterTestBase):
    """Test static filters."""
    
    cls = StaticFilter
    
    repr = "<StaticFilter {shape} ..."
    
    @pytest.fixture
    def filter(self, shape):
        """Generate a random filter for testing."""
        return self.cls(np.random.randn(*shape) + 1j * np.random.randn(*shape))
    