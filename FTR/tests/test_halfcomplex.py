# -*- coding: utf-8 -*-

import numpy as np
import pytest
import random
import copy

from ..libftr._ftr import HalfComplexMapping
from ..utils import unpack_halfcomplex


@pytest.fixture(params=[(0, 0), (0, 3), (3, 1)])
def shape(request):
    """Random integer size n."""
    row_offset, col_offset = request.param
    nn = (random.randint(5, 128) * 2)
    return (nn + row_offset, nn + col_offset)

@pytest.fixture
def halfcomplex_data(shape):
    """Random data."""
    shape_ft = list(shape)
    shape_ft[1] = (shape[1] // 2) + 1
    return np.random.randn(*shape_ft) + 1j * np.random.randn(*shape_ft)
    

def test_hcm_init(shape):
    """Test Halfcomplex mapping init."""
    hcm = HalfComplexMapping(shape)
    assert isinstance(hcm.full_to_halfcomplex, np.ndarray)
    assert isinstance(hcm.halfcomplex_to_full, np.ndarray)

def test_halfcomplex_roundtrip(shape, halfcomplex_data):
    """Test unpack halfcomplex data."""
    hcm = HalfComplexMapping(shape)
    full_data = hcm.unpack(halfcomplex_data, shape)
    round_data = hcm.pack(full_data)
    np.testing.assert_allclose(round_data, halfcomplex_data, atol=1e-7)
    
def test_halfcomplex_compare(shape, halfcomplex_data):
    """Test python halfcomplex unpacking."""
    full_data_python = unpack_halfcomplex(halfcomplex_data, shape)
    hcm = HalfComplexMapping(shape)
    full_data = hcm.unpack(halfcomplex_data, shape)
    assert np.allclose(full_data, full_data_python)
    np.testing.assert_allclose(np.real(full_data), np.real(full_data_python), atol=1e-7)
    np.testing.assert_allclose(np.imag(full_data), np.imag(full_data_python), atol=1e-7)