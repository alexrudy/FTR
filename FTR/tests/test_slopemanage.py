# -*- coding: utf-8 -*-
"""
Tests for slope management.
"""
import pytest
import numpy as np

from ..libftr import _slopemanage
from .. import slopemanage
from ..utils import circle_aperture, remove_piston, complexmp, shapegrid

@pytest.fixture
def shape():
    """Return a random shape."""
    nr = 6
    # nr = np.random.randint(10,20)
    # nc = np.random.randint(10,128)
    return (nr, nr)
    
@pytest.fixture
def aperture(shape):
    """Aperture"""
    nr = min(shape)
    r = (nr // 2) - 1
    return circle_aperture(shape, r).astype(np.int).copy()

@pytest.fixture
def phase(aperture):
    """Random phase."""
    phase = np.random.randn(*aperture.shape)
    phase *= aperture
    return phase

@pytest.fixture
def sx(aperture):
    """Random phase."""
    sx = np.random.randn(*aperture.shape)
    sx *= aperture
    return sx

@pytest.fixture
def sy(aperture):
    """Random phase."""
    sy = np.random.randn(*aperture.shape)
    sy *= aperture
    return sy
    
def test_compare_slopemanage(sx, sy, aperture, shape):
    """Compare slope management results."""
    np.set_printoptions(precision=4, linewidth=120, suppress=True)
    
    sx_c = np.copy(sx)
    sy_c = np.copy(sy)
    _slopemanage.do_slope_management(aperture, sy_c, sx_c)
    sy_p, sx_p = slopemanage.slope_management(aperture, sy, sx)
    np.testing.assert_allclose(sx_c, sx_p, atol=1e-7)
    np.testing.assert_allclose(sy_c, sy_p, atol=1e-7)
    