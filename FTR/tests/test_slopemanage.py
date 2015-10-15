# -*- coding: utf-8 -*-
"""
Tests for slope management.
"""
import pytest
import numpy as np

from .. import slopemanage
from ..utils import circle_aperture, remove_piston, complexmp, shapegrid

@pytest.fixture
def shape():
    """Return a random shape."""
    nr = np.random.randint(10,20)
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
    """Random slopes."""
    sx = np.random.randn(*aperture.shape)
    sx *= aperture
    return sx

@pytest.fixture
def sy(aperture):
    """Random slopes."""
    sy = np.random.randn(*aperture.shape)
    sy *= aperture
    return sy
    
def test_compare_slopemanage(sx, sy, aperture, shape):
    """Compare slope management results."""
    np.set_printoptions(precision=4, linewidth=120, suppress=True)
    
    sx_c = np.copy(sx)
    sy_c = np.copy(sy)
    slopemanage.slope_management_fast(aperture, sy_c, sx_c)
    sy_p, sx_p = slopemanage.slope_management(aperture, sy, sx)
    np.testing.assert_allclose(sx_c, sx_p, atol=1e-7)
    np.testing.assert_allclose(sy_c, sy_p, atol=1e-7)

def test_compare_slopemanage_object(sx, sy, aperture, shape):
    """Compare slope management results."""
    np.set_printoptions(precision=4, linewidth=120, suppress=True)
    
    sx_c = np.copy(sx)
    sy_c = np.copy(sy)
    sm = slopemanage.SlopeManager(aperture)
    sm(sy_c, sx_c)
    sy_p, sx_p = slopemanage.slope_management(aperture, sy, sx)
    np.testing.assert_allclose(sx_c, sx_p, atol=1e-7)
    np.testing.assert_allclose(sy_c, sy_p, atol=1e-7)
    
def test_many_iteration_slopemanage(sx, sy, aperture, shape):
    """Test many iteration slope management."""
    sm = slopemanage.SlopeManager(aperture)
    for i in range(np.random.randint(2,30)):
        sm(sy, sx)
    