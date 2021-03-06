# -*- coding: utf-8 -*-

import pytest
import numpy as np
import random

from .. import utils

@pytest.fixture
def shape():
    """Make a random shape"""
    return (random.randint(10, 128), random.randint(10, 128))

@pytest.fixture
def square_shape():
    """Make a square shape."""
    n = random.randint(10, 128)
    return (n, n)
    
def complex_data(shape):
    """Make complex data of a given shape."""
    return np.random.randn(*shape) + 1j * np.random.randn(*shape)

@pytest.fixture
def random_data(shape):
    """Random data."""
    return np.random.randn(*shape)
    
@pytest.fixture
def aperture(shape):
    """An aperture."""
    ap = np.ones(shape, dtype=np.bool)
    ap[:2,:] = 0.0
    ap[:,:2] = 0.0
    ap[-2:,:] = 0.0
    ap[:,-2:] = 0.0
    return ap

def test_shapegrid(shape):
    """Test the shape grid"""
    yy, xx = utils.shapegrid(shape)
    assert yy.shape == shape
    assert xx.shape == shape

def test_shapegrid_not_centered(shape):
    """Test the fourier transform grid."""
    yy, xx = utils.shapegrid(shape, centered=False)
    
    assert np.allclose(yy[0,:], 0.0)
    assert np.allclose(xx[:,0], 0.0)
    assert np.allclose(yy[6,:], 6.0)
    assert np.allclose(xx[:,6], 6.0)
    
def test_fftgrid(shape):
    """Test the fft grid"""
    fy, fx = utils.fftgrid(shape)
    assert fy.shape == shape
    assert fx.shape == shape
    
def test_remove_piston(shape, random_data):
    """Random data."""
    fap = np.ones(shape, dtype=np.bool)
    d_pr, p = utils.remove_piston(fap, random_data)
    assert np.allclose(np.mean(random_data), p)
    assert np.allclose(np.mean(d_pr), 0.0)
    
def test_remove_piston_ap(shape, random_data, aperture):
    """Test remove piston with an aperture applied."""
    d_pr, p = utils.remove_piston(aperture, random_data)
    assert np.allclose(np.mean(random_data[aperture]), p)
    assert np.allclose(np.mean(d_pr[aperture]), 0.0)
    
def test_remove_piston_known(shape, random_data, aperture):
    """Test remove a known piston shape."""
    random_data -= np.mean(random_data[aperture])
    random_data[aperture] += 10000
    d_pr, p = utils.remove_piston(aperture, random_data)
    assert np.allclose(p, 10000)
    
def test_remove_tiptilt(shape, random_data, aperture):
    """Test removing a known tip-tilt"""
    random_data, _ttx, _tty = utils.remove_tiptilt(aperture, random_data)
    y, x = utils.shapegrid(shape)
    ttx, tty = random.uniform(-10, 10), random.uniform(-10, 10)
    ttx, tty = -5, 10
    tilted_data = utils.apply_tiptilt(aperture, random_data, ttx, tty)
    
    ttr, ttxr, ttyr = utils.remove_tiptilt(aperture, tilted_data)
    assert np.allclose(ttxr, ttx)
    assert np.allclose(ttyr, tty)
    assert np.allclose(ttr, random_data)
    
def test_make_hermitian(square_shape):
    """Test make a hermitian matrix."""
    data = complex_data(square_shape)
    data_h = utils.make_hermitian(data)
    assert utils.is_hermitian(data_h)
    