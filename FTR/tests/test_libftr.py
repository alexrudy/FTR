# -*- coding: utf-8 -*-
"""
Test the fourier transform reconstructor from FTRLib bound via cython.
"""
import pytest
import numpy as np

from .test_reconstructor_base import ReconstructorTestBase
from .test_ftr import FourierTransformReconstructorTestBase, filter_name
from .. import ftr
from ..ftr import mod_hud_filter
from ..utils import circle_aperture, remove_piston, complexmp


class TestFTRLibCython(FourierTransformReconstructorTestBase):
    """Test the Fourier Transform Reconstructor."""
    
    cls = ftr.FastFTReconstructor
    
    repr = "<CFTRBase {shape}>"
    
    def test_compare_to_python(self, ap, filter_name, phase, shape):
        """Test comparison to the python-only filter."""
        reconstructor = self.cls(ap, filter=filter_name)
        reconstructor_py = ftr.FourierTransformReconstructor(ap, filter=filter_name)
        xs, ys = reconstructor.invert(phase)
        phase_rt = reconstructor(xs, ys)
        phase_py = reconstructor_py(xs, ys)
        np.testing.assert_allclose(phase_rt, phase_py, atol=1e-7)
        