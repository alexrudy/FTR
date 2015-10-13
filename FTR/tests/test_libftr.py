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
    
