import pytest
import numpy as np


def compare_attributes_numpy(a, b):
    """Compare numpy array-like attributes."""
    dirb = set(dir(b))
    for attr in dir(a):
        if attr.startswith("__"):
            continue
        attra = getattr(a, attr, None)
        if not isinstance(attra, np.ndarray):
            continue
        assert attr in dirb, "b={0!r} missing attribute {1:s}".format(b, attr)
        attrb = getattr(b, attr)
        np.testing.assert_allclose(attra, attrb)
        dirb.remove(attr)
    
    for attr in dirb:
        if attr.startswith("__"):
            continue
        attrb = getattr(b, attr, None)
        if not isinstance(attrb, np.ndarray):
            continue
        assert attr in dir(a), "a={0!r} missing attribute {1:s}".format(a, attr)
        attra = getattr(a, attr)
        np.testing.assert_allclose(attra, attrb)

class IOTestBase(object):
    """Base class for IO tests."""
    
    cls = None
    
    @pytest.fixture
    def obj(self):
        """Return an instance of cls."""
        return self.cls()
    
    def pytest_generate_tests(self, metafunc):
           if 'format' in metafunc.fixturenames:
               formats = [ str(fmt) for fmt, read_write in self.cls.formats().items() if read_write == (True, True) ]
               metafunc.parametrize("format", formats)
    
    def test_roundtrip(self, obj, tmpdir, format):
        """Test a read-write filter roundtrip."""
        filename = str(tmpdir.join('test-{0}.{1}'.format(obj.__class__.__name__, format)))
        obj.write(filename)
        new_obj = self.cls.read(filename)
        compare_attributes_numpy(obj, new_obj)