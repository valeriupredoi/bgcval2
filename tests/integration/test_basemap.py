"""Deprecate when basemap gets retired."""
import numpy as np

try:
    from mpl_toolkits.basemap import Basemap
except ImportError as exc:
    print("Basemap is no longer supported, please update your scripts.")
    raise exc


def test_basemap_functionality():
    """Test a basic basemap functionality and print deprecation warning."""
    m1 = Basemap(projection='robin', lon_0=-106., resolution='c')
    x1, y1 = m1(np.array([1, 2]), np.array([1, 2]))
    expected_x1 = np.array([27083811.675015, 27176210.433604])
    expected_y1 = np.array([8722331.436056, 8829163.62464])
    np.testing.assert_allclose(x1, expected_x1)
    np.testing.assert_allclose(y1, expected_y1)
