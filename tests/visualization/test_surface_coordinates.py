"""Check native mesh coordinates without requiring a graphics context."""
from __future__ import annotations

import numpy as np
import pytest

from organoid_analysis.microscopy_io.metadata import Spacing
from organoid_analysis.visualization.surface_rendering import build_surface, surface_to_dict


@pytest.mark.parametrize('label', [None, 513])
def test_surface_preserves_asymmetric_zyx_locations_and_sparse_label(label):
    volume = np.zeros((7, 11, 17), np.uint16)
    volume[2:5, 3:7, 9:14] = 1 if label is None else label
    surfaces = build_surface(volume, Spacing(2, 3, 5), threshold=.5,
                             labels=None if label is None else [label])
    mesh = surfaces[0].polydata
    # Level .5 intersects halfway between occupied and background voxel centres.
    expected_xyz_bounds = (8.5 * 2, 13.5 * 2, 2.5 * 3, 6.5 * 3, 1.5 * 5, 4.5 * 5)
    assert tuple(mesh.bounds) == pytest.approx(expected_xyz_bounds, rel=1e-12)
    if label is not None:
        assert set(surface_to_dict(surfaces[0])['polydata']['scalars']) == {label}
