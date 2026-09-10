from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from organoid_analysis.result_export.measurement_tables import export_results


def tables():
    return dict(organoids=pd.DataFrame({'organoid_id': [1]}), cells=pd.DataFrame({'cell_id': [1]}),
                nuclei=pd.DataFrame({'nucleus_id': [1]}), edges=pd.DataFrame({'cell_id_1': [1]}),
                qc_flags=pd.DataFrame({'flag': ['example']}), summary={'status': 'complete'})


def test_invalid_export_never_publishes_complete_summary(tmp_path):
    output = tmp_path / 'result'
    with pytest.raises(ValueError):
        export_results(output, **tables(), label_masks={'organoid': np.full((3, 4, 4), 2**32, np.uint64)},
                       spacing_zyx_um=(1, 1, 1))
    assert not output.exists()


def test_mask_write_failure_leaves_incomplete_marker_not_complete_summary(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError('simulated disk write failure')
    monkeypatch.setattr('organoid_analysis.microscopy_io.tiff_contract.write_labels', fail)
    output = tmp_path / 'result'
    with pytest.raises(OSError, match='disk write'):
        export_results(output, **tables(), label_masks={'organoid': np.ones((3, 4, 4), np.uint32)},
                       spacing_zyx_um=(1, 1, 1))
    assert (output / 'RUN_INCOMPLETE.txt').exists()
    assert not (output / 'summary' / 'analysis_summary.json').exists()


def test_saved_segmentation_write_failure_cannot_be_restored_as_complete(tmp_path, monkeypatch):
    from organoid_analysis.segmentation import cellpose_inference as ci
    def fail(*args, **kwargs):
        raise OSError('archive write failed')
    monkeypatch.setattr(ci.zipfile.ZipFile, 'write', fail)
    with pytest.raises(OSError, match='archive write'):
        ci.save_result(np.ones((3, 4, 4), np.uint32), None, ci.SegmentationConfig(), tmp_path)
    assert ci.list_saved_results(tmp_path) == []
    assert len(list(tmp_path.glob('run-*/RUN_INCOMPLETE.txt'))) == 1


def test_source_hashes_cover_spacing_and_distinct_config_modules():
    import hashlib
    from pathlib import Path

    from organoid_analysis.microscopy_io import tiff_contract
    hashes = tiff_contract.source_code_hashes()
    root = Path(tiff_contract.__file__).resolve().parent.parent
    for relative in ('microscopy_io/metadata.py', 'microscopy_io/voxel_spacing.py',
                     'config.py', 'quantification/multilevel_relationships/config.py'):
        assert hashes[relative] == hashlib.sha256((root / relative).read_bytes()).hexdigest()


def test_saved_masks_are_calibrated_zyx_not_rgb(tmp_path):
    import tifffile

    from organoid_analysis.microscopy_io.tiff_contract import read_tiff
    from organoid_analysis.segmentation import cellpose_inference as ci

    mask = np.arange(48, dtype=np.uint32).reshape(3, 4, 4)
    result = ci.save_result(mask, mask, ci.SegmentationConfig(xy_spacing_um=0.5, anisotropy=3), tmp_path)
    for path in (result.nuclei_mask_path, result.cell_mask_path):
        with tifffile.TiffFile(path) as handle:
            assert handle.series[0].axes == 'ZYX'
            assert handle.pages[0].samplesperpixel == 1
        data, spacing, _ = read_tiff(path)
        np.testing.assert_array_equal(data[0], mask)
        assert spacing == (1.5, 0.5, 0.5)
