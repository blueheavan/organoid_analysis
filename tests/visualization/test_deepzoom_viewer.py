"""Deep-zoom plane viewer: pyramid geometry, contour tiles, budget, HTML contract and a headless render."""
from __future__ import annotations

import base64
import glob
import io
import json
import os
import re
import unittest

import numpy as np
import pytest
from PIL import Image

from organoid_analysis.visualization.deepzoom_viewer import (
    build_deepzoom_payload,
    level_shape,
    max_level_for,
    render_deepzoom_html,
)
from organoid_analysis.visualization.deepzoom_viewer.deepzoom_payload import (
    _VENDOR,
    contour_labels,
    pyramid,
    spec_to_dict,
)
from tests.visualization._optional_deps import HAS_PLAYWRIGHT


def _decode(url: str) -> np.ndarray:
    return np.asarray(Image.open(io.BytesIO(base64.b64decode(url.split(",", 1)[1]))))


def _stack(shape=(5, 300, 520)):
    z, y, x = np.indices(shape)
    image = (100 + 10 * z + (x % 50) * 4).astype(np.uint16)
    labels = np.zeros(shape, np.int32)
    labels[:, 40:120, 60:200] = 3
    labels[1:3, 200:280, 300:500] = 8
    return image, labels


@pytest.mark.parametrize("width, height, top", [(1, 1, 0), (256, 256, 8), (257, 10, 9), (520, 300, 10)])
def test_level_geometry_matches_openseadragon(width, height, top):
    assert max_level_for(width, height) == top
    assert level_shape(width, height, top, top) == (height, width)
    assert level_shape(width, height, 0, top) == (1, 1)


def test_pyramid_halves_by_area_mean_and_keeps_contours_by_max():
    base = np.arange(25, dtype=float).reshape(5, 5)
    levels = pyramid(base, 3, "mean")
    assert [level.shape for level in levels] == [(1, 1), (2, 2), (3, 3), (5, 5)]
    assert levels[2][0, 0] == pytest.approx(base[:2, :2].mean())
    assert levels[2][2, 2] == pytest.approx(24.0)          # edge padding, not zero padding
    contours = np.zeros((5, 5), int)
    contours[4, 4] = 7
    assert pyramid(contours, 3, "max")[0][0, 0] == 7        # a 1-px contour survives to level 0


def test_every_level_is_fully_tiled_and_tiles_reassemble_the_windowed_plane():
    image, labels = _stack()
    spec = build_deepzoom_payload(image, labels, (3, 0.65, 0.65), centre_plane=2, window=(100.0, 336.0))
    plane = spec.planes[2]
    assert len(plane["image"]) == spec.max_level + 1
    for level, tiles in enumerate(plane["image"]):
        rows, cols = level_shape(spec.width, spec.height, level, spec.max_level)
        assert len(tiles) == -(-rows // 256) * -(-cols // 256)
    top = plane["image"][-1]
    mosaic = np.block([[_decode(top[f"{c}_{r}"]) for c in range(3)] for r in range(2)])
    expected = np.clip(np.rint((image[2].astype(float) - 100) / 236 * 255), 0, 255).astype(np.uint8)
    np.testing.assert_array_equal(mosaic, expected)          # PNG is lossless, pixels are native


def test_contour_tiles_mark_inner_boundaries_with_label_colours_and_skip_empty_tiles():
    image, labels = _stack()
    spec = build_deepzoom_payload(image, labels, centre_plane=0)
    overlay = spec.planes[0]["overlay"][-1]
    assert set(overlay) == {"0_0"}                           # label 8 is absent at z 0; other tiles empty
    rgba = _decode(overlay["0_0"])
    edge = contour_labels(labels[0])[:256, :256] > 0
    np.testing.assert_array_equal(rgba[..., 3] > 0, edge)
    assert rgba[40, 100, 3] == 255 and rgba[80, 100, 3] == 0  # boundary row vs interior
    assert spec.planes[0]["ids"] == [[3, pytest.approx(130.0), pytest.approx(80.0)]]


def test_budget_keeps_nearest_planes_and_always_the_centre():
    image, labels = _stack((9, 300, 520))
    full = build_deepzoom_payload(image, labels, centre_plane=4)
    assert full.plane_indices == list(range(9))
    per_plane = full.payload_bytes / 9
    limited = build_deepzoom_payload(image, labels, centre_plane=4, budget_bytes=int(per_plane * 3.5))
    assert limited.plane_indices == [3, 4, 5]
    assert build_deepzoom_payload(image, labels, centre_plane=0, budget_bytes=1).plane_indices == [0]
    assert build_deepzoom_payload(image, labels, centre_plane=8, max_planes=2).plane_indices == [7, 8]


def test_invalid_inputs_are_rejected():
    image, labels = _stack((2, 20, 20))
    with pytest.raises(ValueError, match="must match"):
        build_deepzoom_payload(image, labels[:, :10])
    with pytest.raises(ValueError, match="nonnegative integer"):
        build_deepzoom_payload(image, labels.astype(float))
    with pytest.raises(ValueError, match="centre_plane"):
        build_deepzoom_payload(image, centre_plane=2)
    with pytest.raises(ValueError, match="positive finite"):
        build_deepzoom_payload(image, spacing_zyx_um=(1, 0, 1))
    with pytest.raises(ValueError, match="high > low"):
        build_deepzoom_payload(image, window=(5.0, 5.0))


def test_html_is_self_contained_with_one_script_and_the_vendored_bundle():
    image, labels = _stack((2, 40, 60))
    spec = build_deepzoom_payload(image, labels, (2.0, 0.5, 0.5))
    html = render_deepzoom_html(spec)
    assert html.count("<script") == 1 and html.count("</script>") == 1
    assert "openseadragon 6.1.1" in html and "sourceMappingURL" not in html
    assert not re.search(r"""(src|href)=["']https?://""", html)
    payload = json.loads(json.dumps(spec_to_dict(spec)))
    assert payload["plane_indices"] == [0, 1] and payload["spacing_zyx_um"] == [2.0, 0.5, 0.5]


def test_vendored_bundle_matches_recorded_provenance():
    import hashlib

    digest = hashlib.sha256(_VENDOR.read_bytes()).hexdigest()
    note = (_VENDOR.parent / "openseadragon.PROVENANCE.txt").read_text(encoding="utf-8")
    assert digest in note
    assert "BSD" in note and (_VENDOR.parent / "openseadragon.js.LICENSE.txt").is_file()


def _launch(playwright):
    try:
        return playwright.chromium.launch(headless=True)
    except Exception:
        root = os.environ.get("PLAYWRIGHT_BROWSERS_PATH", "")
        candidates = sorted(glob.glob(os.path.join(root, "chromium-*", "chrome-linux", "chrome"))) if root else []
        if not candidates:
            raise unittest.SkipTest("no Chromium available for Playwright") from None
        return playwright.chromium.launch(headless=True, executable_path=candidates[-1])


@unittest.skipUnless(HAS_PLAYWRIGHT, "playwright not installed; headless deep-zoom render test skipped")
def test_headless_render_readout_and_plane_switch():
    from playwright.sync_api import sync_playwright

    image, labels = _stack((3, 300, 520))
    html = render_deepzoom_html(build_deepzoom_payload(image, labels, (3.0, 0.65, 0.65), centre_plane=1))
    loaded = ("window.__deepzoom && window.__deepzoom.viewer.world.getItemCount() === 2"
              " && window.__deepzoom.viewer.world.getItemAt(0).getFullyLoaded()")
    errors: list[str] = []
    with sync_playwright() as playwright:
        browser = _launch(playwright)
        page = browser.new_page(viewport={"width": 900, "height": 560})
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.set_content(html)
        page.wait_for_function(loaded, timeout=20000)
        position = page.evaluate("""() => { const v = window.__deepzoom.viewer, it = v.world.getItemAt(0);
            const p = it.imageToViewerElementCoordinates(new OpenSeadragon.Point(100.5, 50.5));
            const r = v.container.getBoundingClientRect(); return [p.x + r.left, p.y + r.top]; }""")
        page.mouse.move(*position)
        readout = page.inner_text("#readout")
        assert "(1, 50, 100)" in readout
        assert "y 32.50" in readout and "x 65.00" in readout and "z 3.00" in readout
        page.check("#ids")
        assert page.evaluate("document.querySelectorAll('.idtag').length") == 2
        page.evaluate("const s = document.getElementById('zslider'); s.value = '0'; s.dispatchEvent(new Event('input'))")
        page.wait_for_function(loaded, timeout=20000)
        assert page.inner_text("#zlabel").startswith("z 0 / 2")
        assert page.evaluate("document.querySelectorAll('.idtag').length") == 1
        browser.close()
    assert errors == []
