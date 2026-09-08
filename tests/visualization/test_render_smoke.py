"""Headless render smoke test: the synthetic volume must NOT come out black.

Renders a standalone viewer HTML (via render_viewer_html) in headless Chromium
+ SwiftShader, screenshots the viewer, and asserts via Pillow that a meaningful
non-black, cyan-tinted volume is present. Guards against the silent "all-black
framebuffer" regression (caused previously by the multi-component independent
compositing path).
"""

from __future__ import annotations

import asyncio
import unittest
from pathlib import Path

import numpy as np
from PIL import Image

from organoid_analysis.visualization.volume_viewer.viewer_payload import (
    build_viewer_payload,
    render_viewer_html,
)
from tests.visualization._optional_deps import HAS_PLAYWRIGHT

ARTIFACTS = Path(__file__).parent / "artifacts"


def _sphere_stack(shape=(48, 64, 64)) -> np.ndarray:
    from tests.visualization.test_viewer import _sphere_stack as _mk
    return _mk(shape)


async def _render_and_capture(html_path: Path, png: Path) -> tuple[dict, Path]:
    from playwright.async_api import async_playwright

    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    async with async_playwright() as p:
        browser = await p.chromium.launch(
            headless=True,
            args=[
                "--use-gl=angle",
                "--use-angle=swiftshader",
                "--enable-unsafe-swiftshader",
                "--ignore-gpu-blocklist",
            ],
        )
        page = await browser.new_page(viewport={"width": 1000, "height": 900})
        rendered = asyncio.Event()

        def _on_console(msg):
            if "rendered mode=" in (msg.text or ""):
                rendered.set()

        page.on("console", _on_console)
        await page.goto(html_path.as_uri(), wait_until="load", timeout=60_000)
        try:
            await asyncio.wait_for(rendered.wait(), timeout=25)
        except TimeoutError:
            pass
        await page.wait_for_timeout(2500)
        await page.evaluate("""() => {
          document.getElementById('controls').style.display='none';
          document.getElementById('spacing').style.display='none';
        }""")
        await page.locator("#viewer").first.screenshot(path=str(png))
        await browser.close()
    return _analyze(png), png


def _analyze(png: Path) -> dict:
    img = Image.open(png).convert("RGB")
    a = np.asarray(img).astype(np.int32)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    total = a.shape[0] * a.shape[1]
    nonblack = int(((r + g + b) > 15).sum())
    cyan = int(((g > 30) & (b > 30) & (r < np.maximum(g, b))).sum())
    return {
        "size": [a.shape[1], a.shape[0]],
        "non_black_fraction": float(nonblack) / total,
        "cyan_fraction": float(cyan) / total,
        "max_px": int(a.max()),
    }


@unittest.skipUnless(HAS_PLAYWRIGHT, "playwright not installed; headless render smoke tests skipped")
class RenderSmokeTests(unittest.TestCase):
    def _assert_renders(self, spec, name: str, out_name: str) -> tuple[dict, Path]:
        html = render_viewer_html(spec, inline_vtkjs=True)
        html_path = ARTIFACTS / f"{name}_source.html"
        html_path.write_text(html, encoding="utf-8")
        png = ARTIFACTS / f"{out_name}.png"
        stats, _ = asyncio.run(_render_and_capture(html_path, png))
        print(f"\n  {name}: non_black={stats['non_black_fraction']:.4f} "
              f"cyan={stats['cyan_fraction']:.5f} max_px={stats['max_px']} -> {png}")
        self.assertGreater(
            stats["non_black_fraction"], 0.005,
            f"{name} render came out black: {stats} (see {png})",
        )
        self.assertGreater(
            stats["cyan_fraction"], 0.0,
            f"no cyan pixels in {name}: {stats} (see {png})",
        )
        return stats, png

    def test_synthetic_volume_is_not_black(self) -> None:
        vols = [_sphere_stack((48, 64, 64)).astype(np.float32)]
        spec = build_viewer_payload(
            vols,
            (0.65, 0.65, 2.0),
            channels=None,
        )
        self._assert_renders(spec, "smoke_synthetic", "smoke_synthetic")

    def test_sparse_tiff_like_volume_is_not_black(self) -> None:
        # model the real uploaded TIFF: (Z,Y,X) uint16, ~94.8% zeros, strong tail.
        rng = np.random.default_rng(0)
        shape = (56, 256, 256)
        vol = np.zeros(shape, np.uint16)
        mask = rng.random(shape) < 0.052
        vol[mask] = rng.integers(1000, 55000, size=int(mask.sum())).astype(np.uint16)
        spec = build_viewer_payload(
            [vol.astype(np.float32)], (0.414, 0.414, 1.5), channels=None
        )
        self._assert_renders(spec, "smoke_sparse", "smoke_sparse")

    def test_mask_overlay_renders_without_error(self) -> None:
        from organoid_analysis.visualization.volume_viewer.viewer_payload import spec_to_dict

        vol = _sphere_stack((48, 64, 64)).astype(np.float32)
        labels = np.zeros(vol.shape, dtype=np.uint32)
        # two labelled blobs
        z, y, x = np.mgrid[0:48, 0:64, 0:64]
        labels[(z - 16) ** 2 + (y - 32) ** 2 + (x - 24) ** 2 < 6.0**2] = 1
        labels[(z - 32) ** 2 + (y - 30) ** 2 + (x - 44) ** 2 < 5.0**2] = 2
        spec = build_viewer_payload(
            [vol], (0.65, 0.65, 2.0), channels=None, mask_overlay=labels
        )
        self.assertIsNotNone(spec.mask)
        spec_dict = spec_to_dict(spec)
        html = render_viewer_html(spec, inline_vtkjs=True)
        html_path = ARTIFACTS / "mask_overlay_source.html"
        html_path.write_text(html, encoding="utf-8")
        png = ARTIFACTS / "mask_overlay.png"

        stats, _ = asyncio.run(_render_and_capture(html_path, png))
        print(f"\n  mask_overlay: non_black={stats['non_black_fraction']:.4f} "
              f"max_px={stats['max_px']} -> {png}")
        self.assertGreater(stats["non_black_fraction"], 0.005)
        self.assertIn("mask", spec_dict)
        self.assertEqual(spec_dict["mask"]["shape"], list(labels.shape))


if __name__ == "__main__":
    unittest.main(verbosity=2)
