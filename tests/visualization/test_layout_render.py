"""Objective browser tests for voxel layout and the split 2D/3D viewer."""

from __future__ import annotations

import asyncio
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image

from organoid_analysis.visualization.volume_viewer import (
    ChannelConfig,
    build_viewer_payload,
    render_viewer_html,
)
from tests.visualization._optional_deps import HAS_PLAYWRIGHT

ARTIFACTS = Path(__file__).parent / "artifacts"
ROOT = Path(__file__).parents[4]
SWIFTSHADER = [
    "--use-gl=angle",
    "--use-angle=swiftshader",
    "--enable-unsafe-swiftshader",
    "--ignore-gpu-blocklist",
]


async def _capture(spec, name: str, selectors=("#viewer",), timeout_s=45) -> tuple[dict, dict[str, Path]]:
    from playwright.async_api import async_playwright

    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    outputs = {selector: ARTIFACTS / f"{name}_{selector.strip('#').replace(' ', '_')}.png"
               for selector in selectors}
    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "viewer.html"
        source.write_text(render_viewer_html(spec), encoding="utf-8")
        async with async_playwright() as p:
            browser = await p.chromium.launch(headless=True, args=SWIFTSHADER)
            page = await browser.new_page(viewport={"width": 1100, "height": 760})
            rendered = asyncio.Event()
            errors: list[str] = []
            def on_console(msg):
                if "rendered mode=" in (msg.text or ""):
                    rendered.set()
                if msg.type in ("error", "warning"):
                    errors.append(f"{msg.type}: {msg.text}")
            page.on("console", on_console)
            page.on("pageerror", lambda err: errors.append(f"pageerror: {err}"))
            await page.goto(source.as_uri(), wait_until="load", timeout=90_000)
            try:
                await asyncio.wait_for(rendered.wait(), timeout=timeout_s)
            except TimeoutError as exc:
                raise AssertionError(f"viewer did not render: {errors}") from exc
            await page.wait_for_timeout(1500)
            diag = await page.evaluate("window.__viewer_diag(false)")
            await page.evaluate("""() => {
              document.getElementById('controls').style.display='none';
              document.getElementById('spacing').style.display='none';
            }""")
            for selector, output in outputs.items():
                await page.locator(selector).first.screenshot(path=str(output))
            await browser.close()
    return diag, outputs


def _signal_image(path: Path) -> np.ndarray:
    rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.float32)
    return np.maximum(rgb[..., 1], rgb[..., 2])


def _projected_volume(path: Path, diag: dict, out_shape: tuple[int, int]) -> np.ndarray:
    image = Image.open(path).convert("RGB")
    w, h = image.size
    bounds = diag["bounds"]
    scale = h / (2.0 * diag["camScale"])
    cx = (bounds[0] + bounds[1]) / 2.0
    cy = (bounds[2] + bounds[3]) / 2.0
    left = w / 2.0 + (bounds[0] - cx) * scale
    right = w / 2.0 + (bounds[1] - cx) * scale
    top = h / 2.0 - (bounds[3] - cy) * scale
    bottom = h / 2.0 - (bounds[2] - cy) * scale
    crop = image.crop((round(left), round(top), round(right), round(bottom)))
    crop = crop.resize((out_shape[1], out_shape[0]), Image.Resampling.BILINEAR)
    rgb = np.asarray(crop, dtype=np.float32)
    return np.maximum(rgb[..., 1], rgb[..., 2])


def _ncc(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = (a - a.mean()) / max(a.std(), 1e-9)
    b = (b - b.mean()) / max(b.std(), 1e-9)
    return float(np.mean(a * b))


@unittest.skipUnless(HAS_PLAYWRIGHT, "playwright not installed; headless layout render tests skipped")
class LayoutRenderTests(unittest.TestCase):
    def test_corner_voxel_projects_to_expected_xy(self) -> None:
        volume = np.zeros((32, 64, 64), dtype=np.float32)
        volume[5, 10, 20] = 1000.0
        spec = build_viewer_payload(
            [volume], (1.0, 1.0, 1.0),
            [ChannelConfig(lut="cyan", opacity=1.0, name="corner")],
            render_mode="mip",
        )
        diag, paths = asyncio.run(_capture(spec, "layout_corner"))
        signal = _signal_image(paths["#viewer"])
        mask = signal > 30
        self.assertTrue(mask.any(), "corner voxel did not produce a visible MIP point")
        yy, xx = np.nonzero(mask)
        weights = signal[mask]
        actual_x = float(np.average(xx, weights=weights))
        actual_y = float(np.average(yy, weights=weights))

        h, w = signal.shape
        bounds = diag["bounds"]
        px_per_world = h / (2.0 * diag["camScale"])
        cx = (bounds[0] + bounds[1]) / 2.0
        cy = (bounds[2] + bounds[3]) / 2.0
        expected_x = w / 2.0 + (20.0 - cx) * px_per_world
        expected_y = h / 2.0 - (10.0 - cy) * px_per_world
        error = float(np.hypot(actual_x - expected_x, actual_y - expected_y))
        print(f"\n  corner actual=({actual_x:.2f},{actual_y:.2f}) "
              f"expected=({expected_x:.2f},{expected_y:.2f}) error={error:.2f}px")
        self.assertLess(error, 8.0)

    def test_three_axis_gradient_orientation(self) -> None:
        z, y, x = np.indices((32, 64, 64))
        volume = (z * 100 + y + x / 100.0).astype(np.float32)
        spec = build_viewer_payload([volume], (1, 1, 1), render_mode="mip")
        diag, paths = asyncio.run(_capture(spec, "layout_gradient"))
        observed = _projected_volume(paths["#viewer"], diag, (64, 64))
        reference = volume.max(axis=0)
        ncc_direct = _ncc(observed, reference)
        ncc_y_flipped = _ncc(observed, np.flipud(reference))
        best = max(ncc_direct, ncc_y_flipped)
        orientation = "direct" if ncc_direct >= ncc_y_flipped else "global-y-flip"
        oriented = np.flipud(observed) if orientation == "global-y-flip" else observed
        y_corr = _ncc(oriented.mean(axis=1), np.arange(64))
        x_corr = _ncc(oriented.mean(axis=0), np.arange(64))
        print(f"\n  gradient ncc={best:.4f} orientation={orientation} "
              f"(direct={ncc_direct:.4f}, yflip={ncc_y_flipped:.4f}, "
              f"y_corr={y_corr:.4f}, x_corr={x_corr:.4f})")
        self.assertEqual(orientation, "global-y-flip")
        self.assertGreater(best, 0.75)
        self.assertGreater(y_corr, 0.75)
        self.assertGreater(x_corr, 0.0)

    def test_split_view_and_shared_lut(self) -> None:
        z, y, x = np.indices((32, 64, 64))
        volume = np.exp(-(((z-16)/8)**2+((y-32)/15)**2+((x-32)/15)**2)).astype(np.float32)
        spec = build_viewer_payload([volume], (1, 1, 2))

        async def capture_changed():
            from playwright.async_api import async_playwright
            ARTIFACTS.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp) / "viewer.html"
                source.write_text(render_viewer_html(spec), encoding="utf-8")
                async with async_playwright() as p:
                    browser = await p.chromium.launch(headless=True, args=SWIFTSHADER)
                    page = await browser.new_page(viewport={"width": 1100, "height": 760})
                    await page.goto(source.as_uri(), wait_until="load", timeout=90_000)
                    await page.wait_for_function("window.__viewer_state !== undefined", timeout=45_000)
                    await page.select_option("#lut-select", "magenta")
                    await page.locator("#slice-slider").fill("10")
                    await page.evaluate("""() => {
                      const el=document.getElementById('window-low'); el.value='0.2';
                      el.dispatchEvent(new Event('input',{bubbles:true}));
                    }""")
                    await page.wait_for_timeout(1200)
                    ui_state = await page.evaluate("window.__viewer_ui_state()")
                    await page.evaluate("document.getElementById('controls').style.display='none'")
                    left = ARTIFACTS / "split_2d_magenta.png"
                    right = ARTIFACTS / "split_3d_magenta.png"
                    await page.locator("#slice-canvas").screenshot(path=str(left))
                    await page.locator("#viewer").screenshot(path=str(right))
                    await browser.close()
            return left, right, ui_state

        left, right, ui_state = asyncio.run(capture_changed())
        self.assertEqual(ui_state["channels"][0]["lut"], "magenta")
        self.assertAlmostEqual(ui_state["channels"][0]["lo"], 0.2, places=4)
        self.assertEqual(ui_state["currentSlice"], 10)
        for label, path in (("2D", left), ("3D", right)):
            rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.int32)
            nonblack = float(((rgb.sum(axis=2)) > 15).mean())
            magenta = float(((rgb[..., 0] > 30) & (rgb[..., 2] > 30)
                             & (rgb[..., 1] < np.maximum(rgb[..., 0], rgb[..., 2]))).mean())
            print(f"\n  split {label}: non_black={nonblack:.4f} magenta={magenta:.4f}")
            self.assertGreater(nonblack, 0.005)
            self.assertGreater(magenta, 0.0)

    def test_path_a_quality_and_retina_invariants(self) -> None:
        z, y, x = np.indices((24, 64, 64))
        volumes = []
        for cx in (16, 32, 48):
            volumes.append(np.exp(-(((z-12)/6)**2+((y-32)/10)**2+((x-cx)/8)**2)).astype(np.float32))
        spec = build_viewer_payload(volumes, (0.414, 0.414, 1.5))

        async def inspect():
            from playwright.async_api import async_playwright
            ARTIFACTS.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp) / "viewer.html"
                source.write_text(render_viewer_html(spec), encoding="utf-8")
                async with async_playwright() as p:
                    browser = await p.chromium.launch(headless=True, args=SWIFTSHADER)
                    context = await browser.new_context(
                        viewport={"width": 1100, "height": 760}, device_scale_factor=2
                    )
                    page = await context.new_page()
                    await page.goto(source.as_uri(), wait_until="load", timeout=90_000)
                    await page.wait_for_function("window.__viewer_state !== undefined", timeout=45_000)
                    initial = await page.evaluate("""() => {
                      const s=window.__viewer_state,p=s.prop,m=s.mapper,c=document.querySelector('#viewer canvas');
                      const r=c.getBoundingClientRect();
                      return {implementation:s.implementation,volumeCount:s.volumes.length,
                        components:s.imageData.getPointData().getScalars().getNumberOfComponents(),
                        independent:p.getIndependentComponents(),linear:p.getInterpolationTypeAsString(),
                        shade:p.getShade(),unit:[0,1,2].map(i=>p.getScalarOpacityUnitDistance(i)),
                        auto:m.getAutoAdjustSampleDistances(),sample:m.getSampleDistance(),jitter:s.jittering,
                        dpr:window.devicePixelRatio,ratio:[c.width/r.width,c.height/r.height]};
                    }""")
                    await page.select_option("#quality-select", "smooth")
                    smooth = await page.evaluate("window.__viewer_state.mapper.getSampleDistance()")
                    await page.select_option("#quality-select", "sharp")
                    sharp = await page.evaluate("window.__viewer_state.mapper.getSampleDistance()")
                    path = ARTIFACTS / "multi_component_path_a.png"
                    await page.evaluate("document.getElementById('controls').style.display='none'")
                    await page.locator("#viewer").screenshot(path=str(path))
                    await browser.close()
            return initial, smooth, sharp, path

        initial, smooth, sharp, path = asyncio.run(inspect())
        print(f"\n  path A: {initial}, smooth={smooth}, sharp={sharp}, artifact={path}")
        self.assertEqual(initial["implementation"], "multi-component-single-volume")
        self.assertEqual(initial["volumeCount"], 1)
        self.assertEqual(initial["components"], 3)
        self.assertTrue(initial["independent"])
        self.assertEqual(initial["linear"], "LINEAR")
        self.assertFalse(initial["shade"])
        np.testing.assert_allclose(initial["unit"], [0.414] * 3)
        self.assertFalse(initial["auto"])
        self.assertAlmostEqual(initial["sample"], 0.2)
        self.assertIn(initial["jitter"], ("explicit", "vtkjs-webgl-built-in"))
        self.assertAlmostEqual(smooth, 0.5)
        self.assertAlmostEqual(sharp, 0.1)
        self.assertEqual(initial["dpr"], 2)
        np.testing.assert_allclose(initial["ratio"], [2, 2], atol=0.01)
        rgb = np.asarray(Image.open(path).convert("RGB"), dtype=np.int32)
        self.assertGreater(float((rgb.sum(axis=2) > 15).mean()), 0.005)
        r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
        color_fractions = {
            "cyan": float(((g > 10) & (b > 10) & (r < np.maximum(g, b))).mean()),
            "magenta": float(((r > 10) & (b > 10) & (g < np.maximum(r, b))).mean()),
            "yellow": float(((r > 10) & (g > 10) & (b < np.maximum(r, g))).mean()),
        }
        print(f"  path A component colors: {color_fractions}")
        for fraction in color_fractions.values():
            self.assertGreater(fraction, 0.001)

    def test_real_tiff_mip_matches_numpy_reference(self) -> None:
        import tifffile

        tiff_path = ROOT / "data" / "images" / "ZeroG-Breast-Cancer-Spheroid-C1.tif"
        if not tiff_path.exists():
            self.skipTest(
                "real data not present; restore via scripts/download_data.sh "
                f"(missing: {tiff_path})"
            )
        volume = tifffile.imread(tiff_path)
        self.assertEqual(volume.shape, (56, 512, 512))
        spec = build_viewer_payload(
            [volume.astype(np.float32)],
            (0.414, 0.414, 1.5),
            [ChannelConfig(lut="cyan", opacity=1.0, name="Nuclei")],
            render_mode="mip",
        )
        diag, paths = asyncio.run(_capture(spec, "spheroid_c1_mip", timeout_s=120))
        observed = _projected_volume(paths["#viewer"], diag, (512, 512))

        mip_ref = volume.max(axis=0).astype(np.float32)
        lo, hi = np.percentile(volume, (1, 99))
        reference = np.clip((mip_ref - lo) / max(hi - lo, 1), 0, 1)
        ref_rgb = np.zeros((512, 512, 3), dtype=np.uint8)
        ref_rgb[..., 1] = np.round(reference * 255).astype(np.uint8)
        ref_rgb[..., 2] = ref_rgb[..., 1]
        ref_path = ARTIFACTS / "spheroid_c1_mip.png"
        Image.fromarray(ref_rgb).save(ref_path)

        observed = observed / max(float(observed.max()), 1.0)
        direct = _ncc(observed, reference)
        y_flipped = _ncc(observed, np.flipud(reference))
        best = max(direct, y_flipped)
        orientation = "direct" if direct >= y_flipped else "global-y-flip"
        nonblack = float((observed > (30 / 255)).mean())
        print(f"\n  real TIFF MIP: ncc={best:.4f} orientation={orientation} "
              f"non_black={nonblack:.4f} ref={ref_path} viewer={paths['#viewer']}")
        self.assertEqual(orientation, "global-y-flip")
        self.assertGreater(best, 0.55)
        self.assertGreater(nonblack, 0.005)


if __name__ == "__main__":
    unittest.main(verbosity=2)
