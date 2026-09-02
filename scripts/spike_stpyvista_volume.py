"""Spike: does the HTML serialization chain preserve volume actors?

Tests both serialization backends used by stpyvista:
  1. Panel VTK pane -> save() to HTML
  2. PyVista plotter.export_html() (trame/vtk.js)

For each backend we render:
  - a VOLUME actor (the crux: stpyvista's export chain is claimed not to
    serialise volume actors)
  - a SURFACE (polydata) actor as a control

We then inspect the output HTML for evidence that the volume data made it
through (volume array / scalar file / image data) vs an empty shell.

This is the feasibility proof requested as "Step 0" before any UI work.
"""

from io import StringIO

import pyvista as pv


def probe_volume_actor():
    pl = pv.Plotter(off_screen=True, window_size=(300, 300))
    pl.add_volume(pv.Wavelet(extent=(-8, 8, -8, 8, -8, 8)), mapper="gpu", cmap="viridis")
    pl.camera_position = "iso"
    return pl


def probe_surface_actor():
    pl = pv.Plotter(off_screen=True, window_size=(300, 300))
    pl.add_mesh(pv.Sphere())
    pl.camera_position = "iso"
    return pl


def inspect(label, html: str) -> None:
    size = len(html)
    markers = {
        '"vtkVolume"': "vtkVolume typed",
        "vtkVolume": "vtkVolume",
        'application/vnd.vtk': "mime-vtk",
        '"ArrayName"': "ArrayName",
        "scalars": "scalars",
        "wavelet": "wavelet(probe name)",
        "vtkImageData": "ImageData",
        'polydata': "polydata",
        '"type":"vtkImageData"': "ImageData(typed)",
        '"type":"vtkPolyData"': "PolyData(typed)",
    }
    found = [name for hay, name in markers.items() if hay in html]
    # Count how much raw binary/array data is embedded (indexed_array_buffer / compressed)
    print(f"[{label}] html bytes={size}")
    print(f"    markers found: {found if found else 'NONE'}")
    print(f"    head: {html[:120]!r}")


def test_panel_backend():
    import panel as pn

    print("=== PANEL backend ===")
    for label, make in [
        ("VOLUME", probe_volume_actor),
        ("SURFACE", probe_surface_actor),
    ]:
        try:
            pl = make()
            pane = pn.pane.VTK(pl.ren_win, width=300, height=300)
            with StringIO() as buf:
                pane.save(buf, title="spike")
                html = buf.getvalue()
            inspect(f"PANEL {label}", html)
            # Save for manual inspection
            with open(f"/tmp/spike_panel_{label.lower()}.html", "w") as fh:
                fh.write(html)
            pl.close()
        except Exception as e:  # noqa: BLE001
            print(f"[PANEL {label}] ERROR: {type(e).__name__}: {e}")


def test_trame_export_html():
    print("=== TRAME export_html ===")
    for label, make in [
        ("VOLUME", probe_volume_actor),
        ("SURFACE", probe_surface_actor),
    ]:
        try:
            pl = make()
            raw_html = pl.export_html(filename=None)
            html = raw_html.read().decode("utf-8", errors="ignore")
            inspect(f"TRAME {label}", html)
            pl.close()
        except Exception as e:  # noqa: BLE001
            print(f"[TRAME {label}] ERROR: {type(e).__name__}: {e}")


if __name__ == "__main__":
    test_panel_backend()
    test_trame_export_html()
