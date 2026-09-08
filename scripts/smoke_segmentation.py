from __future__ import annotations

from pathlib import Path
import sys
import tempfile

import tifffile


def _bootstrap() -> None:
    _SRC = Path(__file__).resolve().parent.parent / "src"
    if str(_SRC) not in sys.path:
        sys.path.insert(0, str(_SRC))


_bootstrap()

DATA_IMAGES_DIR = Path(__file__).resolve().parent.parent / "data" / "images"

from organoid_analysis.segmentation.cellpose_inference import (  # noqa: E402
    SegmentationConfig,
    create_model,
    get_accelerator,
    save_result,
    segment_stacks,
)


def main() -> None:
    images_directory = DATA_IMAGES_DIR
    nuclei_path = images_directory / "PDAC-C1.tif"
    cell_path = images_directory / "PDAC-C2.tif"
    if not nuclei_path.exists() or not cell_path.exists():
        print(
            "SMOKE TEST SKIPPED: real microscopy data not found at\n"
            f"  {nuclei_path}\n"
            f"  {cell_path}\n"
            "Restore via:  pixi run data   (or: bash scripts/download_data.sh)",
            file=sys.stderr,
        )
        return
    nuclei_stack = tifffile.imread(nuclei_path)
    cell_stack = tifffile.imread(cell_path)
    z0 = max(nuclei_stack.shape[0] // 2 - 1, 0)
    y0 = max(nuclei_stack.shape[1] // 2 - 128, 0)
    x0 = max(nuclei_stack.shape[2] // 2 - 128, 0)
    z1 = min(z0 + 3, nuclei_stack.shape[0])
    y1 = min(y0 + 256, nuclei_stack.shape[1])
    x1 = min(x0 + 256, nuclei_stack.shape[2])
    nuclei_crop = nuclei_stack[z0:z1, y0:y1, x0:x1]
    cell_crop = cell_stack[z0:z1, y0:y1, x0:x1]

    print(f"Accelerator: {get_accelerator()}")
    print(f"Smoke-test shape: {nuclei_crop.shape}")
    masks = segment_stacks(
        create_model(),
        nuclei_crop,
        cell_crop,
        SegmentationConfig(batch_size=1),
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        result = save_result(*masks, SegmentationConfig(batch_size=1), Path(temporary_directory))
        assert result.archive_path.exists()
        print(f"Nuclei masks: {result.nuclei_count}")
        print(f"Cell masks: {result.cell_count}")
        print("3D segmentation smoke test passed.")


if __name__ == "__main__":
    main()
