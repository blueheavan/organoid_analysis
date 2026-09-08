from __future__ import annotations

import sys
from pathlib import Path

import nbformat
from nbconvert.exporters import PythonExporter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_MARKERS = ("google.colab", "/content/drive", "drive.mount(")


def validate_notebook(path: Path) -> list[str]:
    notebook = nbformat.read(path, as_version=4)
    source = "\n".join(
        "".join(cell.source)
        for cell in notebook.cells
        if cell.cell_type == "code"
    )
    errors = [
        f"{path.name}: contains unsupported Colab reference '{marker}'"
        for marker in FORBIDDEN_MARKERS
        if marker in source
    ]

    script, _ = PythonExporter().from_notebook_node(notebook)
    try:
        compile(script, path.name, "exec")
    except SyntaxError as error:
        errors.append(f"{path.name}: generated Python is invalid: {error}")
    return errors


def main() -> int:
    notebooks = sorted((PROJECT_ROOT / "notebooks").glob("*.ipynb"))
    if not notebooks:
        print("No notebooks found.", file=sys.stderr)
        return 1

    errors = [error for path in notebooks for error in validate_notebook(path)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1

    print(f"Validated {len(notebooks)} notebooks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
