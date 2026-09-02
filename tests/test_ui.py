from pathlib import Path
import subprocess
import sys
import unittest


class UiTests(unittest.TestCase):
    def test_ui_package_resolves_core_analysis_when_ui_path_is_first(self) -> None:
        """Mirror ``streamlit run src/ui/app.py`` path precedence.

        The UI has an ``analysis.py`` helper next to the core ``analysis``
        package. Feature extraction must still reach ``analysis.features``.
        """
        root = Path(__file__).resolve().parents[1]
        script = f'''import sys
sys.path.insert(0, {str(root / "src")!r})
sys.path.insert(0, {str(root / "src" / "ui")!r})
import numpy as np
from ui.vtk_viewer.mask_features import extract_mask_features
labels = np.zeros((4, 6, 6), dtype=np.uint16)
labels[1:3, 2:4, 2:4] = 1
assert len(extract_mask_features(labels)) == 1
'''
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, check=False
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
