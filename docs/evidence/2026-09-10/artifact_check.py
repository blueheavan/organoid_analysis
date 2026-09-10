"""Build/install the current wheel and execute from outside the checkout."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from setuptools.build_meta import build_wheel

root = Path.cwd()
evidence = root / 'docs/evidence/2026-09-10'
work = Path(tempfile.mkdtemp(prefix='organoid-final-artifact-'))
stage, wheel_dir, installed = work / 'source', work / 'wheel', work / 'installed'
wheel_dir.mkdir(parents=True)
shutil.copytree(root / 'src', stage / 'src', ignore=shutil.ignore_patterns('__pycache__', '*.egg-info'))
for name in ('pyproject.toml', 'README.md', 'LICENSE'):
    if (root / name).exists():
        shutil.copy2(root / name, stage / name)
os.chdir(stage)
wheel = wheel_dir / build_wheel(str(wheel_dir))
os.chdir(root)
env = dict(os.environ, UV_CACHE_DIR=str(work / 'uv-cache'))
subprocess.run(['uv', 'pip', 'install', '--python', sys.executable, '--no-deps', '--no-index', '--target', str(installed), str(wheel)], env=env, check=True)
for source in (root / 'src/organoid_analysis').rglob('*'):
    if source.is_file() and '__pycache__' not in source.parts and (source.suffix == '.py' or 'static' in source.parts):
        target = installed / 'organoid_analysis' / source.relative_to(root / 'src/organoid_analysis')
        assert target.read_bytes() == source.read_bytes(), str(source)
print('PASS: installed package Python/static bytes match final source', flush=True)
env['PYTHONPATH'] = str(installed)
probe = """
import organoid_analysis, numpy as np
assert 'installed/organoid_analysis' in organoid_analysis.__file__
from organoid_analysis.visualization.volume_viewer.viewer_payload import build_viewer_payload, render_viewer_html
html=render_viewer_html(build_viewer_payload([np.arange(120,dtype=np.uint16).reshape(4,5,6)],(1,1,2)))
assert len(html)>100000 and 'vtk' in html.lower()
from organoid_analysis.quantification.features import geometry
assert geometry(np.ones((3,4,5),bool),(2,3,4))[0]['volume_um3']==1440
print('PASS: isolated artifact import, offline viewer and analytical voxel volume')
"""
subprocess.run([sys.executable, '-c', probe], cwd=work, env=env, check=True)
for args in [
    ['--version'],
    ['demo', '--out', str(work / 'demo')],
    ['analyze', '--manifest', str(evidence / 'real/manifest.csv'), '--out', str(work / 'real'), '--keep-going'],
]:
    print('COMMAND: python -m organoid_analysis ' + ' '.join(args), flush=True)
    subprocess.run([sys.executable, '-m', 'organoid_analysis', *args], cwd=work, env=env, check=True)
paths = {'work': str(work), 'wheel': str(wheel), 'installed': str(installed), 'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest()}
(evidence / 'artifact-paths.json').write_text(json.dumps(paths, indent=2))
for folder in ('demo', 'real'):
    for source in (work / folder).rglob('*'):
        if source.is_file() and source.suffix in ('.json', '.csv', '.yaml'):
            target = evidence / 'final-artifact' / folder / source.relative_to(work / folder)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
print('PASS: final artifact demo and real-file execution; accuracy is a separate claim', flush=True)
