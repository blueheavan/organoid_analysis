"""Capture source, artifact and evidence identity without self-referential hashes."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

root = Path.cwd()
out = root / 'docs/evidence/2026-09-10'
def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

patch = subprocess.check_output(['git', 'diff', '--binary', 'HEAD'])
(out / 'full-tracked.patch').write_bytes(patch)
if Path('/tmp/organoid-audit-initial.patch').exists():
    shutil.copy2('/tmp/organoid-audit-initial.patch', out / 'initial-tracked.patch')
status = subprocess.check_output(['git', 'status', '--porcelain=v1'], text=True)
paths = [root / 'pyproject.toml', root / 'pixi.lock']
for folder in ('src', 'tests', 'scripts', 'configs'):
    paths.extend(p for p in (root / folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and (p.suffix in ('.py', '.toml', '.yaml', '.yml', '.sh') or 'static' in p.parts))
source = {str(p.relative_to(root)): digest(p) for p in sorted(set(paths))}
excluded = {'docs/evidence/2026-09-10/final-snapshot.json', 'docs/evidence/2026-09-10/full-tracked.patch', 'docs/evidence/2026-09-10/initial-tracked.patch'}
untracked = {}
for name in subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard'], text=True).splitlines():
    if name in excluded or name.startswith('docs/evidence/2026-09-09/'):
        continue
    p = root / name
    if p.is_file():
        untracked[name] = digest(p)
artifact = json.loads((out / 'artifact-paths.json').read_text())
assert digest(artifact['wheel']) == artifact['wheel_sha256']
# Confirm no source drift since the final wheel was built.
installed = Path(artifact['installed']) / 'organoid_analysis'
for p in (root / 'src/organoid_analysis').rglob('*'):
    if p.is_file() and '__pycache__' not in p.parts and (p.suffix == '.py' or 'static' in p.parts):
        assert digest(p) == digest(installed / p.relative_to(root / 'src/organoid_analysis')), str(p)
weight = Path.home() / '.cellpose/models/cpdino-vitb'
record = {
    'captured_utc': datetime.now(timezone.utc).isoformat(),
    'HEAD': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'dirty': bool(status.strip()), 'git_status': status.splitlines(),
    'tracked_patch_sha256': hashlib.sha256(patch).hexdigest(),
    'initial_tracked_patch_sha256': digest(out / 'initial-tracked.patch'),
    'initial_untracked_snapshot_limitation': 'Original untracked contents not independently hashed before test extensions; initial status and fresh baseline tests retained.',
    'source_files_sha256': source,
    'source_manifest_sha256': hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest(),
    'relevant_untracked_files_sha256': untracked,
    'hash_exclusions': sorted(excluded) + ['historical docs/evidence/2026-09-09/*'],
    'artifact': artifact,
    'environment': {'python': sys.version, 'platform': platform.platform(), 'machine': platform.machine(), 'lock_sha256': digest(root / 'pixi.lock')},
    'model': {'name': 'cpdino-vitb', 'path': str(weight), 'sha256': digest(weight) if weight.is_file() else 'unavailable', 'accuracy_validation': 'INSUFFICIENT EVIDENCE'},
    'QA': {'full_pytest': '307 passed, 8 skipped, 8 deselected', 'host_native_render': '8 passed', 'lint': 'PASS', 'mypy': 'FAIL: 193 errors / 37 files', 'format': 'FAIL: 69 files', 'notebooks': 'NOT APPLICABLE', 'artifact': 'PASS execution and source identity', 'biology': 'INSUFFICIENT EVIDENCE', 'verdict': 'NOT READY FOR THE SPECIFIED RESEARCH USE'},
    'separation': 'Tier D; LIMITED INDEPENDENCE',
}
(out / 'final-snapshot.json').write_text(json.dumps(record, indent=2))
print(json.dumps({k: record[k] for k in ['HEAD', 'tracked_patch_sha256', 'source_manifest_sha256']}, indent=2))
print('Final installed artifact/source identity verified; snapshot written.')
