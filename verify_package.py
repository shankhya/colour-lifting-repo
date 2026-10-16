"""Read-only SHA-256 verification of this distribution; Python stdlib only."""
from pathlib import Path
import hashlib

root = Path(__file__).resolve().parent
count = 0
for line in (root / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
    expected, name = line.split('  ', 1)
    target = (root / name).resolve()
    if not target.is_relative_to(root) or not target.is_file():
        raise RuntimeError('Missing or unsafe distribution path: ' + name)
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError('Checksum mismatch: ' + name)
    count += 1
print('Distribution integrity passed:', count, 'files.')
