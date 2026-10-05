"""Fresh verification of this pass and retained prior artifacts."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
commands = [[sys.executable, '-m', 'unittest', 'discover', '-s', str(HERE), '-p', 'test_uncertainty.py'],
            [sys.executable, str(HERE/'run_analysis.py')],
            [sys.executable, str(HERE/'constrained_analysis.py')]]
for command, filename in zip(commands, ['tests.log', 'analysis.log', 'constrained.log']):
    p = subprocess.run(command, capture_output=True, text=True)
    (HERE/filename).write_text(p.stdout+p.stderr, encoding='utf-8')
    assert p.returncode == 0, p.stdout+p.stderr
result = json.loads((HERE/'results.json').read_text())
constrained = json.loads((HERE/'constrained_results.json').read_text())
assert abs(result['scenarios'][-1]['frozen_rule'] - 46/51) < 1e-14
assert abs(result['conditional_envelopes']['frozen_rule']['min']['frozen_rule'] - 31/60) < 1e-14
assert abs(constrained['envelopes']['frozen_rule']['min']['frozen_rule'] - 9/14) < 1e-14
checks = json.loads((HERE/'checks.json').read_text())
checks.update({'unit_tests_passed': 5, 'report_headline_numbers_checked': True,
               'constrained_assignments': constrained['assignments']})
(HERE/'verification.json').write_text(json.dumps(checks, indent=2)+'\n')
files = [p for p in HERE.iterdir() if p.is_file() and p.name != 'artifact_hashes.json']
(HERE/'artifact_hashes.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}, indent=2)+'\n')
print(json.dumps(checks, indent=2))
