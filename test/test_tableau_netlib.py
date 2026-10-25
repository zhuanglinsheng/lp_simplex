"""Replay the bundled Bland instability offline, without generated source files."""
import csv
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]
driver = Path(sys.argv[1]).resolve()
with (root / 'data/netlib/feasible_gurobi_1e-8.csv').open() as stream:
    reference = next(row[3] for row in csv.reader(stream)
                     if row[0] == 'blend' and row[1] == 'OPTIMAL')
with tempfile.TemporaryDirectory(prefix='tableau-blend-') as temporary:
    model = Path(temporary) / 'blend.mps'
    model.write_bytes(gzip.decompress(
        (root / 'data/netlib/raw/feasible/blend.mps.gz').read_bytes()))
    for pricing in ('bland', 'dantzig'):
        for presolve in (True, False):
            for budget in (1, 10000):
                command = [str(driver), '--algorithm', 'tableau', '--criteria',
                           pricing, '--iterations', str(budget), '--expected', reference]
                if not presolve:
                    command.append('--no-presolve')
                result = subprocess.run(command + [str(model)], text=True,
                                        capture_output=True, timeout=30)
                work = re.search(r'^iterations:\s+(\d+)', result.stdout, re.M)
                status = re.search(r'^result:.*code=(\d+)', result.stdout, re.M)
                assert work and status, (command, result.stdout, result.stderr)
                if budget == 1:
                    assert result.returncode != 0 and int(status[1]) == 3 and int(work[1]) == 1, (command, result.stdout, result.stderr)
                else:
                    assert result.returncode == 0 and int(status[1]) == 0 and 0 < int(work[1]) <= budget and '[PASS]' in result.stdout, (command, result.stdout, result.stderr)
print('8 Tableau pricing/presolve/budget scenarios passed')
