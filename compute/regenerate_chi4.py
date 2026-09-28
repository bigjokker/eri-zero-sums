"""Generate, validate and evaluate chi_-4 zeros without replacing archived inputs.

Only COMPLETE.json identifies a successful run. Every stage propagates failure;
partial outputs remain in their new directory for diagnosis.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFS = [ROOT / 'data/chi4_zeros_mp_snapshot.txt',
                ROOT / 'data/chi4_zeros_mp_6k8k.txt',
                ROOT / 'data/chi4_zeros_to60.txt']


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_stage(command, output: Path, errors: Path, stdin=None) -> None:
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    with output.open('w', encoding='utf-8') as stdout, errors.open('w', encoding='utf-8') as stderr:
        subprocess.run([str(c) for c in command], cwd=ROOT, env=env,
                       stdin=stdin, stdout=stdout, stderr=stderr, check=True)


def run_pipeline(output: Path, refs, gp=None, reuse=None, tmax=10000,
                 precision=38, subdivision=32) -> Path:
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    zero_file = output / 'zeros.txt'
    commands = []
    try:
        if reuse is not None:
            shutil.copyfile(reuse, zero_file)
        else:
            if not gp:
                raise ValueError('PARI/GP not found; pass --gp PATH or use --reuse-zeros LIST')
            script = output / 'generate.gp'
            script.write_text(
                'default(parisize, "4G");\n'
                f'default(realprecision, {precision});\n'
                'L = lfuncreate(-4);\n'
                f'z = lfunzeros(L, {tmax}, {subdivision});\n'
                'for(i=1,#z, print(z[i]));\nquit;\n', encoding='utf-8')
            commands.append([str(gp), '-q'])
            with script.open('r', encoding='utf-8') as stdin:
                run_stage(commands[-1], zero_file, output / 'gp.err', stdin)
        commands.append([sys.executable, ROOT / 'compute/validate_zeros.py', zero_file, *refs])
        run_stage(commands[-1], output / 'validation.txt', output / 'validation.err')
        theorem_file = output / 'theoremA_q4.txt'
        commands.append([sys.executable, ROOT / 'compute/theoremA_q4.py', zero_file,
                         '--output', theorem_file])
        run_stage(commands[-1], output / 'theorem.log', output / 'theorem.err')
        if not theorem_file.is_file() or not theorem_file.stat().st_size:
            raise RuntimeError('theorem generator did not produce a nonempty output')
        inputs = [Path(r).resolve() for r in refs]
        scripts = [ROOT / 'compute' / name for name in
                   ('regenerate_chi4.py', 'validate_zeros.py', 'theoremA_q4.py', 'mobius_ei.py')]
        manifest = dict(completed_utc=datetime.now(timezone.utc).isoformat(),
                        python=sys.version, commands=[[str(c) for c in cmd] for cmd in commands],
                        mode='reuse' if reuse is not None else 'PARI/GP',
                        generation=None if reuse is not None else
                        dict(tmax=tmax, precision=precision, subdivision=subdivision),
                        reused_source=str(Path(reuse).resolve()) if reuse is not None else None,
                        sha256={str(p): sha256(p) for p in
                                [zero_file, theorem_file, *inputs, *scripts]},
                        qualification='Numerical count and overlap checks; not a rigorous completeness certificate.')
        (output / 'COMPLETE.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    except Exception as exc:
        (output / 'FAILED.json').write_text(json.dumps(dict(error=str(exc)), indent=2) + '\n', encoding='utf-8')
        raise
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gp', default=shutil.which('gp'), help='PARI/GP executable path')
    parser.add_argument('--reuse-zeros', type=Path, help='Validate/evaluate an existing list without GP')
    parser.add_argument('--reference', type=Path, action='append', help='Overlap reference (repeatable)')
    parser.add_argument('--output-dir', type=Path, help='New directory; existing paths are rejected')
    parser.add_argument('--tmax', type=int, default=10000)
    parser.add_argument('--precision', type=int, default=38)
    parser.add_argument('--subdivision', type=int, default=32)
    args = parser.parse_args()
    if args.tmax < 2 or args.precision < 20 or args.subdivision < 1:
        parser.error('require tmax >= 2, precision >= 20 and subdivision >= 1')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = args.output_dir or ROOT / 'data/regenerated' / (stamp + '-' + uuid.uuid4().hex[:8])
    try:
        run_pipeline(output, args.reference or DEFAULT_REFS, args.gp, args.reuse_zeros,
                     args.tmax, args.precision, args.subdivision)
    except Exception as exc:
        print(f'FAILED: {exc}\nInspect {output}', file=sys.stderr)
        return 1
    print(f'Completed: {output.resolve()}\nReview COMPLETE.json and validation.txt before selecting these outputs.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
