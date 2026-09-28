"""Failure-path checks for staged regeneration and cache invalidation, using temporary data."""
from pathlib import Path
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np

import paperB_q4
import regenerate_chi4 as pipeline


def main():
    canonical = pipeline.ROOT / 'data/chi4_zeros_pari10000_clean.txt'
    original_hash = pipeline.sha256(canonical)
    source = pipeline.ROOT / 'data/chi4_zeros_to60.txt'
    reference = pipeline.ROOT / 'data/chi4_zeros_mp_snapshot.txt'
    real_run = subprocess.run
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for stage in ('gp', 'validation', 'theorem'):
            output = root / stage

            def runner(command, **kwargs):
                if command[0] == 'fake-gp':
                    kwargs['stdout'].write(source.read_text(encoding='utf-8'))
                    if stage == 'gp':
                        raise subprocess.CalledProcessError(1, command)
                    return subprocess.CompletedProcess(command, 0)
                if stage == 'theorem' and Path(command[1]).name == 'theoremA_q4.py':
                    # Even a nonempty partial result cannot create a success manifest.
                    Path(command[-1]).write_text('partial theorem output', encoding='utf-8')
                    raise subprocess.CalledProcessError(1, command)
                return real_run(command, **kwargs)

            refs = [root / 'absent-reference.txt'] if stage == 'validation' else [reference]
            with patch('regenerate_chi4.subprocess.run', side_effect=runner):
                try:
                    pipeline.run_pipeline(output, refs, gp='fake-gp')
                except subprocess.CalledProcessError:
                    pass
                else:
                    raise AssertionError(f'{stage} failure accepted')
            assert (output / 'FAILED.json').exists()
            assert not (output / 'COMPLETE.json').exists()
            if stage != 'theorem':
                assert not (output / 'theoremA_q4.txt').exists()

        existing = root / 'existing'
        existing.mkdir()
        sentinel = existing / 'COMPLETE.json'
        sentinel.write_text('previous result', encoding='utf-8')
        try:
            pipeline.run_pipeline(existing, [reference], reuse=source)
        except FileExistsError:
            pass
        else:
            raise AssertionError('existing output directory accepted')
        assert sentinel.read_text(encoding='utf-8') == 'previous result'

        data = root / 'cache'
        data.mkdir()
        zeros = data / canonical.name
        zeros.write_bytes(source.read_bytes())
        g, pair = paperB_q4.pair_sums_chi4(data)
        with patch('paperB_q4.eri', side_effect=AssertionError('valid cache unnecessarily recomputed')):
            cached_g, cached_pair = paperB_q4.pair_sums_chi4(data)
        assert np.array_equal(cached_g, g) and np.array_equal(cached_pair, pair)
        cache = data / 'eri_list_chi4.csv'
        cache.write_text(cache.read_text(encoding='utf-8').replace(repr(float(pair[0])), '999999'), encoding='utf-8')
        _, restored = paperB_q4.pair_sums_chi4(data)
        assert np.array_equal(restored, pair), 'tampered cache reused'
        zeros.write_text('6.0209499046978\n10.2437703041669\n', encoding='utf-8')
        changed_g, changed_pair = paperB_q4.pair_sums_chi4(data)
        assert len(changed_g) == 2 and changed_g[0] != g[0]
        assert changed_pair[0] != pair[0], 'changed source reused old values'
    assert pipeline.sha256(canonical) == original_hash, 'canonical zero data modified'
    print('Regeneration checks passed: GP/validation/theorem failures, stale outputs, cache provenance.')


if __name__ == '__main__':
    main()
