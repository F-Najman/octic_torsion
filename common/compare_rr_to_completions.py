#!/usr/bin/env python3
"""Compare exact differential-space kernels with completion restrictions."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys


def matrices(source, pattern):
    out = {}
    for m in re.finditer(pattern + r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);', source, re.S):
        groups = m.groups()
        key = tuple(map(int, groups[:-3]))
        nr, nc = map(int, groups[-3:-1])
        entries = ast.literal_eval(groups[-1])
        assert len(entries) == nr * nc and set(entries) <= {0, 1}
        assert key not in out
        out[key] = ([sum(entries[i * nc + j] << j for j in range(nc)) for i in range(nr)], nc)
    return out


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            k = row.bit_length()
            if k not in pivots:
                pivots[k] = row
                break
            row ^= pivots[k]
    return len(pivots)


def apply(row, matrix):
    result = 0
    while row:
        bit = row & -row
        result ^= matrix[bit.bit_length() - 1]
        row ^= bit
    return result


def main():
    rr, completion = map(Path, sys.argv[1:3])
    old, new = rr.read_text(), completion.read_text()
    for name in ('CartierMatrix', 'ImageBasis'):
        assert matrices(old, r'\b' + name) == matrices(new, r'\b' + name)
    degrees = ast.literal_eval(re.search(r'\bds:=(\[.*?\]);', new)[1])
    assert degrees == ast.literal_eval(re.search(r'\bds:=(\[.*?\]);', old)[1])
    kernels = matrices(old, r'kernelspaces\[<(\d+),(\d+)>\]')
    restrictions = matrices(new, r'basisjets\[(\d+)\]')
    assert kernels
    checked = 0
    for (i, m), (kernel, dimension) in kernels.items():
        full, width = restrictions[i,]
        length = m * degrees[i - 1]
        assert length <= width and len(full) == dimension
        cut = [row & ((1 << length) - 1) for row in full]
        assert all(apply(v, cut) == 0 for v in kernel)
        assert rank(kernel) + rank(cut) == dimension
        checked += 1
    result = {'differential_space_reference': str(rr), 'completion_export': str(completion),
              'reference_sha256': hashlib.sha256(rr.read_bytes()).hexdigest(),
              'completion_sha256': hashlib.sha256(completion.read_bytes()).hexdigest(),
              'identical_Cartier_matrix_and_image': True,
              'exact_kernel_equalities_checked': checked,
              'boundary': 'Only exported RR kernels are compared; this is not the exhaustive divisor-rank replay.'}
    print(json.dumps(result, indent=2))
    print('RR_COMPLETION_CROSS_CHECK_DONE')


if __name__ == '__main__':
    main()
