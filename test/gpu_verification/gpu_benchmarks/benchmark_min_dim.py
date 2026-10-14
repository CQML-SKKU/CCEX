#!/usr/bin/env python3
"""Where backend=auto should switch to the GPU (gpu_options.min_dim, A3 of (internal notes)).

The P1 inputs of check_input_forms (NV + P1 10 ppm, 15N), cut to the S bath sites nearest the NV with
order S-1, so that every cluster size appears with C(S,k) clusters (with order S the expansion is exact
and only the one full cluster is computed): cCCE on 7 sites gives dimensions 4, 16, ..., 4096 (7, 21,
35, 35, 21, 7 clusters), gCCE on 6 sites 12, 48, ..., 3072; the variants use one site less. Each condition runs
with backend=cpu and with backend=cuda (forced), profile on; the [profile] records give, per dimension,
the time of the build, set-up and evolution summed over ranks and the number of clusters. Printed per
condition: the time per cluster on each backend and their ratio by dimension, the smallest dimension
from which the GPU is faster at every larger one, and the largest |CPU - CUDA| of the coherence. A
condition whose CPU and CUDA results do not match gets no crossover, and the run then exits 1.
usage: benchmark_min_dim.py --cpu BIN --cuda BIN --mpirun MPIRUN --output DIR [--only NAME ...]"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'check_gpu'))
from benchmark import INPUTS, REPO, compare, parse_profiles, FILE_KEYS, get, put

# name, method folder, nstate file, sites, npulse, nstep, ranks
CONDITIONS = [
    ('cce_single',   'cce',  1, 7, 1, 40, 1),
    ('cce_ensemble', 'cce',  0, 6, 1, 40, 1),
    ('cce_pulse4',   'cce',  1, 6, 4, 40, 1),
    ('cce_steps400', 'cce',  1, 6, 1, 400, 1),
    ('cce_ranks8',   'cce',  1, 6, 1, 40, 8),
    ('gcce_single',  'gcce', 1, 6, 1, 40, 1),
    ('gcce_ensemble', 'gcce', 0, 5, 1, 40, 1),
    ('gcce_pulse4',  'gcce', 1, 5, 4, 40, 1),
    ('gcce_steps400', 'gcce', 1, 5, 1, 400, 1),
    ('gcce_ranks8',  'gcce', 1, 5, 1, 40, 8),
]


def prepare(method, nstate, sites, npulse, nstep, directory):
    """the input, cut to the nearest sites, with every path absolute"""
    source = INPUTS/f'Diamond_NV_P1_10ppm_{method}'/'inputfiles'/f'ccein_nstate{nstate}.json'
    config = json.loads(source.read_text())
    for key in FILE_KEYS:
        if get(config, key): config = put(config, {key: str((source.parent.parent/get(config, key)).resolve())})
    config['bathfile'] = [str((source.parent.parent/p).resolve()) for p in config['bathfile']]
    origin = ([float(v) for v in Path(config['Qubit']['file']).read_text().split()[:3]] if get(config, 'Qubit.file')
              else config['Qubit']['qubit'][0]['xyz'])
    rows = Path(config['bathfile'][0]).read_text().splitlines()[1:]
    near = sorted(sorted(enumerate(rows), key=lambda r: math.dist([float(v) for v in r[1].split()[:3]], origin))[:sites])
    directory.mkdir(parents=True, exist_ok=False)
    (directory/'bath').write_text(f'{sites+1}.0 0 0 0\n'+'\n'.join(r for _, r in near)+'\n')
    radius = max(math.dist([float(v) for v in r.split()[:3]], origin) for _, r in near)+1
    config = put(config, {'bathfile': [str(directory/'bath')], 'order': sites-1, 'rbath': radius, 'rdip': 2*radius,
                          'pulse.count': npulse, 'nstep': nstep, 'profile': True, 'savemode': 'coherence',
                          'max_clusters_per_order': None})
    return config


# The profile backends a run asked for with `backend` may report: backend=cpu also computes small cCCE clusters in the
# CPU batch (A30, profile backend=cpu-batch); any other backend only itself (review 20261013-180049 G3).
ACCEPTED = {'cpu': {'cpu', 'cpu-batch'}}


def check_backends(directory, backend, seen):
    other = set(seen) - ACCEPTED.get(backend, {backend})
    if other:
        raise SystemExit(f'{directory}: asked for {backend}, profiled on {", ".join(sorted(other))}')


def run(binary, config, directory, backend, ranks, mpirun, timeout):
    directory.mkdir(parents=True, exist_ok=False)
    config = dict(config, backend=backend, outfile=str(directory/'coherence'))
    (directory/'input.json').write_text(json.dumps(config, indent=2)+'\n')
    command = [mpirun, '-n', str(ranks), str(binary), '-f', str(directory/'input.json')]
    start = time.perf_counter()
    with (directory/'run.log').open('w') as log:
        subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=True)
    wall = time.perf_counter()-start
    by_dim = {}
    for p in parse_profiles((directory/'run.log').read_text()):
        d = by_dim.setdefault(int(p['dim']), {'clusters': 0, 'seconds': 0.0, 'backends': set()})
        d['clusters'] += int(p['clusters'])
        d['seconds'] += float(p['build_s'])+float(p['setup_s'])+float(p['evolution_s'])
        d['backends'].add(p['backend'])
    # a run that profiled nothing, or ran on another backend than asked, is no measurement (review 7eed38 R3)
    if not by_dim:
        raise SystemExit(f'{directory}: no [profile] records in run.log')
    check_backends(directory, backend, {b for v in by_dim.values() for b in v['backends']})
    return wall, {dim: dict(v, backends=sorted(v['backends'])) for dim, v in sorted(by_dim.items())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cpu', type=Path, required=True); ap.add_argument('--cuda', type=Path, required=True)
    ap.add_argument('--mpirun', required=True); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--only', nargs='+', choices=[c[0] for c in CONDITIONS])   # a misspelt name stops here (review 1bd0ff R5)
    ap.add_argument('--timeout', type=float, default=7200)
    args = ap.parse_args()
    # the binaries as the caller named them: the runs start in the repository (review fa7f39 R2)
    args.cpu, args.cuda = args.cpu.resolve(), args.cuda.resolve()
    for b in (args.cpu, args.cuda):
        if not b.is_file():
            sys.exit(f'no binary at {b}')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    summary = []; invalid = []
    for name, method, nstate, sites, npulse, nstep, ranks in CONDITIONS:
        if args.only and name not in args.only: continue
        config = prepare(method, nstate, sites, npulse, nstep, out/name/'input')
        row = dict(name=name, method=method, nstate=nstate, sites=sites, npulse=npulse, nstep=nstep, ranks=ranks,
                   cpu_binary=str(args.cpu), cuda_binary=str(args.cuda))
        for backend, binary in (('cpu', args.cpu), ('cuda', args.cuda)):
            row[backend+'_wall_s'], row[backend] = run(binary, config, out/name/backend, backend, ranks, args.mpirun, args.timeout)
        try:
            row['difference'] = compare(out/name/'cpu', out/name/'cuda', 1e-8, 2e-6)
        except ValueError as e:     # a CUDA result that does not match gives no crossover (review R3)
            row['difference'] = {'error': str(e)}; invalid.append(name)
        # both backends have to have done the same work: the same dimensions with the same cluster counts
        counts = {b: {dim: v['clusters'] for dim, v in row[b].items()} for b in ('cpu', 'cuda')}
        if counts['cpu'] != counts['cuda']:
            row['difference'] = dict(row.get('difference') or {}, work=f'cpu {counts["cpu"]} against cuda {counts["cuda"]}')
            if name not in invalid: invalid.append(name)
        ratio = {dim: row['cpu'][dim]['seconds']/row['cuda'][dim]['seconds']
                 for dim in row['cpu'] if dim in row['cuda'] and row['cuda'][dim]['seconds'] > 0}
        row['speedup_by_dim'] = ratio
        dims = sorted(ratio)
        row['gpu_faster_from'] = None if name in invalid else \
            next((d for i, d in enumerate(dims) if all(ratio[e] > 1 for e in dims[i:])), None)
        summary.append(row)
        print(f'{name}: wall cpu {row["cpu_wall_s"]:.1f} s, cuda {row["cuda_wall_s"]:.1f} s; GPU faster from dim '
              f'{row["gpu_faster_from"]}; difference {row["difference"]}')
        for dim in dims:
            c, g = row['cpu'][dim], row['cuda'][dim]
            print(f'    dim {dim:5d}  clusters {c["clusters"]:3d}  per cluster cpu {c["seconds"]/c["clusters"]:.3e} s'
                  f'  cuda {g["seconds"]/g["clusters"]:.3e} s ({",".join(g["backends"])})  cpu/cuda {ratio[dim]:.2f}')
        (out/'summary.json').write_text(json.dumps(summary, indent=1, default=str)+'\n')
    if invalid:
        sys.exit(f'CPU and CUDA results differ beyond the tolerance in: {", ".join(invalid)}')


if __name__ == '__main__':
    main()
