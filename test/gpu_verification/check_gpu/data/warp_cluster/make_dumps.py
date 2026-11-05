#!/usr/bin/env python3
"""The CPU-dumped clusters that test_warp_cluster_cuda compares the warp kernel against (A40, gpu_options.warp_max_dim).

  python3 make_dumps.py BIN [OUT]                    the fixture files below
  python3 make_dumps.py --step0 CLUSTERS.BIN [OUT]   p1_production_cr60_gcce2_single_hahn.bin

BIN is a CPU build (bin/main.out, any preset). Each run, one rank, writes the first 8 clusters of every dimension the warp
kernel has (CCEX_WARP_CLUSTER_DUMP, CCEX_WARP_CLUSTER_DUMP_N=24, small_cluster_batch.h), 40 time points (more than the 32
lanes of a warp, so that a lane takes two points). The inputs:
  gcce_natab_*  test/code_verification/check_input_forms/inputs/Diamond_NV_natab_gcce (13C bath, NV qubit): dimension 6, 12
                and, with order 3, 24;
  gcce_p1n14_*  .../Diamond_NV_P1_10ppm_gcce with its P1 nitrogen made 14N (spin 1): dimension 18;
  gcce_hbn_*    example/tutorial/advanced/06_hBN/gCCE2_Hahn_states (until A66 07_hBN_VB_natab/single-sample_gCCE; 10B, 11B, 14N; the tutorial's quantity dm): 9, 12, 21, 27;
  cce_natab_*   .../Diamond_NV_natab_cce with order 4: bath dimension 8 and 16.
CLUSTERS.BIN is the A40 step-0 dump of the production P1 CR60 sigma4 gCCE2 input (rank 0, 2 states, dimension 12, Hahn,
151 points; idm4-gpu job 698, (internal notes)): its first 64 clusters are rewritten
in the format of the dump hook (single state, coherence).
"""
import json, os, shutil, struct, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '../../../../..'))
FORMS = os.path.join(ROOT, 'test/code_verification/check_input_forms/inputs')
HBN = os.path.join(ROOT, 'example/tutorial/advanced/06_hBN')

# name: (input folder, input file, changes)
NATAB_G = ('Diamond_NV_natab_gcce', 'ccein_nstate1.json')
P1_G = ('Diamond_NV_P1_10ppm_gcce', 'ccein_nstate1.json')
NATAB_C = ('Diamond_NV_natab_cce', 'ccein_nstate1.json')
CASES = {
    'gcce_natab_single_hahn':        (NATAB_G, dict(order=3, pulse={'count': 1})),
    'gcce_natab_single_cpmg2_dm':    (NATAB_G, dict(order=3, pulse={'count': 2}, quantity='dm')),
    'gcce_natab_single_ramsey_exact': (NATAB_G, dict(order=2, pulse={'count': 0}, evolution_update='exact')),
    'gcce_natab_ensemble_hahn':      (NATAB_G, dict(order=3, pulse={'count': 1}, nstate=0)),
    'gcce_natab_ensemble_cpmg3_dm':  (NATAB_G, dict(order=2, pulse={'count': 3}, nstate=0, quantity='dm')),
    'gcce_p1n14_single_hahn':        (P1_G, dict(order=1, pulse={'count': 1}, p1n14=True)),
    'gcce_p1n14_ensemble_cpmg2':     (P1_G, dict(order=1, pulse={'count': 2}, nstate=0, p1n14=True)),
    'gcce_hbn_single_hahn_dm':       ('hbn', dict(order=2, rbath=30, rdip=8, nstate=1, seed=20261012)),
    'gcce_hbn_ensemble_hahn':        ('hbn', dict(order=2, rbath=30, rdip=8, nstate=0, quantity='coherence')),
    'cce_natab_single_hahn':         (NATAB_C, dict(order=4, pulse={'count': 1})),
    'cce_natab_single_cpmg3':        (NATAB_C, dict(order=4, pulse={'count': 3})),
    'cce_natab_ensemble_hahn':       (NATAB_C, dict(order=4, pulse={'count': 1}, nstate=0)),
    'cce_natab_ensemble_cpmg2_exact': (NATAB_C, dict(order=4, pulse={'count': 2}, nstate=0, evolution_update='exact')),
    'cce_natab_ensemble_ramsey':     (NATAB_C, dict(order=4, pulse={'count': 0}, nstate=0)),
}


def lenient(text):
    import re
    text = re.sub(r'(#|//)[^\n]*', '', text)
    return json.loads(re.sub(r',\s*([}\]])', r'\1', text))


def fixture(binary, out):
    work = tempfile.mkdtemp(prefix='warp_cluster_')
    for name, (where, change) in CASES.items():
        change = dict(change)
        p1n14 = change.pop('p1n14', False)
        if where == 'hbn':
            # the tutorial's relative paths: ../bath/... from gCCE2_Hahn_states/
            run = os.path.join(work, name, 'gCCE2_Hahn_states')
            os.makedirs(run)
            shutil.copytree(os.path.join(HBN, 'bath/flat'), os.path.join(work, name, 'bath/flat'))
            shutil.copy(os.path.join(HBN, 'bath/hBN_natab_gyro'), os.path.join(work, name, 'bath'))
            cfg = lenient(open(os.path.join(HBN, 'gCCE2_Hahn_states/ccein.jsonc')).read())
            cfg['configuration']['bathfile'] = ['../bath/flat/bath_hBN_VB_flat_1']   # -I on the command line until A63
        else:
            folder, file = where
            run = os.path.join(work, name)
            shutil.copytree(os.path.join(FORMS, folder, 'inputfiles'), os.path.join(run, 'inputfiles'))
            cfg = json.load(open(os.path.join(FORMS, folder, 'inputfiles', file)))
        cfg.update(change)
        for key in ('nstate', 'seed'):   # the configuration block's (A64)
            if key in change: cfg['configuration'][key] = cfg.pop(key)
        n = cfg['configuration'].pop('nstate', None)   # configuration.s_range [1, N] (A94)
        if n: cfg['configuration']['s_range'] = [1, n]
        if 'evolution_update' in change:   # the execution block's (A65)
            cfg.setdefault('execution', {})['evolution_update'] = cfg.pop('evolution_update')
        if p1n14:
            d = cfg['Defect']['species'][0]
            d['types'], d['spins'], d['gyros'] = ['14N'], [1.0], [1.9337792]
            for h in d['hyperfine']: h[1] = '14N'
        cfg['nstep'] = 40
        cfg['output'] = dict(cfg.get('output', {}), file='./out')   # output.file (A65)
        json.dump(cfg, open(os.path.join(run, 'in.json'), 'w'), indent=1)
        target = os.path.join(out, name + '.bin')
        if os.path.exists(target): os.remove(target)
        env = dict(os.environ, CCEX_WARP_CLUSTER_DUMP=target, CCEX_WARP_CLUSTER_DUMP_N='8')
        subprocess.run(['mpirun', '-n', '1', binary, '-f', 'in.json'], cwd=run, env=env, check=True,
                       stdout=subprocess.DEVNULL)
        print(name, os.path.getsize(target), 'bytes')
    shutil.rmtree(work)


def step0(clusters, out, n=64):
    # step-0 record ((internal notes) readRecord): dim, qdim, npulse, nstep,
    # recursive; deltat; fraction[npulse+1]; sequence index[npulse+1]; U[npulse] (qdim x qdim); psia, psib; H, psi0; nstep CCEX values
    f = open(clusters, 'rb').read(); o = 0; taken = 0
    with open(os.path.join(out, 'p1_production_cr60_gcce2_single_hahn.bin'), 'wb') as w:
        while o < len(f) and taken < n:
            dim, qd, npulse, nstep, rec = struct.unpack_from('5i', f, o); o += 20
            dt, = struct.unpack_from('d', f, o); o += 8
            frac = struct.unpack_from('%dd' % (npulse + 1), f, o); o += 8 * (npulse + 1)
            o += 4 * (npulse + 1)
            pulses = f[o:o + 16 * npulse * qd * qd]; o += 16 * npulse * qd * qd
            psia = f[o:o + 16 * qd]; o += 16 * qd
            psib = f[o:o + 16 * qd]; o += 16 * qd
            body = f[o:o + 16 * (dim * dim + dim)]; o += 16 * (dim * dim + dim)
            cpu = f[o:o + 16 * nstep]; o += 16 * nstep
            distinct = sorted(set(frac), key=frac.index)
            same = [distinct.index(x) for x in frac]
            # mode gCCE, dim, qdim, single state, coherence, intervals, distinct, nstep; deltat; fractions; same;
            # pulses; psia; psib; qubit psi0 (not used for a single state: zeros)
            w.write(struct.pack('8i', 1, dim, qd, 0, 0, npulse + 1, len(distinct), nstep) + struct.pack('d', dt))
            w.write(struct.pack('%dd' % len(distinct), *distinct) + struct.pack('%di' % (npulse + 1), *same))
            w.write(pulses + psia + psib + bytes(16 * qd) + body + cpu)
            taken += 1
    print('p1_production_cr60_gcce2_single_hahn', taken)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--step0':
        step0(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else HERE)
    elif len(sys.argv) > 1:
        fixture(os.path.abspath(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else HERE)
    else:
        sys.exit(__doc__)
