#!/usr/bin/env python3
"""
Standalone xTB + g-xTB(COSMO) + CREST workflow.

Only obabel, xtb and crest are required. RDKit is optional: when available
(for example in the psi4 environment) it is used to write final.sdf;
otherwise a built-in V2000 SDF writer preserves the input topology.

Input: a conformer ensemble in SDF format (for example Ex91_confs.sdf).
The topology template is the first molecule of the ensemble SDF (RDKit reads
the first record), or a separate single-molecule SDF supplied with -template.

Steps:
  1. obabel splits the SDF into individual CONF_XXXXX.xyz files.
  2. GFN2-xTB geometry optimization in water (ALPB water).
  3. g-xTB single point energy in water (ddCOSMO):
         xtb CONF.xyz -c Q -u U --sp --gxtb --cosmo water
  4. CREST cregen clustering with rthr = 0.5 A (default).
  5. Final SDF rebuilt from the topology template, with
     Energy_xTB / Rel_Energy_gxTB / Boltzmann_Weight properties.

Before starting, the script checks that RDKit, obabel, xtb and crest exist.

Flat output directory:

  <workdir>/
  ├── 00_README.txt
  ├── 0_input/      input SDF + topology template
  ├── 1_xtbopt/     CONF_00001.xyz, .xtbopt.xyz, .xtb.log,
  │                 ensemble.xtbopt.xyz
  ├── 2_gxtb/       CONF_00001.gxtb.log, CONF_00001_gxtb.xyz,
  │                 ensemble_gxtb.xyz, energy_summary.csv
  ├── 3_cluster/    crest.log, crest_ensemble.xyz
  ├── 4_final/      final.sdf, final.csv
  └── logs/

Usage (no Schrödinger environment needed):

  source /public/gkxiao/bin/xtb_env.sh
  python3 xtb_gxtb_cluster.py Ex91_confs.sdf -charge 0 -nproc 36 -rthr 0.5
"""

import argparse
import concurrent.futures
import csv
import glob
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False

HARTREE_TO_KCAL = 627.5094740631
KCAL_PER_MOL_K = 0.001987204259
TEMPERATURE = 298.15
JOB_TIMEOUT = 6 * 3600

XTB_METHODS = {
    'gfn2': ['--gfn', '2'],
    'gfn1': ['--gfn', '1'],
    'gfnff': ['--gfnff'],
}
NO_SOLVENT = ('none', 'no', 'off', 'gas', 'gasphase', '')

RE_TOTAL_ENERGY_BOX = re.compile(
    r'\|\s*TOTAL\s+ENERGY\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh',
    re.IGNORECASE)
RE_TOTAL_ENERGY_SUMMARY = re.compile(
    r'::\s*total\s+energy\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh',
    re.IGNORECASE)
RE_TOTAL_ENERGY_LOOSE = re.compile(
    r'total\s+energy\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh', re.IGNORECASE)

logger = logging.getLogger('xtb_gxtb_cluster')

# Prefer the standalone xTB 6.7.1 build (supports --gxtb); PATH fallback.
TOOL_BIN_DIR = os.environ.get('XTB_BIN_DIR',
                              '/public/gkxiao/software/xtb-6.7.1/bin')


def findTool(name):
    env_key = name.upper().replace('-', '_') + '_EXE'
    if os.environ.get(env_key):
        return os.environ[env_key]
    candidate = os.path.join(TOOL_BIN_DIR, name)
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return shutil.which(name)


def checkEnvironment():
    """Check that RDKit, obabel, xtb and crest are available.

    :return: dict with the resolved tool paths.
    :rtype: dict
    """
    tools = {}
    missing = []
    for name in ('obabel', 'xtb', 'crest'):
        path = findTool(name)
        if path:
            tools[name] = path
        else:
            missing.append(name)
    if not HAVE_RDKIT:
        missing.append('RDKit (python -c "import rdkit")')
    if missing:
        raise RuntimeError(
            'Missing required tools: {0}.\n'
            'Activate the proper environment first, e.g.:\n'
            '  source /public/gkxiao/bin/mamba_init.sh\n'
            '  mamba activate psi4\n'
            '  source /public/gkxiao/bin/xtb_env.sh'.format(
                ', '.join(missing)))
    return tools


def readXyz(xyz_file):
    with open(xyz_file) as fh:
        lines = fh.read().splitlines()
    natoms = int(lines[0].split()[0])
    elements = []
    coords = []
    for line in lines[2:2 + natoms]:
        fields = line.split()
        elements.append(fields[0])
        coords.append((float(fields[1]), float(fields[2]), float(fields[3])))
    if len(coords) != natoms:
        raise ValueError('Truncated XYZ: {0}'.format(xyz_file))
    return elements, coords


def readMultiXyz(xyz_file):
    with open(xyz_file) as fh:
        lines = fh.read().splitlines()
    conformers = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        natoms = int(lines[i].split()[0])
        comment = lines[i + 1]
        elements = []
        coords = []
        for line in lines[i + 2:i + 2 + natoms]:
            fields = line.split()
            elements.append(fields[0])
            coords.append((float(fields[1]), float(fields[2]),
                           float(fields[3])))
        energy = None
        for token in comment.replace('=', ' ').split():
            try:
                energy = float(token)
                break
            except ValueError:
                continue
        conformers.append({'elements': elements, 'coords': coords,
                           'energy': energy})
        i += natoms + 2
    return conformers


def writeCoordsXyz(path, elements, coords, comment):
    with open(path, 'w') as fh:
        fh.write('{0}\n'.format(len(elements)))
        fh.write(comment + '\n')
        for elem, (x, y, z) in zip(elements, coords):
            fh.write('{0:<3s} {1:14.8f} {2:14.8f} {3:14.8f}\n'.format(
                elem, x, y, z))


class SdfTemplate:
    """Minimal V2000 SDF template reader/writer (topology kept as-is)."""

    def __init__(self, path):
        with open(path) as fh:
            lines = fh.read().splitlines()

        end = lines.index('$$$$') if '$$$$' in lines else len(lines)
        block = lines[:end]
        if len(block) < 4:
            raise ValueError('Not a valid SDF: {0}'.format(path))

        self.title = block[0]
        self.program = block[1]
        self.comment = block[2]
        self.counts = block[3]
        try:
            self.natoms = int(self.counts[0:3])
            self.nbonds = int(self.counts[3:6])
        except ValueError:
            raise ValueError('Only V2000 SDF is supported: {0}'.format(path))

        self.atom_lines = block[4:4 + self.natoms]
        self.bond_lines = block[4 + self.natoms:4 + self.natoms + self.nbonds]
        if len(self.atom_lines) != self.natoms:
            raise ValueError('Truncated atom block: {0}'.format(path))

        self.elements = []
        for line in self.atom_lines:
            symbol = line[31:34].strip() if len(line) > 31 else ''
            if not symbol:
                symbol = line.split()[3]
            self.elements.append(symbol)

    def write(self, fh, coords, properties, title=None):
        if len(coords) != self.natoms:
            raise ValueError('Atom count mismatch when writing SDF.')

        fh.write((title or self.title) + '\n')
        fh.write('  xtb_gxtb_cluster          3D\n')
        fh.write(self.comment + '\n')
        fh.write(self.counts + '\n')

        for line, (x, y, z) in zip(self.atom_lines, coords):
            fh.write('{0:10.4f}{1:10.4f}{2:10.4f}{3}\n'.format(
                x, y, z, line[30:]))

        for line in self.bond_lines:
            fh.write(line + '\n')

        fh.write('M  END\n')
        for key, value in properties:
            fh.write('>  <{0}>\n{1}\n\n'.format(key, value))
        fh.write('$$$$\n')


class XtbGxtbClusterWorkflow:
    def __init__(self, infile, charge, mult=1, xtb_method='gfn2',
                 xtb_opt_level='tight', xtb_solvent='water',
                 gxtb_solvent='water', nproc=8, rthr=0.5, ewin=10.0,
                 template=None, workdir=None):
        self.infile = os.path.abspath(infile)
        self.charge = int(charge)
        self.mult = int(mult)
        self.xtb_method = xtb_method
        self.xtb_opt_level = xtb_opt_level
        self.xtb_solvent = xtb_solvent
        self.gxtb_solvent = gxtb_solvent
        self.nproc = max(1, int(nproc))
        self.rthr = float(rthr)
        self.ewin = float(ewin)
        self.template_arg = os.path.abspath(template) if template else None
        self.base = os.path.splitext(os.path.basename(infile))[0]
        self.workdir = os.path.abspath(workdir or
                                       (self.base + '_xtb_gxtb'))
        self.template = None
        self.records = []

    # ------------------------------------------------------------------ setup
    def run(self):
        self.tools = checkEnvironment()
        self.prepareWorkdir()
        self.setupLogging()
        logger.info('Tools   : obabel=%s | xtb=%s | crest=%s | RDKit=yes',
                    self.tools['obabel'], self.tools['xtb'],
                    self.tools['crest'])
        logger.info('Input   : %s', self.infile)
        logger.info('Workdir : %s', self.workdir)
        logger.info('xTB     : %s / %s / ALPB %s', self.xtb_method,
                    self.xtb_opt_level, self.xtb_solvent)
        logger.info('g-xTB   : xtb --sp --gxtb --cosmo %s',
                    self.gxtb_solvent or 'none')
        logger.info('Cluster : CREST rthr=%.3f A, ewin=%.2f kcal/mol',
                    self.rthr, self.ewin)

        self.prepareInputs()
        self.splitSdf()
        self.runXtb()
        self.runGxtb()
        self.runCrest()
        self.writeFinal()
        self.writeReadme()
        logger.info('Done: %s',
                    os.path.join(self.workdir, '4_final', 'final.sdf'))

    def prepareWorkdir(self):
        if os.path.isdir(self.workdir) and os.listdir(self.workdir):
            raise RuntimeError(
                'Work directory exists and is not empty: {0}\n'
                'Use another -workdir or remove it first.'.format(
                    self.workdir))
        for sub in ('0_input', '1_xtbopt', '2_gxtb', '3_cluster', '4_final',
                    'logs'):
            os.makedirs(os.path.join(self.workdir, sub), exist_ok=True)

    def setupLogging(self):
        logger.setLevel(logging.INFO)
        logger.handlers = []
        fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        stream = logging.StreamHandler()
        stream.setFormatter(fmt)
        fh = logging.FileHandler(
            os.path.join(self.workdir, 'logs', 'workflow.log'))
        fh.setFormatter(fmt)
        logger.addHandler(stream)
        logger.addHandler(fh)

    def prepareInputs(self):
        target = os.path.join(self.workdir, '0_input',
                              os.path.basename(self.infile))
        shutil.copy(self.infile, target)

        # Topology template: use the user-provided -template if given,
        # otherwise the first molecule of the ensemble itself.
        if self.template_arg:
            self.template_path = os.path.join(self.workdir, '0_input',
                                              'topo.sdf')
            shutil.copy(self.template_arg, self.template_path)
            source = 'user template'
        else:
            self.template_path = target
            source = 'first molecule of the ensemble'
        self.template = SdfTemplate(self.template_path)
        logger.info('Topology template: %s (%d atoms, %s)',
                    self.template_path, self.template.natoms, source)
        logger.info('Charge %d, multiplicity %d', self.charge, self.mult)

    def splitSdf(self):
        obabel = self.tools['obabel']
        opt_dir = os.path.join(self.workdir, '1_xtbopt')
        pattern = os.path.join(opt_dir, 'CONF_.xyz')
        proc = subprocess.run(
            [obabel, '-isdf', self.infile, '-oxyz', '-O', pattern, '-m'],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError('obabel failed: {0}'.format(proc.stderr))

        files = glob.glob(os.path.join(opt_dir, 'CONF_*.xyz'))
        if not files:
            raise RuntimeError('obabel produced no conformers.')

        def number(path):
            match = re.search(r'CONF_(\d+)\.xyz$', path)
            return int(match.group(1)) if match else 0

        files.sort(key=number)
        for idx, old in enumerate(files, start=1):
            new = os.path.join(opt_dir, 'CONF_{0:05d}.xyz'.format(idx))
            if old != new:
                os.replace(old, new)
            elements, _ = readXyz(new)
            if len(elements) != self.template.natoms:
                raise RuntimeError(
                    'Atom count mismatch: {0} has {1} atoms, template has '
                    '{2}.'.format(new, len(elements), self.template.natoms))
            self.records.append({'index': idx, 'xyz': new,
                                 'elements': elements, 'status': 'pending'})
        logger.info('Split %d conformers into 1_xtbopt/', len(self.records))

    def jobLayout(self, njobs):
        workers = max(1, min(njobs, self.nproc // 2))
        threads = max(1, self.nproc // workers)
        return workers, threads

    def runParallel(self, tasks, workers, threads, stage):
        total = len(tasks)
        done = [0]

        def execute(task):
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = str(threads)
            env['MKL_NUM_THREADS'] = str(threads)
            try:
                with open(task['log'], 'w') as fh:
                    proc = subprocess.run(task['cmd'], cwd=task['cwd'],
                                          stdout=fh,
                                          stderr=subprocess.STDOUT,
                                          env=env, timeout=JOB_TIMEOUT)
                return proc.returncode
            except subprocess.TimeoutExpired:
                return -9
            except Exception as err:  # noqa: BLE001
                logger.error('Failed to start job: %s', err)
                return -1

        with concurrent.futures.ThreadPoolExecutor(
                max_workers=workers) as pool:
            futures = {pool.submit(execute, t): t for t in tasks}
            for fut in concurrent.futures.as_completed(futures):
                task = futures[fut]
                rc = fut.result()
                done[0] += 1
                idx = task['record']['index']
                if rc == 0:
                    logger.info('[%d/%d] CONF_%05d %s done', done[0], total,
                                idx, stage)
                else:
                    logger.warning('[%d/%d] CONF_%05d %s rc=%s', done[0],
                                   total, idx, stage, rc)

    # ------------------------------------------------------------------- xTB
    def validateUhf(self):
        """Convert multiplicity 2S+1 to the unpaired-electron count used by
        xtb (-u) and CREST (-uhf): uhf = mult - 1, and check that the electron
        count is consistent with it."""
        uhf = self.mult - 1
        table = Chem.GetPeriodicTable()
        zsum = sum(table.GetAtomicNumber(e) for e in self.template.elements)
        nelec = zsum - self.charge
        if nelec < 0:
            raise RuntimeError(
                'Charge {0} is too large: Z_sum={1}, N_elec={2}.'.format(
                    self.charge, zsum, nelec))
        if nelec % 2 != uhf % 2:
            raise RuntimeError(
                'Charge/multiplicity mismatch: {0} electrons with '
                'multiplicity {1} (uhf={2}).'.format(nelec, self.mult, uhf))
        self.uhf = uhf
        logger.info('Spin: multiplicity=%d -> uhf=%d '
                    '(Z_sum=%d, charge=%d, N_elec=%d)',
                    self.mult, uhf, zsum, self.charge, nelec)

    def runXtb(self):
        xtb = self.tools['xtb']
        self.xtb_exe = xtb
        self.validateUhf()
        logger.info('xTB executable: %s', xtb)

        workers, threads = self.jobLayout(len(self.records))
        logger.info('xTB optimization: %d workers x %d threads', workers,
                    threads)
        opt_dir = os.path.join(self.workdir, '1_xtbopt')
        tasks = []
        for rec in self.records:
            namespace = 'CONF_{0:05d}'.format(rec['index'])
            cmd = [xtb, os.path.basename(rec['xyz']), '--opt',
                   self.xtb_opt_level, '-c', str(self.charge), '-u',
                   str(self.uhf), '--namespace', namespace]
            cmd += XTB_METHODS[self.xtb_method]
            if self.xtb_solvent.lower() not in NO_SOLVENT:
                cmd += ['--alpb', self.xtb_solvent]
            tasks.append({'record': rec, 'cmd': cmd, 'cwd': opt_dir,
                          'log': os.path.join(opt_dir,
                                              namespace + '.xtb.log')})
        self.runParallel(tasks, workers, threads, 'xTB')

        good = 0
        for rec in self.records:
            namespace = 'CONF_{0:05d}'.format(rec['index'])
            opt_xyz = os.path.join(opt_dir, namespace + '.xtbopt.xyz')
            log = os.path.join(opt_dir, namespace + '.xtb.log')
            energy = self.parseXtbEnergy(log)
            if os.path.isfile(opt_xyz) and energy is not None:
                elements, coords = readXyz(opt_xyz)
                if elements != rec['elements']:
                    logger.warning('Element mismatch: %s', opt_xyz)
                    continue
                rec['xtb_energy'] = energy
                rec['coords'] = coords
                rec['status'] = 'xtb_ok'
                good += 1
            else:
                self.logTail(log)
        logger.info('xTB optimized %d/%d', good, len(self.records))
        if good == 0:
            raise RuntimeError('All xTB optimizations failed.')

        with open(os.path.join(opt_dir, 'ensemble.xtbopt.xyz'), 'w') as fh:
            for rec in self.records:
                if rec['status'] != 'xtb_ok':
                    continue
                fh.write('{0}\nenergy: {1:.10f}\n'.format(
                    len(rec['elements']), rec['xtb_energy']))
                for elem, (x, y, z) in zip(rec['elements'], rec['coords']):
                    fh.write('{0:<3s} {1:14.8f} {2:14.8f} {3:14.8f}\n'.format(
                        elem, x, y, z))

    # ----------------------------------------------------------------- g-xTB
    def checkGxtbSupport(self):
        try:
            proc = subprocess.run([self.xtb_exe, '--help'],
                                  capture_output=True, text=True, timeout=60)
            blob = (proc.stdout or '') + (proc.stderr or '')
            if 'gxtb' not in blob.lower():
                raise RuntimeError(
                    'The xtb at {0} does not support --gxtb (needs 6.7+ '
                    'with tblite; set XTB_BIN_DIR).'.format(self.xtb_exe))
        except RuntimeError:
            raise
        except Exception as err:  # noqa: BLE001
            logger.warning('Could not probe --gxtb support: %s', err)

    def runGxtb(self):
        self.checkGxtbSupport()
        targets = [r for r in self.records if r['status'] == 'xtb_ok']
        workers, threads = self.jobLayout(len(targets))
        logger.info('g-xTB single points: %d workers x %d threads', workers,
                    threads)

        gxtb_dir = os.path.join(self.workdir, '2_gxtb')
        tmp_root = tempfile.mkdtemp(prefix='gxtb_tmp_', dir=gxtb_dir)
        tasks = []
        for rec in targets:
            work = os.path.join(tmp_root, 'CONF_{0:05d}'.format(rec['index']))
            os.makedirs(work, exist_ok=True)
            shutil.copy(os.path.join(self.workdir, '1_xtbopt',
                                     'CONF_{0:05d}.xtbopt.xyz'.format(
                                         rec['index'])),
                        os.path.join(work, 'CONF.xyz'))
            cmd = [self.xtb_exe, 'CONF.xyz', '-c', str(self.charge), '-u',
                   str(self.uhf), '--sp', '--gxtb']
            if self.gxtb_solvent.lower() not in NO_SOLVENT:
                cmd += ['--cosmo', self.gxtb_solvent]
            tasks.append({'record': rec, 'cmd': cmd, 'cwd': work,
                          'log': os.path.join(work, 'gxtb.log')})
        self.runParallel(tasks, workers, threads, 'g-xTB')

        good = 0
        for task in tasks:
            rec = task['record']
            energy = self.parseGxtbEnergy(task['log'])
            if energy is None:
                self.logTail(task['log'])
                continue
            rec['gxtb_energy'] = energy
            rec['status'] = 'gxtb_ok'
            elements, coords = readXyz(os.path.join(task['cwd'], 'CONF.xyz'))
            rec['coords'] = coords
            tag = 'energy: {0:.10f} | xtb --gxtb --sp | cosmo: {1}'.format(
                energy, self.gxtb_solvent or 'none')
            writeCoordsXyz(
                os.path.join(gxtb_dir,
                             'CONF_{0:05d}_gxtb.xyz'.format(rec['index'])),
                elements, coords, tag)
            shutil.copy(task['log'],
                        os.path.join(gxtb_dir,
                                     'CONF_{0:05d}.gxtb.log'.format(
                                         rec['index'])))
            good += 1
        shutil.rmtree(tmp_root, ignore_errors=True)
        logger.info('g-xTB succeeded for %d/%d', good, len(targets))
        if good == 0:
            raise RuntimeError('All g-xTB calculations failed.')

        ok = sorted([r for r in targets if r['status'] == 'gxtb_ok'],
                    key=lambda r: r['index'])
        emin = min(r['gxtb_energy'] for r in ok)
        with open(os.path.join(gxtb_dir, 'ensemble_gxtb.xyz'), 'w') as fh:
            for rec in ok:
                fh.write('{0}\nenergy: {1:.10f}\n'.format(
                    len(rec['elements']), rec['gxtb_energy']))
                for elem, (x, y, z) in zip(rec['elements'], rec['coords']):
                    fh.write('{0:<3s} {1:14.8f} {2:14.8f} {3:14.8f}\n'.format(
                        elem, x, y, z))
        with open(os.path.join(gxtb_dir, 'energy_summary.csv'), 'w',
                  newline='') as fh:
            wr = csv.writer(fh)
            wr.writerow(['conformer', 'gxtb_energy_hartree',
                         'relative_kcal_mol'])
            for rec in ok:
                wr.writerow([
                    'CONF_{0:05d}'.format(rec['index']),
                    '{0:.12f}'.format(rec['gxtb_energy']),
                    '{0:.3f}'.format((rec['gxtb_energy'] - emin) *
                                     HARTREE_TO_KCAL)
                ])

    # ----------------------------------------------------------------- CREST
    def runCrest(self):
        crest = self.tools['crest']
        cluster_dir = os.path.join(self.workdir, '3_cluster')
        conformers = readMultiXyz(
            os.path.join(self.workdir, '2_gxtb', 'ensemble_gxtb.xyz'))
        if not conformers:
            raise RuntimeError('No conformers to cluster.')
        ref_idx = min(range(len(conformers)),
                      key=lambda i: conformers[i]['energy'])
        ref = conformers[ref_idx]
        writeCoordsXyz(os.path.join(cluster_dir, 'reference.xyz'),
                       ref['elements'], ref['coords'],
                       'energy: {0:.10f}'.format(ref['energy']))

        cmd = [crest, 'reference.xyz', '--cregen',
               os.path.join('..', '2_gxtb', 'ensemble_gxtb.xyz'),
               '-chrg', str(self.charge), '-uhf', str(self.uhf),
               '-rthr', str(self.rthr), '-ewin', str(self.ewin)]
        logger.info('CREST: %s', ' '.join(cmd))
        with open(os.path.join(cluster_dir, 'crest.log'), 'w') as fh:
            proc = subprocess.run(cmd, cwd=cluster_dir, stdout=fh,
                                  stderr=subprocess.STDOUT,
                                  env=os.environ.copy(), timeout=JOB_TIMEOUT)
        if proc.returncode != 0:
            self.logTail(os.path.join(cluster_dir, 'crest.log'))
            raise RuntimeError('CREST failed with code {0}'.format(
                proc.returncode))
        if not os.path.isfile(os.path.join(cluster_dir,
                                           'crest_ensemble.xyz')):
            raise RuntimeError('CREST did not produce crest_ensemble.xyz')
        logger.info('CREST retained %d conformers',
                    len(readMultiXyz(os.path.join(cluster_dir,
                                                  'crest_ensemble.xyz'))))

    # ----------------------------------------------------------------- final
    def writeFinal(self):
        xyz = os.path.join(self.workdir, '3_cluster', 'crest_ensemble.xyz')
        conformers = readMultiXyz(xyz)
        if not conformers:
            raise RuntimeError('Clustered ensemble is empty.')
        if any(c['energy'] is None for c in conformers):
            raise RuntimeError('Clustered conformers lack energies.')

        energies = [c['energy'] for c in conformers]
        emin = min(energies)
        rt = KCAL_PER_MOL_K * TEMPERATURE
        weights = [math.exp(-(e - emin) * HARTREE_TO_KCAL / rt)
                   for e in energies]
        total_w = sum(weights)
        order = sorted(range(len(conformers)), key=lambda i: energies[i])

        final_path = os.path.join(self.workdir, '4_final', 'final.sdf')
        if HAVE_RDKIT:
            self.writeFinalSdfRdkit(final_path, conformers, order, energies,
                                    emin, weights, total_w)
            logger.info('Final SDF written with RDKit')
        else:
            with open(final_path, 'w') as fh:
                for rank, idx in enumerate(order, start=1):
                    conf = conformers[idx]
                    if conf['elements'] != self.template.elements:
                        raise RuntimeError(
                            'Element order mismatch in clustered conformer %d'
                            % rank)
                    dE = (conf['energy'] - emin) * HARTREE_TO_KCAL
                    weight = weights[idx] / total_w
                    properties = [
                        ('Energy_xTB', '{0:.10f}'.format(conf['energy'])),
                        ('Rel_Energy_gxTB', '{0:.6f}'.format(dE)),
                        ('Boltzmann_Weight', '{0:.6f}'.format(weight)),
                    ]
                    self.template.write(
                        fh, conf['coords'], properties,
                        title='{0}_conf{1:05d}'.format(self.base, rank))
            logger.info('Final SDF written with the built-in SDF writer')

        with open(os.path.join(self.workdir, '4_final', 'final.csv'), 'w',
                  newline='') as fh:
            wr = csv.writer(fh)
            wr.writerow(['rank', 'gxtb_energy_hartree', 'relative_kcal_mol',
                         'boltzmann_weight'])
            for rank, idx in enumerate(order, start=1):
                wr.writerow([
                    rank, '{0:.10f}'.format(conformers[idx]['energy']),
                    '{0:.6f}'.format((conformers[idx]['energy'] - emin) *
                                     HARTREE_TO_KCAL),
                    '{0:.6f}'.format(weights[idx] / total_w)
                ])
        logger.info('Wrote %d conformers to %s', len(order), final_path)

    def writeFinalSdfRdkit(self, path, conformers, order, energies, emin,
                           weights, total_w):
        """Write final.sdf with RDKit, preserving the template topology."""
        template = Chem.SDMolSupplier(self.template_path, removeHs=False)[0]
        if template is None:
            raise RuntimeError('Could not read template SDF with RDKit.')
        natoms = template.GetNumAtoms()
        writer = Chem.SDWriter(path)
        try:
            for rank, idx in enumerate(order, start=1):
                conf = conformers[idx]
                if len(conf['coords']) != natoms:
                    raise RuntimeError('Atom count mismatch in cluster %d' %
                                       rank)
                mol = Chem.Mol(template)
                mol_conf = mol.GetConformer()
                for atom_idx, (x, y, z) in enumerate(conf['coords']):
                    mol_conf.SetAtomPosition(atom_idx, Point3D(x, y, z))
                dE = (conf['energy'] - emin) * HARTREE_TO_KCAL
                weight = weights[idx] / total_w
                mol.SetProp('_Name', '{0}_conf{1:05d}'.format(self.base, rank))
                mol.SetProp('Energy_xTB', '{0:.10f}'.format(conf['energy']))
                mol.SetProp('Rel_Energy_gxTB', '{0:.6f}'.format(dE))
                mol.SetProp('Boltzmann_Weight', '{0:.6f}'.format(weight))
                writer.write(mol)
        finally:
            writer.close()

    def writeReadme(self):
        with open(os.path.join(self.workdir, '00_README.txt'), 'w') as fh:
            fh.write('xTB + g-xTB(COSMO) + CREST 独立流程（无构象搜索）\n')
            fh.write('====================================================\n')
            fh.write('输入文件 : {0}\n'.format(os.path.basename(self.infile)))
            fh.write('总电荷   : {0}\n'.format(self.charge))
            fh.write('xTB 优化 : {0} / {1} / ALPB {2}\n'.format(
                self.xtb_method, self.xtb_opt_level, self.xtb_solvent))
            fh.write('g-xTB    : xtb --sp --gxtb --cosmo {0}\n'.format(
                self.gxtb_solvent or 'none'))
            fh.write('聚类     : CREST rthr={0} A, ewin={1} kcal/mol\n'.format(
                self.rthr, self.ewin))
            fh.write('完成时间 : {0}\n\n'.format(
                time.strftime('%Y-%m-%d %H:%M:%S')))
            fh.write('0_input/   输入 SDF 与拓扑模板\n')
            fh.write('1_xtbopt/  xTB 优化（平铺）\n')
            fh.write('2_gxtb/    g-xTB 单点（COSMO 水）\n')
            fh.write('3_cluster/ CREST 聚类\n')
            fh.write('4_final/   final.sdf / final.csv\n')
            fh.write('logs/      运行日志\n')

    @staticmethod
    def parseXtbEnergy(log_file):
        if not os.path.isfile(log_file):
            return None
        with open(log_file, errors='ignore') as fh:
            text = fh.read()
        matches = re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh', text)
        return float(matches[-1]) if matches else None

    @staticmethod
    def parseGxtbEnergy(log_file):
        if not os.path.isfile(log_file):
            return None
        with open(log_file, errors='ignore') as fh:
            text = fh.read()
        for pattern in (RE_TOTAL_ENERGY_BOX, RE_TOTAL_ENERGY_SUMMARY,
                        RE_TOTAL_ENERGY_LOOSE):
            matches = pattern.findall(text)
            if matches:
                return float(matches[-1])
        return None

    @staticmethod
    def logTail(log_file):
        if not os.path.isfile(log_file):
            return
        with open(log_file, errors='ignore') as fh:
            lines = fh.readlines()
        logger.warning('Last lines of %s:\n%s', log_file,
                       ''.join(lines[-10:]))


def getParser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help='Conformer ensemble SDF.')
    parser.add_argument('-charge', type=int, required=True,
                        help='Total formal charge.')
    parser.add_argument('-mult', type=int, default=1,
                        help='Spin multiplicity 2S+1 (default 1). '
                        'xtb/CREST receive uhf = mult - 1.')
    parser.add_argument('-xtb_method', default='gfn2',
                        choices=sorted(XTB_METHODS),
                        help='xTB Hamiltonian for optimization.')
    parser.add_argument('-xtb_opt_level', default='tight',
                        choices=('crude', 'loose', 'normal', 'tight',
                                 'vtight'),
                        help='xTB optimization convergence level.')
    parser.add_argument('-xtb_solvent', default='water',
                        help="ALPB solvent for xTB; 'none' for gas phase.")
    parser.add_argument('-gxtb_solvent', default='water',
                        help="ddCOSMO solvent for g-xTB; 'none' for gas.")
    parser.add_argument('-nproc', type=int,
                        default=min(8, os.cpu_count() or 1),
                        help='Total CPU cores for xTB/g-xTB stages.')
    parser.add_argument('-rthr', type=float, default=0.5,
                        help='CREST clustering RMSD threshold (default 0.5 A).')
    parser.add_argument('-ewin', type=float, default=10.0,
                        help='CREST clustering energy window (kcal/mol).')
    parser.add_argument('-template', default=None,
                        help='Optional topology SDF. Default: first molecule '
                        'of the input ensemble.')
    parser.add_argument('-workdir', default=None,
                        help='Output directory (default <base>_xtb_gxtb).')
    return parser


def main():
    parser = getParser()
    opts = parser.parse_args()
    if not os.path.isfile(opts.infile):
        parser.error('Input file not found: {0}'.format(opts.infile))
    if not opts.infile.lower().endswith(('.sdf', '.sd')):
        parser.error('This standalone workflow expects an SDF input.')
    if opts.template and not os.path.isfile(opts.template):
        parser.error('Template file not found: {0}'.format(opts.template))
    if opts.mult < 1:
        parser.error('-mult must be positive.')
    if opts.nproc < 1:
        parser.error('-nproc must be positive.')
    if opts.rthr <= 0:
        parser.error('-rthr must be positive.')

    workflow = XtbGxtbClusterWorkflow(
        infile=opts.infile, charge=opts.charge, mult=opts.mult,
        xtb_method=opts.xtb_method, xtb_opt_level=opts.xtb_opt_level,
        xtb_solvent=opts.xtb_solvent, gxtb_solvent=opts.gxtb_solvent,
        nproc=opts.nproc, rthr=opts.rthr, ewin=opts.ewin,
        template=opts.template, workdir=opts.workdir)
    try:
        workflow.run()
    except RuntimeError as err:
        print('ERROR: {0}'.format(err), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
