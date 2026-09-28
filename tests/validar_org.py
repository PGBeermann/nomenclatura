# Validación nombre→estructura con OPSIN y RDKit para el motor IupacOrganica.
import subprocess, sys, os
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
JAR = os.environ.get('OPSIN_JAR', 'opsin_cli.jar')
def opsin(names):
    out = subprocess.run(['java', '-jar', JAR, '-osmi'], input='\n'.join(names) + '\n', capture_output=True, text=True).stdout.split('\n')
    out = [o for o in out if not o.startswith(('Picked', 'Run the'))]
    return out[:len(names)]
def can(s):
    m = Chem.MolFromSmiles(s) if s else None
    return Chem.MolToSmiles(m) if m else None
def validar(filas, cols, etiquetas):
    names = [f[c] for c in cols for f in filas]
    idx = [(i, c) for c in cols for i in range(len(filas))]
    out = opsin(names)
    bad = 0; n = 0
    for (i, c), nm, got in zip(idx, names, out):
        if not nm: continue
        n += 1
        if can(got) != can(filas[i][0]):
            bad += 1; print('MISMATCH', etiquetas[c], filas[i][0], '|', nm, '|', got)
    return n, bad
if __name__ == '__main__':
    src = sys.argv[1]
    if src == 'gen':
        N = sys.argv[2]; seed = sys.argv[3] if len(sys.argv) > 3 else '1'
        p = subprocess.run(['php', 'volcar_org.php', N, seed], capture_output=True, text=True)
        filas = [l.split('\t') for l in p.stdout.strip().split('\n') if l]
        print(p.stderr.strip().split('\n')[-1])
    else:
        smis = [l.strip() for l in open(src) if l.strip()]
        p = subprocess.run(['php', 'nombrar_org.php'], input='\n'.join(smis) + '\n', capture_output=True, text=True)
        filas = []
        for s, r in zip(smis, p.stdout.rstrip('\n').split('\n')):
            f = r.split('\t')
            if f[0] == 'ERR': print('ERR', s, f[1]); continue
            filas.append([s, f[0], f[2], f[3], f[4], f[1]])
    n, bad = validar(filas, [1, 2, 3, 4], {1: 'pin', 2: 'trad', 3: 'sis79', 4: 'clase'})
    print('estructuras', len(filas), 'nombres', n, 'discrepancias', bad)
