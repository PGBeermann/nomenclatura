# Estructuras arbitrarias (como las que podría dibujar un estudiante) con grupos funcionales:
# round-trip con OPSIN (3 estilos + clase funcional) e invariancia ante el orden de los átomos.
import random, subprocess, sys
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
from validar_org import validar, can
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
GR = ['O', 'O', 'N', 'NC', 'N(C)C', '=O', '=O', 'C=O', 'C(=O)O', 'C(=O)OC', 'C(=O)OCC', 'C(N)=O', 'C(=O)NC', 'C#N', 'OC', 'OCC', 'Cl', 'Br', 'F', 'I', '[N+](=O)[O-]', 'c1ccccc1']
def aleatoria():
    rw = Chem.RWMol(); ring = random.random() < 0.5; k = 0
    if ring:
        k = random.randint(3, 8); arom = k == 6 and random.random() < 0.5
        for i in range(k): a = Chem.Atom(6); rw.AddAtom(a)
        for i in range(k): rw.AddBond(i, (i + 1) % k, Chem.BondType.SINGLE)
        if arom:
            for i in range(0, 6, 2): rw.GetBondBetweenAtoms(i, i + 1).SetBondType(Chem.BondType.DOUBLE)
    else:
        rw.AddAtom(Chem.Atom(6))
    n = random.randint(0, 12)
    for _ in range(n):
        cand = [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 6]
        a = random.choice(cand); b = rw.AddAtom(Chem.Atom(6)); rw.AddBond(a, b, Chem.BondType.SINGLE)
        try: Chem.SanitizeMol(rw.GetMol())
        except Exception: rw.RemoveAtom(b)
    smi = Chem.MolToSmiles(rw.GetMol())
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    # insaturaciones
    for _ in range(random.choice([0, 0, 1, 2])):
        bs = [b for b in m.GetBonds() if b.GetBondType() == Chem.BondType.SINGLE and not b.GetIsAromatic()
              and b.GetBeginAtom().GetTotalNumHs() > 0 and b.GetEndAtom().GetTotalNumHs() > 0]
        if not bs: break
        b = random.choice(bs); t = random.choice([Chem.BondType.DOUBLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
        old = b.GetBondType(); b.SetBondType(t)
        try: Chem.SanitizeMol(m)
        except Exception: b.SetBondType(old)
    s = Chem.MolToSmiles(m)
    # grupos
    for _ in range(random.choice([0, 1, 1, 2, 2, 3])):
        mm = Chem.MolFromSmiles(s)
        cand = [a.GetIdx() for a in mm.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() >= 1]
        if not cand: break
        g = random.choice(GR); a = random.choice(cand)
        if g == '=O' and mm.GetAtomWithIdx(a).GetTotalNumHs() < 2: continue
        rwm = Chem.RWMol(mm)
        frag = Chem.MolFromSmiles(g[1:] if g.startswith('=') else g)
        off = rwm.GetNumAtoms()
        combo = Chem.RWMol(Chem.CombineMols(rwm, frag))
        combo.AddBond(a, off, Chem.BondType.DOUBLE if g.startswith('=') else Chem.BondType.SINGLE)
        try:
            mol = combo.GetMol(); Chem.SanitizeMol(mol); s = Chem.MolToSmiles(mol)
        except Exception: pass
    return s
N = int(sys.argv[1]); smis = list({aleatoria() for _ in range(N)})
res = subprocess.run(['php', 'nombrar_org.php'], input='\n'.join(smis) + '\n', capture_output=True, text=True).stdout.rstrip('\n').split('\n')
filas = []; errs = {}
for s, r in zip(smis, res):
    f = r.split('\t')
    if f[0] == 'ERR': errs[f[1]] = errs.get(f[1], 0) + 1; continue
    if len(f) < 5: print('SALIDA RARA', s, r); continue
    filas.append([s, f[0], f[2], f[3], f[4], f[1]])
n, bad = validar(filas, [1, 2, 3, 4], {1: 'pin', 2: 'trad', 3: 'sis79', 4: 'clase'})
print('estructuras', len(smis), 'nombradas', len(filas), 'nombres', n, 'discrepancias', bad)
for e, c in sorted(errs.items(), key=lambda x: -x[1]): print('  fuera de dominio', c, e)
# invariancia
inp = []; exp = []
for f in filas[:600]:
    m = Chem.MolFromSmiles(f[0])
    for _ in range(3): inp.append(Chem.MolToSmiles(m, doRandom=True, canonical=False, kekuleSmiles=random.random() < 0.5) if True else ''); exp.append((f[1], f[5]))
out = subprocess.run(['php', 'nombrar_org.php'], input='\n'.join(inp) + '\n', capture_output=True, text=True).stdout.rstrip('\n').split('\n')
inv = 0
for e, g, s in zip(exp, out, inp):
    gg = g.split('\t')
    if (gg[0], gg[1] if len(gg) > 1 else '') != e: inv += 1; print('NO-INVARIANTE', s, e, gg[:2])
print('ordenaciones aleatorias', len(inp), 'no invariantes', inv)
